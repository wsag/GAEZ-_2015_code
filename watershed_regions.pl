#!/usr/bin/perl -w

#######################################################################
#
#	Finds watershed and river network regions of about the same area
#	defined by a given thresholds. Those primarily search for large
#	river tributaries to be used as regions.
#
#	Written by Dr. A. Prusevich (alex.proussevitch@unh.edu)
#
#	September 2015
#	Modified	November 2015	Added locations of region in/outflows
#			June 2017	Some useful modification to be used in other projects
#
#######################################################################

use strict;
use File::Basename;
use File::Path;
use Geo::GDAL;
use Geo::Proj4;
use Math::Trig qw/pi/;;
use PDL;
use PDL::IO::FlexRaw;
use PDL::NetCDF;
use PDL::NiceSlice;
use Fcntl;
use Inline qw/Pdlpp/;
use RIMS;		### WSAG UNH module

my $init_path	= '/net/nfs/zero/data3/WBM_TrANS/wbm_trans/';
my ($init_file, $io_file)	= ($init_path.'wbm_path.init', $init_path.'wbm_io.pl');
{
  local @ARGV	= ($init_file);
  require $io_file;
}

use vars qw(*OLDERR);		# To avoid silly message-
open OLDERR, ">&STDERR";	# "No UNIDATA NC_GLOBAL:Conventions attribute"
STDOUT->autoflush(1);			# Disable buffering

#######################################################################
#############   Files, Directories and other Inputs     ###############

	### Mandatory input files
my $base_dir	= '/net/nfs/zero/home/WBM_TrANS/data/watershed_regions/hyd_dir_5m_ng/';
my $networkFile	= '/net/nfs/squam/raid/data/GLEAM/Network/hyd_dir_5m_ng.asc';			# Network file
# my $base_dir	= '/net/nfs/zero/home/WBM_TrANS/data/watershed_regions/HiMAT_full_210_Subset/';
# my $networkFile= '/net/nfs/zero/home/WBM_TrANS/data/HiMAT_full_210_Subset.asc';		# Network file
my $upsTrsSrcDir= $base_dir.'upstream_transect/';						# Source dir for transects

	### Optional input files (leave blank '' if you do not want it)
my $basinIDFile = '/net/nfs/squam/raid/data/GLEAM/Network/hyd_dir_5m_ng_IDs.asc';		# Basin ID File
my $upsAreaFile	= '/net/nfs/squam/raid/data/GLEAM/Network/hyd_dir_5m_ng_upstrArea.asc';		# Upstream area file
# my $upsAreaFile= '/net/nfs/zero/home/WBM_TrANS/data/HiMAT_full_210_Subset_upstrArea.asc';	# Upstream area file
# my $basinIDFile = '/net/nfs/zero/home/WBM_TrANS/data/HiMAT_full_210_Subset_IDs.asc';		# Basin ID File
my $mouthFile	= '/net/nfs/squam/raid/data/GLEAM/Network/hyd_dir_5m_ng_Mouth_ID.csv';		# File containing mouth coords
my $attr_file	= '';	# Used for basin name only to build $mouthFile, if not given above

	### Output files
my $regionFile	= $base_dir.'hyd_dir_5m_ng_Subset_regions_0.02M-0.2M.asc';
my $regAttFile	= $base_dir.'hyd_dir_5m_ng_Subset_regions_0.02M-0.2M.csv';

	###  Executable utilityfiles
my $cellTblCode	= "/net/nfs/zero/data3/WBM_TrANS/cell_table_pp.pl -v -sub -1 $networkFile";
my $rvrMouthCode= '/net/nfs/zero/data3/WBM_TrANS/river_mouth_locations.pl';
my $upstreamCode= '/net/home/cv/alexp/perl/wbm/upstream_transect.pl';
# my $upstreamCode= '/net/nfs/zero/data3/WBM_TrANS/upstream_transect.pl';

	###  Sub-basin area thresholds to be used in this pre-processing
my @threshold	= (20000, 200000);		# Region area thresholds (Min, Max), km2
# my @threshold	= (50000, 500000);		# Region area thresholds (Min, Max), km2	690 regions
# my @threshold	= (500000, 1500000);		# Region area thresholds (Min, Max), km2	WAY TOO LARGE... No use
# my @threshold	= (10000, 100000);		# Region area thresholds (Min, Max), km2

#######################################################################
##################     Build Input Files      #########################

(my $cellTableFile = $networkFile) =~ s/\.\w+$/.csv/;
 my $run_cell_tbl  = 0;

unless (-e $basinIDFile) {
  ($basinIDFile = $networkFile) =~ s/\.\w+$/_IDs.asc/;
  $run_cell_tbl++ unless -e $basinIDFile;
}
unless (-e $upsAreaFile) {
  ($upsAreaFile = $networkFile) =~ s/\.\w+$/_upstrArea.asc/;
  $run_cell_tbl++ unless -e $upsAreaFile;
}
system $cellTblCode if $run_cell_tbl;

#######################################################################
##################     Read Input Data        #########################

my  $extent	= get_extent( $networkFile);
my  $network	= $$extent{mask};
my  $basinID	= read_raster($basinIDFile);
my  $upstrArea	= read_raster($upsAreaFile);
my  $cellArea	= cell_area(lonLat($extent), $extent);	### Cell area in km^2
my ($ma, @mouth)= mouth_attr( $networkFile, $mouthFile, $basinIDFile, $attr_file, $rvrMouthCode);

	### Initializations
my  $regions	= zeroes(long,$$extent{ncols},$$extent{nrows})->copybad($$extent{mask});;

unless (-e $upsTrsSrcDir) { mkpath($upsTrsSrcDir,0,0775) or die "Cannot create-\n$upsTrsSrcDir\n"; }
my $regID	= 0;
my %regName;

#######################################################################
##################     First Cut with Coastal Regions     #############

print "Processing initial step...\n";

my @todoList;
my @count	= (0,0,0,0);
my %seaDone;

foreach my $rec (@mouth) {
  my $id	= $$rec[$$ma{BasinID}];
  my $area	= $$rec[$$ma{BasinArea}];
  my $seaCode	= $$rec[$$ma{SeaBasinCode}];
     $seaCode  += 1000 if $$rec[$$ma{Endorheic}];	# We assume, max(SeaBasinCode) < 1000		
  my($x, $y)	=($$rec[$$ma{X}], $$rec[$$ma{Y}]);

	###############################################################
	### Coastal sea-basin regions (plus small endorheric watersheds)
  if ($area < $threshold[0]) {
    my $numb;
		### Naming of Sea-basin regions
    # if ($seaCode > 0) {
      if (exists $seaDone{$seaCode}) {
	$numb =  $seaDone{$seaCode};
      }
      else {
	$regName{++$regID}{Name}	= $seaCode == 0 ? 'Unmatched coastal small basins' :
		$seaCode == 1000 ? 'Unmatched endorheric small basins' :
		$$rec[$$ma{SeaBasinName}].', sea-basin'.($$rec[$$ma{Endorheic}]? ' (Endorheic)' : '');
	$regName{  $regID}{Numb}	= 1;
	$seaDone{$seaCode}		= $regID;
	$numb = $regID;
	$count[1]++;
      }
    # }
		### Naming of small endorheric watershed regions
    # else {
	# $regName{++$regID}{Name}	= $$rec[$$ma{BasinName}].', small endorheric watershed';
	# $regName{  $regID}{Numb}	= 1;
	# $numb = $regID;
	# $count[0]++;
    # }
    $regName{$numb}{Area}	+= $area;

    $regions	+= ($basinID == $id)->long * $numb;
    push @{$regName{$numb}{xyOut}}, [$x, $y];
  }

	###############################################################
	### Coastal regions within size range
  elsif ($area < $threshold[1]) {
    my $suffix	= $seaCode == 0 ? 'endorheric' : 'coastal';
    $regName{++$regID}{Name}	= $$rec[$$ma{BasinName}].", whole $suffix watershed";
    $regName{  $regID}{Numb}	= 1;
    $regName{  $regID}{Area}	= $area;
    $count[2]++;

    $regions	+= ($basinID == $id)->long * $regID;
    push @{$regName{$regID}{xyOut}}, [$x, $y];
  }

	###############################################################
	### Large watersheds for  processing in the next step
  else {
    push @todoList, [$x, $y];
    $count[3]++;
  }
}

print "\tCoastal         regions = $count[0]\n";
print "\tEndorheric      regions = $count[1]\n";
print "\tWhole watershed regions = $count[2]\n";
print "\tLarge watersheds to do  = $count[3]\n";

#######################################################################
##################     Processing Large Watersheds     ################
my $stepCount = 0;
my %basinReg;

while (@todoList) {
  my @todoListNew;
  printf "Processing step # %d: %d large (sub)watersheds...\n", ++$stepCount, scalar(@todoList);
  my $todoCount	= 0;

  foreach my $colRow (@todoList) {
    printf "   (Sub)watershed %d of (%d):\n", ++$todoCount, scalar(@todoList);
       @count		= (0,0,0);
    my $upstrFile	= $upsTrsSrcDir. "colRow_$$colRow[0]_$$colRow[1].csv";
    my $upstreamMask	= $network->upstreamMask(@$colRow);
    my $upArea		= $upstrArea * $upstreamMask;
    my $bsnID		= $basinID->at(@$colRow);
#     my @splitList	= ([@$colRow, $upstrArea->at(@$colRow)]);	# Node has (col,row,upstream_area)

    my $subNetwork	= $network * $upstreamMask;
#        $subNetwork(@$colRow) .= 0;
#     my $subNetworkCut	= $subNetwork->copy;

	### Build transect
    system(     "$upstreamCode -ct $cellTableFile -o $upstrFile $$colRow[0] $$colRow[1]")  unless -e $upstrFile;
    die "Failed: $upstreamCode -ct $cellTableFile -o $upstrFile $$colRow[0] $$colRow[1]\n" unless -e $upstrFile;
    my($hdr, @upstrTable) = read_table($upstrFile);
	### Transect cell indices
    my $indTrs	= pdl(map([$$_[$$hdr{X}],$$_[$$hdr{Y}]], @upstrTable));

	### Mask out main stem of the river
    $upArea->indexND($indTrs)	.= 0;
	### Search tributaries of needed size
    do {
      my @mouthInd	= unflat($network->dims, $upArea->flat->maximum_ind);
      my $area		= $upstrArea->at(@mouthInd);

	###############################################################
	### Large watersheds for  processing in the next loop
      if ($area >= $threshold[1]) {
# 	my @split	= flowto($network,@mouthInd);
	my $tributary	= $network->upstreamMask(@mouthInd);
# 	   $tributary(@split) .= 1;
	   $upArea     *=!$tributary;
	   $subNetwork *=!$tributary;

# 	push @splitList,   [@split, $upstrArea->at(@split)];
	push @todoListNew, \@mouthInd;
	$count[2]++;
      }
	###############################################################
	### Tributary regions within size range
      elsif ($area >= $threshold[0]) {
	$basinReg{ $bsnID}++;
	$regName{++$regID}{Name}	= $mouth[$bsnID-1][$$ma{BasinName}].', whole tributary';
	$regName{  $regID}{Numb}	= $basinReg{$bsnID};
	$regName{  $regID}{Area}	= $area;

# 	my @split	= flowto($network,@mouthInd);
	my $tributary	= $network->upstreamMask(@mouthInd);
# 	   $tributary(@split) .= 1;
	   $upArea     *=!$tributary;
	   $subNetwork *=!$tributary;
	   $regions    += $tributary->long * $regID;
	push @{$regName{$regID}{xyOut}}, \@mouthInd;

# 	push @splitList,   [@split, $upstrArea->at(@split)];
	$count[1]++;
      }
    } while($upArea->max > $threshold[0]);
#     @splitList = sort {$$b[2] <=> $$a[2]} @splitList;

	###############################################################
	### Regions along main stem, but outside of large tributaries

    my($trRow,$trRowPrev)	= (0, 0);
    my $nowMask	= $subNetwork->upstreamMask($upstrTable[0][$$hdr{X}], $upstrTable[0][$$hdr{Y}]);
    my $nowArea	=($cellArea * $nowMask)->sum;

    while ($nowArea > $threshold[1]) {
      my $startMask = $nowMask->copy;
      my $startArea = $nowArea;

		### Search for the region of right size
      do {		  $trRow++;
	$nowMask	= $subNetwork->upstreamMask($upstrTable[$trRow][$$hdr{X}], $upstrTable[$trRow][$$hdr{Y}]);
	$nowArea	=($cellArea * $nowMask)->sum;
      } while ($startArea-$nowArea < $threshold[0]);

      $basinReg{ $bsnID}++;
      $regName{++$regID}{Name}	= $mouth[$bsnID-1][$$ma{BasinName}].', main-stem section';
      $regName{  $regID}{Numb}	= $basinReg{$bsnID};
      $regName{  $regID}{Area}	=($cellArea * ($startMask - $nowMask))->sum;

      $regions    += ($startMask - $nowMask)->long * $regID;
      push @{$regName{$regID}{xyOut}}, [$upstrTable[$trRowPrev][$$hdr{X}], $upstrTable[$trRowPrev][$$hdr{Y}]];
      push @{$regName{$regID}{xyIn}},  [$upstrTable[$trRow    ][$$hdr{X}], $upstrTable[$trRow    ][$$hdr{Y}]];
      $trRowPrev = $trRow;
      $count[0]++;
    }
		### Last upstream region
    $basinReg{ $bsnID}++;
    $regName{++$regID}{Name}	= $mouth[$bsnID-1][$$ma{BasinName}].', last upstream section';
    $regName{  $regID}{Numb}	= $basinReg{$bsnID};
    $regName{  $regID}{Area}	= $nowArea;

    $regions    += $nowMask->long * $regID;
    push @{$regName{$regID}{xyOut}},[$upstrTable[$trRowPrev][$$hdr{X}], $upstrTable[$trRowPrev][$$hdr{Y}]];
    $count[0]++;

    print "\tMain stem       regions = $count[0]\n";
    print "\tWhole tributary regions = $count[1]\n";
    print "\tLarge tributaries to do = $count[2]\n";
  }
  @todoList = @todoListNew;
}

#######################################################################
#############     Find Downstream Region ID       #####################

foreach my $regID (keys %regName) {
  my $downID = $regions->at(flowto($network,@{$regName{$regID}{xyOut}[0]}));
  $regName{$regID}{downstrReg} = ($downID eq 'BAD' || $downID == $regID) ? '' : $downID;
}

#######################################################################
#############     Save Data and Attributes       ######################

write_gridascii($regionFile,  $regions, $extent, {FORMAT => '%d'});
write_attrib(   $regAttFile, \%regName);

#######################################################################

print "\nDone ($regID regions):\n\t$regionFile\n\t$regAttFile\n";

close OLDERR;
exit;

#######################################################################
######################  Functions  ####################################

sub write_attrib
{
  my ($file, $regName) = @_;

  open (FILE,">$file") or die "Couldn't open $file, $!";
	### Print Header
    print FILE "ID	Name	Area_km2	DownstrID	xyIn	xyOut\n";
	### Print Data
    foreach my $i (sort {$a <=> $b} keys(%$regName)) {
     (my $name = $$regName{$i}{Name}) =~ s/GHAAS/Unnamed/;
      printf FILE "%d\t%s\t%.1f\t%s\t%s\t%s\n", $i,
	sprintf("%s: Region %d ($i)", $name, $$regName{$i}{Numb}), $$regName{$i}{Area}, $$regName{$i}{downstrReg},
	join(';',map(join(':',@$_),@{$$regName{$i}{xyIn}} )),
	join(';',map(join(':',@$_),@{$$regName{$i}{xyOut}}));
    }
  close FILE;
}

#######################################################################

sub mouth_attr
{
  unless (-e $mouthFile) {
    ($mouthFile = $networkFile) =~ s/\.\w+$/_Mouth_ID.csv/;
    return read_table($mouthFile) if -e $mouthFile;

	### Build mouth attributes file
    my $attr_option = -e $attr_file ? "-att $attr_file" : '';
    system "$rvrMouthCode -b $basinIDFile $attr_option $networkFile";
  }

  return read_table($mouthFile);
}

#######################################################################

sub flowto
{
  my ($dir, $x, $y) = @_;
  my @direction	= ([0,1,1,0,-1,-1,-1,0,1], [0,0,1,1,1,0,-1,-1,-1]);
  my %log2	= (-1=>0,0=>0,1=>1,2=>2,4=>3,8=>4,16=>5,32=>6,64=>7,128=>8);

  my $dir_to	= $dir->at($x, $y);
  my $x_to	= $x + $direction[0][$log2{$dir_to}];
  my $y_to	= $y + $direction[1][$log2{$dir_to}];

  return $x_to, $y_to;
}

#######################################################################

sub unflat
{
  my ($dimX, $dimY, $ind) = @_;

  my $y = int($ind / $dimX);
  my $x =     $ind - $dimX * $y;

  return $x, $y;
}

#######################################################################

sub script_dir
{
  $0=~/^(.+[\\\/])[^\\\/]+[\\\/]*$/;
  return $1 || "./";
}

#######################################################################
###################  PDL::PP Functions  ###############################

__DATA__
__Pdlpp__

#######################################################################

pp_addhdr('
  #include <unistd.h>       /* we need defs of XXXX */
  #include <stdio.h>
  #include <math.h>

  static void upstr(PDL_Byte * flowDir, PDL_Long * stack, PDL_Long n_size, PDL_Long m_size, PDL_Long N, PDL_Long M)
  {
    long from[2][8] = { {1,1,0,-1,-1,-1,0,1} , {0,1,1,1,0,-1,-1,-1} };
    long dir, ind, xx, yy;
    long NN = N;	long count = 1;
    long MM = M;	long pos   = 0;
    stack[0] = N + M*n_size;

    while (pos < count) {
      MM = stack[pos] / n_size;
      NN = stack[pos] - n_size*MM;
      pos++;
      for (dir=0; dir<8; dir++) {
	xx  = NN - from[0][dir];
	yy  = MM - from[1][dir];
	if (xx<0 || yy<0 || xx==n_size || yy==m_size) continue;
	ind = xx + yy*n_size;
	if (flowDir[ind] == (0x01<<dir)) stack[count++] = ind;
      }
    }
  }
');

#######################################################################

pp_def('upstreamMask', HandleBad => 1,
  Pars => 'byte flowDir(n,m);
    int N(); int M();
    byte [o] mask(n,m);',
  Code => '
    int ind;
    int n_size = $SIZE(n);    int NN = $N();
    int m_size = $SIZE(m);    int MM = $M();
    int *myStack;	myStack = malloc(n_size*m_size*sizeof *myStack);

    loop(n,m) %{ $mask() = 0; %}  	//	Initialization of the output arrays
    for (ind=0; ind < n_size*m_size; ind++) {
      myStack[ind] = -1;
    }

    upstr($P(flowDir),myStack,n_size,m_size,NN,MM);
    ind = 0;
    while (myStack[ind] != -1 && ind < n_size*m_size) {
      MM = myStack[ind] / n_size;
      NN = myStack[ind] - n_size*MM;
      ind++;
      $mask(n=>NN,m=>MM) = 1;
    }
');

#######################################################################
pp_done();
