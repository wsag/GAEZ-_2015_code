#!/usr/bin/perl -w

#######################################################################
#
#	GAEZ crop data pre-processing for use by WBM.
#	CDL_crops.pl is used as a template for this script.
#
#	Written by Dr. A. Prusevich (alex.proussevitch@unh.edu)
#
#	August, 2018
#		Modified-
#	January, 2021	- Added patching of “PlantDayFile” and “SeasonLengthFile” files
#
#######################################################################

use strict;
use Benchmark;
use Cwd   qw/abs_path/;
use File::Basename;
use File::Path;
use File::Temp;
use Math::Trig qw/pi/;
use Math::VecStat;
use PDL;
use PDL::Image2D;
use PDL::IO::FlexRaw;
use PDL::Math;
use PDL::NetCDF;
use PDL::NiceSlice;
use Geo::GDAL;
use Fcntl;
use Inline qw/Pdlpp/;
use Time::JulianDay;
use Time::DaysInMonth;
use version;
use RIMS;		### WSAG UNH module
use RIMS::DataSet;	### WSAG UNH module

my ($init_file, $io_file) = ('/net/home/cv/alexp/perl/wbm/wbm_dev/wbm_path.init', '/net/home/cv/alexp/perl/wbm/wbm_dev/wbm_io.pl');
{
  local @ARGV	= ($init_file);
  require $io_file;
}			### wbm_path.init must be in the same directory
my %path	= read_init($init_file);
my $PATH	= get_file_path();
my $time_start	= time();

use vars qw(*NEWERR *OLDERR);	# To avoid silly warning messages from GDAL, such as
open NEWERR, ">/dev/null";	# "No UNIDATA NC_GLOBAL:Conventions attribute"
open OLDERR, ">&STDERR";	# "No UNIDATA NC_GLOBAL:Conventions attribute"

set_autopthread_targ(4);	# Number of CPU threads
set_autopthread_size(1);	# Piddle size lower limit (in Meg) for multi-threading

#######################################################################
###########	Input		#######################################

my $dir_root	= '/net/nfs/squam/raid/data/GLEAM';
my $dir_IN	= "$dir_root/GAEZ_harvArea_original_2015";
my $dir_in	= "$dir_root/GAEZ_harvArea_2015";		# Folder with pre-split CropsNES and Othercereals
my $dir_out	= "$dir_root/GAEZ_crop_frac_2015";
my $dir_PD_in	= '/net/nfs/bog/raid/data/dgrogan/MIRCA_data/WBM_format/Start_day';
my $dir_SL_in	= '/net/nfs/bog/raid/data/dgrogan/MIRCA_data/WBM_format/Season_length';
my $dir_PD_out	= "$dir_root/GAEZ_crop_PlantDay_2015";
my $dir_SL_out	= "$dir_root/GAEZ_crop_SeasonLength_2015";
my $ID_file	= "$dir_root/GAEZ_to_MIRCA_mapping.csv";

unless (-e $dir_out) { mkpath($dir_out,0,0775) or die "Cannot create-\n$dir_out\n"; }

			### Generic dataset metadata
my %meta = ('Var_Scale' => 1, 'Var_Offset' => 0, 'Processing' => '', 'Projection' => 'epsg:4326');

#######################################################################
###########	Read/Process Input Files		###############

my @GAEZ_files	= <$dir_IN/*.tif>;
my $extent	= get_extent($GAEZ_files[0]);
my $cell_area	= cell_area(lonLat($extent), $extent);		# in km2
my $kha_frac	= 1000 / (100*$cell_area);
my $zeroes	= zeroes(double,$$extent{ncols},   $$extent{nrows})   ->copybad($$extent{mask});
my($hdr, @IDs)	= read_table($ID_file);
my %ID_GZ	= read_table_hash($ID_file,'GAEZ_crop');
my %ID_MC	= read_table_hash($ID_file,'MIRCA2000_crop');

	### Read Input DataSets
print "\n\tReading GAEZ and MIRCA data...\n";
my %GZ_data	= read_GZ(@GAEZ_files);
my %MC_data	= read_MC(%ID_MC);

	### Pre-split CropsNES and Othercereals
print "\tPre-splitting GAEZ data...\n";
pre_split('CropsNES',	 'Others_perennial','Others_annual');
pre_split('Othercereals','Rye',		    'Millet');

	### Split GAEZ data by MIRCA subgroups. And find GAEZ cropland area
print "\tSplitting GAEZ data to subgroups...\n";
my @GZ_crops	= sort keys %GZ_data;
foreach my $crop (@GZ_crops) {
  print "\tSplitting GAEZ $crop...\n";
  my $MC_crop	= $ID_GZ{$crop}{MIRCA2000_crop};
  my $MC_cropID	= $ID_GZ{$crop}{MIRCA2000_crop_ID};
  my %MC_N	=(Irrigated => $ID_GZ{$crop}{Irr_Subgroups}, Rainfed => $ID_GZ{$crop}{Rfd_Subgroups});

		### Patch MIRCA crop files for pixels where GAEZ exists
  foreach my $type ('Irrigated','Rainfed') {
    my $time	= time();
    my $TYPE	= uc($type);
    my $tp	= $type eq 'Irrigated' ? 'Irr' : 'Rf';
    my $MC_ID	= $type eq 'Irrigated' ? $MC_cropID : $MC_cropID+26;
    my @GZ_file	= map "$dir_out/GAEZ_$crop\_$type\_sub$_.tif", 1 .. $MC_N{$type};
    if (-e $GZ_file[0]) {
      printf "\t\tGAEZ %9s $crop crop data is already done\n", $type;
      next;
    }
		### Must use copies, because different GAEZ crops are using same MIRCA crop
    my %MC_CROP		= map(($_ => $MC_data{$MC_crop}{$type}{$_}->copy), 1..$MC_N{$type});
       $MC_CROP{TT}	= $MC_data{$MC_crop}{$TYPE}->copy;
		### Nearest neighbor search
    my $GZ_mask		=  $GZ_data{   $crop}{$type}	> 0;
    my $MC_mask		=  $MC_data{$MC_crop}{$TYPE}	> 0;
    my $diff_mask	= ($GZ_mask - $MC_mask)		> 0;
    my $idxMC		= whichND($MC_mask);
    my $idxMsk		= whichND($diff_mask);
    my $idx		= $MC_mask->nearest_neighbor($idxMsk);
		### Nearest neighbor patching
   map $MC_CROP{$_}->indexND($idxMsk) .= $MC_data{$MC_crop}{$type}{$_}->indexND($idx),  1..$MC_N{$type};
       $MC_CROP{TT}->indexND($idxMsk) .= $MC_data{$MC_crop}{$TYPE}    ->indexND($idx);
    foreach my $sub (1..$MC_N{$type}) {
      my $file_PD = "$dir_PD_out/GAEZ_$crop\_$type\_sub$sub.tif";
      my $file_SL = "$dir_SL_out/GAEZ_$crop\_$type\_sub$sub.tif";
      my $data_PD = read_GDAL($extent,\%meta,0,"NETCDF:$dir_PD_in/$tp"."Crop$MC_ID\_sub$sub\_start_day.nc:start_day",        1, 0);
      my $data_SL = read_GDAL($extent,\%meta,0,"NETCDF:$dir_SL_in/$tp"."Crop$MC_ID\_sub$sub\_season_length.nc:season_length",1, 0);
      my $DATA_PD = $data_PD->copy;
      my $DATA_SL = $data_SL->copy;
	 $DATA_PD->indexND($idxMsk) .= $data_PD->indexND($idx);
	 $DATA_SL->indexND($idxMsk) .= $data_SL->indexND($idx);
      write_tif($extent, 'Int16', $file_PD, $DATA_PD, {CO=>{COMPRESS=>'LZW'}});
      write_tif($extent, 'Int16', $file_SL, $DATA_SL, {CO=>{COMPRESS=>'LZW'}});
    }
    printf "%8d pixels patched for %9s $MC_crop in MIRCA to split $crop GAEZ data\n", $idx->dim(1), $type;

		### Split GAEZ data by patched MIRCA data
    foreach my $sub (1 .. $MC_N{$type}) {
      my $data = condition_slice($MC_CROP{TT},$MC_CROP{$sub}/$MC_CROP{TT} * $GZ_data{$crop}{$type}, 0)->hclip(1);
		### Save split GAEZ data
      write_tif($extent, 'Float32', $GZ_file[$sub-1], $data, {CO=>{COMPRESS=>'LZW'}});
    }
    printf "\t\tDone! - %d hours, %d minutes, and %d seconds\n", time_used($time,time());
  }
  delete  $GZ_data{$crop};	# Free some memory
}
	### Write cropland area totals
my %GZ_area = GZ_Area();
write_tif($extent, 'Float32', "$dir_out/GAEZ_Irrigated.tif", $GZ_area{Irrigated}, {CO=>{COMPRESS=>'LZW'}});
write_tif($extent, 'Float32', "$dir_out/GAEZ_Rainfed.tif",   $GZ_area{Rainfed},   {CO=>{COMPRESS=>'LZW'}});
write_tif($extent, 'Float32', "$dir_out/GAEZ_Total.tif",     $GZ_area{Total},     {CO=>{COMPRESS=>'LZW'}});

#######################################################################
printf "\n\nTime used for the job - %d hours, %d minutes, and %d seconds\n", time_used($time_start,time());
print  "\n\tAll Done!\n\n";

close OLDERR;	close NEWERR;
exit;

#######################################################################
######################  Functions  ####################################

sub read_GZ
{
  my @files = @_;
  my %data;
  map { $data{$1}{$2} = $kha_frac * read_raster($_) if m/HarvArea_([&\w]+)_(Irrigated|Rainfed).tif/ } @files;
# foreach my $crop (keys %data) { foreach my $type ('Irrigated','Rainfed') {
#   printf "%20s %9s = %d\n", $crop,$type,$data{$crop}{$type}->ngood;
# }} die;

  return %data;
}

#######################################################################

sub read_MC
{
  my %ID_MC = @_;
  my %MC_data;
  foreach my $crop (keys %ID_MC) { foreach my $type ('Irrigated','Rainfed') {
    my $TYPE	= uc($type);
    my $tp	= $type eq 'Irrigated' ? 'Irr' : 'Rf';
    my %tp	= ('Irr' => 'Irr', 'Rf' => 'Rfd');		# Bad naming or directory and files...
    my $id	= $type eq 'Irrigated' ? $ID_MC{$crop}{MIRCA2000_crop_ID} : $ID_MC{$crop}{MIRCA2000_crop_ID} + 26;
    foreach my $sub (1 .. $ID_MC{$crop}{$tp{$tp}.'_Subgroups'}) {
      $MC_data{$crop}{$type}{$sub} = read_GDAL($extent, \%meta, 0,
	'NETCDF:/net/nfs/bog/raid/data/dgrogan/MIRCA_data/WBM_format/'.$tp."Crop_fraction/$tp".
		"Crop$id\_sub$sub\_fraction.nc:crop$id\_subcrop$sub\_fraction", 1, 0);
      $MC_data{$crop}{$TYPE} += $MC_data{$crop}{$type}{$sub};
  }}}
  return %MC_data;
}

#######################################################################

sub pre_split
{
  my($crop,  $crop1, $crop2)	=  @_;
  my($crop_1,$crop_2)		=  map $crop."_$_", $crop1,$crop2;
     $crop_1			=~ s/Others_//;
     $crop_2			=~ s/Others_//;

  foreach my $type ('Irrigated','Rainfed') {
    my $file1	= "$dir_in/GAEZAct2015_HarvArea_$crop_1\_$type.tif";
    my $file2	= "$dir_in/GAEZAct2015_HarvArea_$crop_2\_$type.tif";
    if (-e $file1) {
      $GZ_data{$crop_1}{$type}	= $kha_frac * read_raster($file1);
      $GZ_data{$crop_2}{$type}	= $kha_frac * read_raster($file2);
      print "\t\tPre-split crop $crop\_$type is read\n";
      next;
    }

    my $TYPE	= uc($type);
    my $MC_sum	= ($MC_data{$crop1}{$TYPE} + $MC_data{$crop2}{$TYPE});
		### Must use copies, because different GAEZ crops are using same MIRCA crop
    my %MC_CROP		= map(($_ => $MC_data{$_}{$TYPE}->copy), $crop1,$crop2);
		### Nearest neighbor search
    my $GZ_mask		=  $GZ_data{$crop }{$type}	> 0;
    my $MC_mask		=  $MC_sum			> 0;
    my $diff_mask	= ($GZ_mask - $MC_mask)		> 0;
    my $idxMC		= whichND($MC_mask);
    my $idxMsk		= whichND($diff_mask);
    my $idx		= $MC_mask->nearest_neighbor($idxMsk);

    $MC_CROP{$crop1}->indexND($idxMsk) .= $MC_data{$crop1}{$TYPE}->indexND($idx);
    $MC_CROP{$crop2}->indexND($idxMsk) .= $MC_data{$crop2}{$TYPE}->indexND($idx);
    $MC_sum	    ->indexND($idxMsk) .= $MC_sum		 ->indexND($idx);
		### Pre-split here
    $GZ_data{$crop_1}{$type}	= condition_slice($MC_sum, $MC_CROP{$crop1}/$MC_sum * $GZ_data{$crop}{$type}, 0);
    $GZ_data{$crop_2}{$type}	= condition_slice($MC_sum, $MC_CROP{$crop2}/$MC_sum * $GZ_data{$crop}{$type}, 0);
		### Save pre-split data
    write_tif($extent, 'Float32', $file1,  $GZ_data{$crop_1}{$type}/$kha_frac,{CO=>{COMPRESS=>'LZW'}});
    write_tif($extent, 'Float32', $file2,  $GZ_data{$crop_2}{$type}/$kha_frac,{CO=>{COMPRESS=>'LZW'}});

    printf "%8d pixels patched for %9s $crop1/$crop2 in MIRCA to split $crop GAEZ data\n", $idxMsk->dim(1), $type;
  }

  delete $GZ_data{$crop};
}

#######################################################################

sub GZ_Area
{
  my %GZ_area;

  foreach my $crop (@GZ_crops) {
    my %MC_N	= (Irrigated => $ID_GZ{$crop}{Irr_Subgroups}, Rainfed => $ID_GZ{$crop}{Rfd_Subgroups});
    foreach my $type ('Irrigated','Rainfed') {
      $GZ_area{$type} += cat(map(read_raster("$dir_out/GAEZ_$crop\_$type\_sub$_.tif"),1..$MC_N{$type}))->mv(2,0)->maximum;
  } }
  $GZ_area{Total}	= $GZ_area{Irrigated} + $GZ_area{Rainfed};
  $GZ_area{Total}	= $GZ_area{Total}	->hclip(1);
  $GZ_area{Irrigated}	= $GZ_area{Irrigated}	->hclip(1);
  $GZ_area{Rainfed}	= $GZ_area{Rainfed}	->hclip(1);

  return %GZ_area;
}

#######################################################################

sub MIRCA_list
{
  my $file = '/net/nfs/zero/data3/WBM_TrANS/spreadsheets/MIRCA_complete_landCoverParameters.csv';
  my($hdr,@data) = read_table($file);
  return grep m/\w+/, map $$_[$$hdr{CropFracFile}], @data;
}

#######################################################################
###################  PDL::PP Functions  ###############################

__DATA__
__Pdlpp__

#######################################################################

pp_def('nearest_neighbor', HandleBad => 1,
  Pars => 'int mask(n,m);
    int     indexIN( l,k);
    int [o] indexOUT(l,k);',

  Code => '
    int s_max  = $SIZE(n) > $SIZE(m) ? $SIZE(n) : $SIZE(m);	// Maximum search distance
    int i, m, ii, jj, iii, jjj, dx, dy, rad, go;
    double RAD;

	//	Initialization of the output array
    loop(l,k) %{ $indexOUT() = $indexIN(); %}

	//	Search for nearest neighbor

    for(i=0; i<$SIZE(k); i++) {
      ii = $indexIN(l=>0, k=>i);
      jj = $indexIN(l=>1, k=>i);
      go = 1;
      for(rad=1; rad<s_max && go; rad++) {
	dy = 0;
	for(dx=rad; dx>=0 && go; dx--) {		// Check distance on the first quadrant only
	  iii = ii + dx;
	  for(m=dy; m<=rad && go; m++) {
	    RAD = sqrt(pow((double)m - 0.5, 2) + pow((double)dx - 0.5, 2));
	    if (RAD <= (double)rad) {
	      dy = m;
	      jjj = jj + dy;				// Upper-right quadrant
	      if (jjj >= 0 && jjj < $SIZE(m) && iii >= 0 && iii < $SIZE(n) &&
			$ISGOOD($mask(n=>iii, m=>jjj)) && $mask(n=>iii, m=>jjj) > 0) {
		$indexOUT(l=>0, k=>i)	= iii;
		$indexOUT(l=>1, k=>i)	= jjj;
		go		= 0;
		break;
	      }
	      jjj = jj - dy;				// Lower-right quadrant
	      if (jjj >= 0 && jjj < $SIZE(m) && iii >= 0 && iii < $SIZE(n) &&
			$ISGOOD($mask(n=>iii, m=>jjj)) && $mask(n=>iii, m=>jjj) > 0) {
		$indexOUT(l=>0, k=>i)	= iii;
		$indexOUT(l=>1, k=>i)	= jjj;
		go		= 0;
		break;
	      }
	      iii = ii - dx;				// Lower-left quadrant
	      if (jjj >= 0 && jjj < $SIZE(m) && iii >= 0 && iii < $SIZE(n) &&
			$ISGOOD($mask(n=>iii, m=>jjj)) && $mask(n=>iii, m=>jjj) > 0) {
		$indexOUT(l=>0, k=>i)	= iii;
		$indexOUT(l=>1, k=>i)	= jjj;
		go		= 0;
		break;
	      }
	      jjj = jj + dy;				// Upper-left quadrant
	      if (jjj >= 0 && jjj < $SIZE(m) && iii >= 0 && iii < $SIZE(n) &&
			$ISGOOD($mask(n=>iii, m=>jjj)) && $mask(n=>iii, m=>jjj) > 0) {
		$indexOUT(l=>0, k=>i)	= iii;
		$indexOUT(l=>1, k=>i)	= jjj;
		go		= 0;
		break;
	      }
	    }
	    else { break; }
	  }
	}
      }
      if ($indexOUT(l=>0, k=>i)==$indexIN(l=>0, k=>i) && $indexOUT(l=>1, k=>i)==$indexIN(l=>1, k=>i)) {
	printf( "\nNearest neighbor search in Pdlpp has FAILED. Aborting...\n\n" );
	exit(0xFF);		// It is exit with BAD status
      }
    }
');

#######################################################################
pp_done();
