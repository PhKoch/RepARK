#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  RepARK.pl
#
#        USAGE:  ./RepARK.pl  
#
#  DESCRIPTION:  RepARK.pl is a wrapper script for constructing a repeat library from sequencing reads.
#
#      OPTIONS:  see below
# REQUIREMENTS:  jellyfish, velvet
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Philipp Koch, Bryan Downie     
#      COMPANY:  Leibniz Institute for Age Research - Fritz Lipmann Institute
#      VERSION:  1.2
#===============================================================================

use strict;
#use warnings;

use Getopt::Long;
use File::Spec;
use File::Basename;

my $jellyfish_path = ""; 	# example: "/bin/jellyfish-1.1.6/"
my $velvet_path = ""; 		# where your velveth/velvetg binary is located
				### can left empty, if programs are in PATH

my $RepARK_prefix = dirname(File::Spec->rel2abs($0));

# Set default values if not in config file
my $usage;
my $jellyfish_hash_size = 100000000;
my $jellyfish_kmer_size = 31;
my $thread_count = 1;
my $manual_threshold;
my $assembler = "velvet"; 
my $prefix = "RepARK_working";
my @libfiles;
push (@libfiles, "pe400.fq");
#push (@libfiles, "more.fa");		# here you can add more default read libraries
#push (@libfiles, "evenmore.fq");	# here you can add more default read libraries

my @libs;

my ($DEBUG, $NOJF);

	  
my $ret = GetOptions(	'o=s' => \$prefix,
		'l=s' => \@libs,
#		'a=s' => \$assembler,
		's=i' => \$jellyfish_hash_size,
		'k=i' => \$jellyfish_kmer_size,
		'p=i' => \$thread_count,
		't=i' => \$manual_threshold,
		'd'   => \$DEBUG,
		'nojf'=> \$NOJF,
		'h'   => \$usage,
		'help'=> \$usage);


# if no read libraries were provided via command line, use the default library 
if (!@libs){
	@libs = @libfiles;
		
}

# check for file existence
foreach  (@libs){
	unless (-e $_){
		print  "File $_ not found!\n\n"; 
		$usage = 1;
	}
}

unless ($ret) {$usage = 1;}

if ($usage)  {
	print  "Unknown option: @_\n" if ( @_ );
	print  "Usage: RepARK.pl [ options ]\n";
	print  "Options:\n"; 
	print  "  -o  output dir [RepARK_working]\n";
	print  "  -l  library file - can be used multiple times [pe400.fq]\n";
#	print  "  -a  assembler  \n";
	print  "  -s  jellyfish hash size [100000000] \n";
	print  "  -k  jellyfish kmer size [31]\n";
	print  "  -p  number of threads used by jellyfish [1]\n";
	print  "  -t  set manually a threshold, if automatic prediction is not working (this step will then be skipped)\n";
	print  "  -d  debug mode gives some more information\n";
	print  "  -n --nojf  skip Jellyfish computation (jf_RepARK.kmers|histo must exist in the working dir)\n";
	print  "  -h --help  this help\n";
	print  "\nIf no options are provided, the script looks for the demo data pe400.fq and creates a repeat library based on that.\n";
	print  "\nversion 1.2\n";
	exit;
}	  


my $cur_lib;
my @libraries;
my %config_hash;
my @command_list;


mkdir $prefix;
foreach my $lib (@libs) { 
	symlink "../$lib", "$prefix/$lib";
}


chdir $prefix;


## Jellyfish
if ($NOJF) {
	if (-e "jf_RepARK.kmers") {print "using existing kmers file \"jf_RepARK.kmers\"\n" if ($DEBUG);}
#		else {print STDERR "file \"jf_RepARK.kmers\" not found\n"; }
		
	if (-e "jf_RepARK.histo") {print "using existing histogram file \"jf_RepARK.histo\"\n" if ($DEBUG);}
#		else {print STDERR "file \"jf_RepARK.histo\" not found\n"; }
}
else {
	run_jellyfish($jellyfish_hash_size, $jellyfish_kmer_size, $thread_count, @libs);
}

# check for jf files 
unless (-e "jf_RepARK.kmers"){
	print STDERR "file \"jf_RepARK.kmers\" not found. Please remove the -n/--nojf flag or check the jellyfish command.\n";
	exit;
}
unless (-e "jf_RepARK.histo"){
	print STDERR "file \"jf_RepARK.histo\" not found. Please remove the -n/--nojf flag or check the jellyfish command.\n";
	exit;
}


## Use manually set threshold or calculate theshold
if ($manual_threshold){
	print "Theshold provided by user: $manual_threshold\n";
	open F, "jf_RepARK.kmers";
	open G, ">jf_RepARK.repeat.kmers";
	while (my $count = <F>) { 
		my $seq = <F>;
		chomp $count;
		$count =~ s/^>//; 
		if ($count >= $manual_threshold) { 
			print G ">$count\n$seq";
		}
	}
	close F;
	close G;
}
else {
	run_calc_threshold();
}


## Assembly
if ($assembler eq "velvet") { 
	my $velvet_cmd = "${velvet_path}velveth velvet_repeat_lib 29 -fasta jf_RepARK.repeat.kmers > velvet.log 2>&1 ; ${velvet_path}velvetg velvet_repeat_lib -cov_cutoff auto -exp_cov auto -scaffolding no >> velvet.log 2>&1 ";
	print  "Assembling: $velvet_cmd\n";
	system $velvet_cmd;
	symlink "velvet_repeat_lib/contigs.fa", "repeat_lib.fasta"; 
}
else {
	print "Assembler \"$assembler\" not known. Assembly step skipped.\n";
}

chdir "..";

print "Done. Check the directory $prefix for results.\n";
		

############################################################################################
#
# Subroutines from here on


sub run_jellyfish {
	my ($jf_hash_size, $jf_kmer_size, $threads, @read_set) = @_;

	my $jellyfish_command = "${jellyfish_path}jellyfish count -s $jf_hash_size -C -m $jf_kmer_size -t $thread_count -o jellyfish ";
	$jellyfish_command .= join " ", @read_set;
	print  "Counting kmers with: $jellyfish_command\n";
	system $jellyfish_command;

	if (-e "jellyfish_1") {
		my $jf_merge = "${jellyfish_path}jellyfish merge -s $jf_hash_size -o mergedjellyfish.db jellyfish_*";
		print  "Merging dbs with: $jf_merge\n";
		system ($jf_merge);
		if (-e "mergedjellyfish.db"){
			rename ("mergedjellyfish.db", "jf_RepARK.db");
			unlink glob "jellyfish_*" or warn "Could not unlink file: $!";
		}
	}
	elsif (-e "jellyfish_0") { 
		rename ("jellyfish_0", "jf_RepARK.db");
	}
	else {
		print "Jellyfish did not properly finish. Please check the command and rerun.\n";
		exit;		
	}
	
	my $jf_dump = "${jellyfish_path}jellyfish dump jf_RepARK.db > jf_RepARK.kmers";
	my $jf_histo = "${jellyfish_path}jellyfish histo jf_RepARK.db > jf_RepARK.histo";
	system $jf_dump;
	system $jf_histo;
}

# Run the loop to calculate threshold. It probably runs the subroutine &calc_theshold several times.
sub run_calc_threshold {
#	my $low;
	my $high;
	my $peak;
	my $coverage;
	for (my $cov_value = 10; $cov_value < 200; $cov_value += 20) {
		#my $calc_thresholds_cmd = "$RepARK_prefix/CalculateThresholds.pl -c $cov_value jf_RepARK.histo"; 
		#print  "Calculating k-mer frequency threshold\n" if ($DEBUG);
		my @lines = &calc_threshold("jf_RepARK.histo", $cov_value);
		foreach my $line (@lines) { 
			if ($line =~ /Peak: (\d+)/) { $peak = $1; }
		#	if ($line =~ /Left Threshold: (\d+)/) { $low = $1; }
			if ($line =~ /Right Threshold: (\d+)/) { $high = $1; }
		}
#		if (abs($cov_value - $peak) >= 10) {
#			$cov_value = $peak - 20;
#		}
#		elsif ($high) {
#			$coverage = $cov_value;
#			$cov_value = 200;
#		}

		if (!$high){
			next;
		}
		else {
			last;	
		}

	}
	unless ($high) { 
		print  "Couldn't calculate k-mer thresholds. Check your histogram file to make sure it makes sense (peak and valley)\n";
		exit;
	}
	$high *= 2;
	print "Theshold used: $high\n";
	open F, "jf_RepARK.kmers";
	open G, ">jf_RepARK.repeat.kmers";
	while (my $count = <F>) { 
		my $seq = <F>;
		chomp $count;
		$count =~ s/^>//; 
		if ($count >= $high) { 
			print G ">$count\n$seq";
		}
	}
	close F;
	close G;
}


#===============================================================================
#  DESCRIPTION: Performs a linear fit to the ascending/descending portions of
#		k-mer frequency histogram curve to predict frequencies at which
#		error/repetitive k-mers will occur. Provides lower and upper bounds.
#  
#  RETURNS:	"right threshold"
#
#  USAGE: 	calc_threshold (histogram_file,expected_coverage)
#	 	calc_threshold (histogram_file,expected_coverage,sensitivity,errorsensitivity)
#
#   		with
#		sensitivity: Sensitivity of repeat detection (default is 0.01)	
#		errorsensitivity: Sensitivity of error detection (default is 0.1)	
#
#  Note on sensitivity: 
#		* To aggressively include all "good" kmers, set sensitivity to 0.
#		* To exclude as many repetitive kmers as possible, set sensitivity to 1.
#
#===============================================================================
#sub calc_threshold ($$;$$)
sub calc_threshold {
	my $file = $_[0];
	# Coverage is necessary for smoothing of data. Smoothing window size determined later.
	my $coverage = $_[1];
	my $sensitivity = 0.01;
		if ($_[2]) {$sensitivity = $_[2];}
	my $errorsensitivity = 0.1;
		if ($_[3]) {$errorsensitivity = $_[3];}

	my $return = "";

	my @raw;
	my @x_values;
	my $kmer_count = 0;
	for (my $i = 1; $i < 10000; $i++) {
		$x_values[$i] = $i;
		$raw[$i] = 0;
	}

	open FILE, $file or die "Couldn't open file $file";
	my @line;
	while (<FILE>) { 
		@line = split /\s+/;
		$raw[$line[0]] = $line[1];
		$kmer_count += $line[1] * $line[0];
	}
	close FILE;

	$return .= "Using coverage: $coverage\n";
	print "\tUsing coverage: $coverage\n" if ($DEBUG);

	# Remove the last entry. Jellyfish puts everything > 10000 into one bin which causes a peak
	# at the end of the curve.
	pop @raw;
	pop @x_values;

	# Smoothing window is 1 for each 20x coverage above 10.
	# e.g. coverage 30 -> smoothing window 1, coverage 70 -> window 3
	my $smooth_param = int($coverage/10);
	if ($smooth_param < 3) { $smooth_param = 1; }
	elsif ($smooth_param < 5) { $smooth_param = 3; }
	elsif (($smooth_param %2 ) != 1) { $smooth_param += 1; }

	# Offset is to exclude the lowest values after taking the derivative.
	my $smooth_offset = int($smooth_param/2);

	# Valley assigned elsewhere. This is probably obsolete.
	my $valley = 0;
	my $peak = 0;
	my $best_peak_val = 0;
	my $is_ascending = 0;
	my $last_value = $raw[1];
	for (my $i = 2; $i < $#raw; $i++) {
		if ($is_ascending) { 
			if ($raw[$i] > $best_peak_val) { 
				$peak = $i;
				$best_peak_val = $raw[$i];
			}
		}
		elsif ($last_value < $raw[$i]) {
			$is_ascending = 1;
		}
		else { $last_value = $raw[$i]; }
	}


	# Do smoothing;
	my @smoothed_raw = smooth_array($smooth_param, @raw);
	my @first_deriv;
	for (my $i =  1+$smooth_offset; $i < $#smoothed_raw - 1; $i++) {
		if ($smoothed_raw[$i-1]) { 
			$first_deriv[$i] = ($smoothed_raw[$i] - $smoothed_raw[$i-1])/($x_values[$i] - $x_values[$i-1]);
		}
	}

	$last_value = $raw[$peak];
	for (my $i = $peak; $i > 0; $i--) { 
		if ($smoothed_raw[$i] > $last_value) { 
			$valley = $i + 1;
			$i = 0;
		}
		$last_value = $smoothed_raw[$i];
	}
	$return .= "Peak: $peak\nValley: $valley\n";
	print "\tPeak: $peak\n\tValley: $valley\n" if ($DEBUG);


	# Guess genome size based on the total number of k-mers divided by the peak frequency.
	# If there are excessive error k-mers, The genome size will be over-estimated.  
	my $genome_size = int($kmer_count/$peak);
	
	$return .= "Predicted genome size: $genome_size\n";
	print "\tPredicted genome size: $genome_size\n" if ($DEBUG);
	
	# Smooth again (we need the third derivative.
	my @first_deriv_smoothed = smooth_array($smooth_param,@first_deriv);

	my @second_deriv;
	for (my $i = 2+$smooth_offset; $i < $#first_deriv_smoothed - 2; $i++) {
		if ($first_deriv_smoothed[$i-1]) { 
			$second_deriv[$i] = ($first_deriv_smoothed[$i] - $first_deriv_smoothed[$i-1])/($x_values[$i] - $x_values[$i-1]);
		}
	}

	# Final smoothing
	my @second_deriv_smoothed = smooth_array($smooth_param,@second_deriv);
	my $left_regression_end;
	# TO CHANGE: 0 should be $valley here
	for (my $i = 0; $i < $peak; $i++) {
		# Left regression start/end is between the valley and peak
		if (($second_deriv_smoothed[$i] * $second_deriv_smoothed[$i-1]) < 0) { 
			$left_regression_end = $x_values[$i+1];
		}
	}

	# Smooth again so that we know when the rate of curving changes (3rd derivative)
	my @third_deriv;
	for (my $i = 3+$smooth_offset; $i < $#second_deriv_smoothed - 3; $i++) {
		if ($second_deriv_smoothed[$i-1]) { 
			$third_deriv[$i] = ($second_deriv_smoothed[$i] - $second_deriv_smoothed[$i-1])/($x_values[$i] - $x_values[$i-1]);
		}
	}

	my @third_deriv_smoothed = smooth_array($smooth_param,@third_deriv);

	# These are backup values in case the next section fails.
	# But doesn't fit the logic exactly. Will need to adjust this.
	# Should include whether or not the curve is up or down (second derivative)
	my $left_regression_start = $valley;
	$left_regression_end = $peak;
	for (my $i = $valley + 1; $i < $peak; $i++) { 
		if ($third_deriv_smoothed[$i] == 0) { $third_deriv_smoothed[$i] = $third_deriv_smoothed[$i-1]; }
		elsif (($third_deriv_smoothed[$i] * $third_deriv_smoothed[$i-1]) < 0) { 
			if ($third_deriv_smoothed[$i] < 0) {
				if (!defined($left_regression_start)) {
					# This is where the errors end
					$left_regression_start = $x_values[$i];
					if ($left_regression_start > $peak) {
						$left_regression_start = $valley;
					}
				}
				else { $left_regression_end = $i - 1; }
			}
		}
	}
	$return .= "Left regression start: $left_regression_start\nLeft regression end: $left_regression_end\n";
	print "\tLeft regression start: $left_regression_start\n\tLeft regression end: $left_regression_end\n" if ($DEBUG);
	
	# Determine window for right regression fit.
	my $right_regression_start;
	my $right_regression_end;
	for (my $i = $peak; $i < $#third_deriv_smoothed - 3; $i++) { 
		if ($third_deriv_smoothed[$i] == 0) { $third_deriv_smoothed[$i] = $third_deriv_smoothed[$i-1]; }
		elsif (!$right_regression_start) { 
			if (($third_deriv_smoothed[$i] * $third_deriv_smoothed[$i-1]) < 0) { 
				$right_regression_start = $i;
			}
		}
		else {
			if (($third_deriv_smoothed[$i] * $third_deriv_smoothed[$i-1]) < 0) { 
				$right_regression_end = $i - 1;
				if ($right_regression_end == ($right_regression_start + 1)) { $right_regression_end++; }
				last;
			}
		}
	}

	$return .= "Right regression start: $right_regression_start\nRight regression end: $right_regression_end\n";
	print "\tRight regression start: $right_regression_start\n\tRight regression end: $right_regression_end\n" if ($DEBUG);

	my $total = 0;
	my $count = 0;
	# Take the average of the slope between each point between left regression start/end
	for (my $i = $left_regression_start; $i <= $left_regression_end; $i++) { 
		$total += $raw[$i] - $raw[$i-1];
		$count++;
	}
	
	my $leftslope;
	if ($count == 0){ die "ERROR: Please check the k-mer histogram. Is there no clear peak near your expected genome coverage?";}
	else { $leftslope = $total/$count;}
	
	# Determine the intercept of the left regression.
	$total = 0;
	$count = 0;
	for (my $i = $left_regression_start; $i <= $left_regression_end; $i++) { 
		$count++;
		$total += int($raw[$i] - ($i * $leftslope));
	}
	
	my $left_intercept;
	if ($count == 0){ die "ERROR: Please check the k-mer histogram. Is there no clear peak near your expected genome coverage?";}
	else { $left_intercept = int($total/$count);}


	# Take the average of the slope between each point between right regression start/end
	$total = 0;
	$count = 0;
	my $tot;
	for (my $i = $right_regression_start; $i <= $right_regression_end; $i++) { 
		$count++;
		$total += $raw[$i] - $raw[$i-1];
		$tot = $raw[$i] - $raw[$i-1];
	}

	my $rightslope;
	if ($count == 0){ die "ERROR: Please check the k-mer histogram. Is there no clear peak near your expected genome coverage?";}
	else { $rightslope = $total/$count;}


	# Determine the intercept of the right regression.
	$total = 0;
	$count = 0;
	if ($right_regression_end < ($right_regression_start + $smooth_param)) { 
		$right_regression_end = $right_regression_start + $smooth_param;
	}
	for (my $i = $right_regression_start; $i <= $right_regression_end; $i++) { 
		$count++;
		$total += int($raw[$i] - ($i * $rightslope));
	}

	my $right_intercept;
	if ($count == 0){ die "ERROR: Please check the k-mer histogram. Is there no clear peak near your expected genome coverage?";}
	else { $right_intercept = int($total/$count);}


	# Predict where errors start after adjusting for error sensitivity
	my $error_sum;
	my $extrapolated_sum;
	my $left_threshold;
	my $percent;
	for (my $i = $left_regression_start; $i > 0; $i--) { 
		$error_sum = $raw[$i];
		$extrapolated_sum = int(($leftslope * $x_values[$i]));# + $leftintercept);
		if ($extrapolated_sum < 0) { $extrapolated_sum = 0; }
		if ($error_sum == 0) { 
			$left_threshold = $x_values[$i];
			print "($left_threshold)($x_values[$i])($i)\n";
			last;
		}
		else { 
			$percent = $extrapolated_sum/$error_sum;
			if ($percent > $errorsensitivity) {
				$left_threshold = $x_values[$i];
				last;
			}
		}
	}
	
	$error_sum = 0;
	$extrapolated_sum = 0;
	$percent = 0;

	my $right_threshold = 1 + int((0- $right_intercept)/$rightslope);
	#$x0_intercept = 1 + int((0- $right_intercept)/$rightslope);

	# Predict where errors start after adjusting for sensitivity
	#if ($sensitivity > 0.9) { $sensitivity = 0.9; }
	#for ($i = $x0_intercept ; $i  >= $right_regression_start; $i--) {
	#	if ($percent > $sensitivity) {
	#		$right_threshold = $i;
	#		last;
	#	}
	#		
	#	$actual_val = $raw[$i];
	#	$extrapolated_val = int(($rightslope * $x_values[$i]) + $right_intercept);
	#	$percent = $extrapolated_val/$actual_val;
	#}
#	print "---------$right_threshold\n";
	unless ($right_threshold) { 
		print "Couldn't calculate right threshold. Try a different sensitivity value\n";
		return; 
		$right_threshold = $right_regression_start; 
		}
	unless ($left_threshold) { 
		$left_threshold = $left_regression_end; 
		print "Couldn't calculate left threshold. Try adjusting the coverage or check your histogram file!\n";
		return;
	}

	if ($right_threshold < $right_regression_end) { 
		$right_threshold = $right_regression_end;
	}
	if ($valley) { 
		$left_threshold = $valley;
	}
	#$right_threshold *=2;

	$return .= "Left Threshold: $left_threshold\n";
	$return .= "Right Threshold: $right_threshold\n";
	print "\tLeft Threshold: $left_threshold\n" if ($DEBUG);
	print "\tRight Threshold: $right_threshold\n" if ($DEBUG);

	
	return $return;
	
}


### subs for calc_threshold

# averages values in a array with its neighbors (window size determined by $number)
sub smooth_array {
	my $number = shift;
	my @array = @_;
	my @return_array;
	
	unless (($number %2) == 1) { print STDERR "Array smoothing requires odd integer.\n"; exit (0); }
	my $div = int($number/2);
	my $count;
	my $value;
	for (my $i = $div; $i < $#array; $i++) {
		$count = 0;
		$value = 0;
		next if !defined($array[$i]);
		for (my $j = $i - $div; $j <= $i + $div; $j++) {
			if ($j < 0) { $j = 0; }
			if ($j <= $#array) { 
				if ($array[$j]) { 
					$value += $array[$j];
					$count++;
				}
			}
		}
		if ($count) { 
			$return_array[$i] = int($value/$count);
		}
	}
	return @return_array;
}

# Add an array together (better as eval)
sub sum {
	my @array = @_;

	my $val = 0;
	
	foreach my $num (@array) {
		if ($num =~ /\D/) { return "error"; }
		$val += $num;
	}

	return $val;
}
