#!/opt/perl/bin/perl

use strict;
use warnings;
use Carp;

# $LastChangedDate: 2013-09-16 08:33:47 -0500 (Mon, 16 Sep 2013) $
# $LastChangedRevision: 57 $
# $LastChangedBy: konganti $

# Most of these modules come with default perl installation, but I am being
# extra cautious.

check_and_load_modules(['Getopt::Long', 'Pod::Usage', 'File::Basename']);

# Declare initial global variables

my ($map, $ped, $from_bp, $bp_interval_string, $kb_interval_string,
    $mb_interval_string, $no_of_markers, $output_path, $quiet, $help,
    $upper_interval, $lower_interval, $which_metric, $chr, $max_interval,
    $plink_interval_string, $fastphase, $overlap, $stop, $filename, $all_chrs,
    $chrs, $max_coord, $max_allowed_kb_limit, $max_allowed_mb_limit,
    $user_defined_upper_interval, $user_defined_lower_interval, %map_data,
    %ped_data);

my $is_valid_option = GetOptions ('help|?' => \$help,
                                  'quiet' => \$quiet,
                                  'stop' => \$stop,
                                  'all_chrs' => \$all_chrs,
                                  'ped=s' => \$ped,
                                  'map=s' => \$map,
                                  'chr=s' => \$chr,
                                  'no_of_markers=i' => \$no_of_markers,
                                  'bp_interval=s' => \$bp_interval_string,
                                  'kb_interval=s' => \$kb_interval_string,
                                  'mb_interval=s' => \$mb_interval_string,
                                  'output_path=s' => \$output_path,
                                  'fastphase:s' => \$fastphase,
                                  'overlap=s' => \$overlap,
                                  'filname=s' => \$filename);

# Verify if option supplied is valid or not

verify_options($is_valid_option, $help);
verify_input_files([$map, $ped]);

my $checked_output_path = check_output_path($output_path);

# Create input and ouput directories for fastPHASE

my $fastphase_in_dir = $checked_output_path . 'fastPHASE_corrected_input/';
my $fastphase_out_dir = $checked_output_path . 'fastPHASE_output/';

execute_system_command(0,
                       "\nfastPHASE input files will be stored at:\n$fastphase_in_dir\n\n");
execute_system_command(0,
                       "\nMaking directory $fastphase_in_dir\n\n"),
system("mkdir $fastphase_in_dir")
    if (!-d $fastphase_in_dir);


execute_system_command(0,
                       "\nfastPHASE output files will be created at:\n$fastphase_out_dir\n\n"),
execute_system_command(0,
                       "\nMaking directory $fastphase_out_dir\n\n"),
system("mkdir $fastphase_out_dir")
    if (!-d $fastphase_out_dir && defined ($fastphase));

# Check for plink installation
my $plink_bin_path = `which plink`;
chomp $plink_bin_path;

execute_system_command(0,
                           "\nLooking for plink on your system ... \n");
execute_system_command(0, "\nFound plink at:\n$plink_bin_path\n");

confess error("Cannot find plink on your system.\nYour \${PATH}: $ENV{'PATH'}")
    if (!$plink_bin_path);

execute_system_command(0,
                           "\nLooking for fastPHASE on your system ... \n");

my $fastphase_cmd = `which fastPHASE`;
chomp $fastphase_cmd;

if ($fastphase_cmd) {
    execute_system_command(0,
                               "\nFound fastPHASE at:\n$fastphase_cmd\n\n");
}
else {
    confess error("Cannot find fastPHASE on your system.\nYour \${PATH}: $ENV{'PATH'}");
}

# Check if filename contains any spaces.
if (defined $filename) {
    $filename =~ s/\s+/\_/g;
    $filename =~ s/\W+/\_/g;
    $filename .= '_' if ($filename !~ m/\_$/);
}
else {
    $filename = 'plink_fastphase_out_';
}

# Store all chr's from map file
store_all_chrs();

##############################################
# All checks complete.                       #
# Run plink to produces fastPHASE output.    #
##############################################


############## This following section of code is buggy. Please do not use
############## -a option with this script
if (defined($all_chrs) && !defined($chr)) {
    foreach my $chr_id (keys %$chrs) {
        $chr = $chr_id;
        print warning("Found chr0 in map file, plink may not properly produce
                      haplotype resulting in wrong analysis downstream
                      Skipping chr0 ..."),
        next
        if ($chr == 0);
        check_intervals();
        #$user_defined_upper_interval = $upper_interval;
        #$user_defined_lower_interval = $lower_interval;
        run_plink_and_convert_plinks_fastphase();
    }
}
#################### Buggy section ends ##################################
elsif (!defined($all_chrs) && defined($chr)) {
    confess error("Found chr0 in map file, plink may not properly produce
                  haplotype resulting in wrong analysis downstream")
        if ($chr == 0);

    check_intervals();
    $user_defined_upper_interval = $upper_interval;
    $user_defined_lower_interval = $lower_interval;
    run_plink_and_convert_plinks_fastphase();
}
elsif (defined($all_chrs) && defined($chr)) {
    confess error("You cannot define both --chr and --all_chrs at the same time!");
}
else {
    confess error("Neither --chr nor --all_chrs is defined!");
}

execute_system_command(0,
                       "\nDone!\n");

# Clearing the shell prompt
print "\n";


##################################### Functions ################################

# Check for over the limit interval error

sub check_intervals {

    $max_coord = index_map_data_and_get_max_coord();
    $max_allowed_kb_limit = $max_coord /  1000; #sprintf("%d", ($max_coord / 1000)) + 1;
    $max_allowed_mb_limit = $max_coord / 1000000; #sprintf("%d", ($max_coord / 1000000)) + 1;

    confess error("You cannot specify --no_of_markers and inteval at the same time")
        if ((interval() == 1) && defined($no_of_markers) && ($no_of_markers ne ''));

    confess error("Only one interval can be defined at a time!\n\n")
        if (interval() > 1);

    confess error("Maximum allowed for chr$chr (plink's --to-kb) : $max_allowed_kb_limit")
        if (($which_metric eq 'kb') && ($lower_interval > $max_allowed_kb_limit));

    confess error("Maximum allowed for chr$chr (plink's --to-mb) : $max_allowed_mb_limit")
        if (($which_metric eq 'mb') && ($lower_interval > $max_allowed_mb_limit));

    confess error("\nMaximum allowed for chr$chr (plink's --to-bp): $max_coord\n")
        if (($which_metric eq 'bp') && ($lower_interval > $max_coord));

    confess error("Neither --no_of_markers nor window size (bp, kb, mb) defined!")
        if ((!defined $no_of_markers || $no_of_markers eq '') && !$which_metric);
    return;
}

# Store all chr numbers from the map file

sub store_all_chrs {
    open(MAP, "<", $map) || confess error("Cannot open $map for reading: $!");
    while (my $line = <MAP>) {
        chomp $line;
        my ($chr_no, $snp_id, $genetic_distance, $chr_pos) = split /\s+/, $line;
        $chr_no = lc($chr_no);
        if (!exists $chrs->{$chr_no}) {
            $chrs->{$chr_no} = 1;
        }
    }
    close MAP;
    return;
}

# Store and return hash ref to map data.
sub index_map_data_and_get_max_coord {
    open(MAP, "<", $map) || confess error("Cannot open $map for reading: $!");
    my $line_no = 0;

    while (my $line = <MAP>) {
        chomp $line;
        my ($chr_no, $snp_id, $genetic_distance, $chr_pos) = split /\s+/, $line;
        $chr_no = lc($chr_no);
        $chr = lc($chr);
        $map_data{$chr_pos} = $line_no if ($chr eq $chr_no);
        $line_no++;
    }

    confess error("Cannot find chr $chr in $map")
        if (!exists $chrs->{$chr});

    close MAP;
    return (sort {$a <=> $b} (keys %map_data))[-1];
}

# Store PED genotype column for each individual id
sub recode_fastPHASE_from_plink{
    my $file_in = shift;
    my $file_out = shift;
    my $marker_end = shift;
    my $marker_beg = shift;
    my $total_markers = shift;

    open (PED, "<", $ped) || confess error("Cannot open $ped for reading: $!");
    open (PLINK_FASTPHASE, "<", $file_in) ||
        confess error("Cannot open $file_in for reading: $!");
    open (FORMATTED_FASTPHASE, '>', $file_out) ||
        confess error("Cannot open $file_out for writing: $!");

    my @coords;
    my $individual_count = <PLINK_FASTPHASE>;
    my $snp_count = <PLINK_FASTPHASE>;

    while (my $line = <PLINK_FASTPHASE>) {
        chomp $line;
        if ($line =~ m/^P\s+?(.+)/) {
            if ($marker_end) {
                if ($marker_end > $total_markers) {
                    my $marker_count = $total_markers - $marker_beg;
                    $marker_end = $total_markers;
                    print FORMATTED_FASTPHASE "$individual_count$marker_count\n";
                    $marker_end--; # adjusting for perl array index
                }
                else {
                    print FORMATTED_FASTPHASE "$individual_count$no_of_markers\n";
                    $marker_end--; # adjusting for perl array index
                }

                print FORMATTED_FASTPHASE "P ";
                my @all_coords = split/\s+/, $1;

                for (my $i=$marker_beg; $i<=$marker_end; $i++) {
                    push (@coords, $all_coords[$i]);

                    if ($i == $marker_end) {
                        print FORMATTED_FASTPHASE $all_coords[$i];
                    }
                    else {
                        print FORMATTED_FASTPHASE "$all_coords[$i] ";
                    }
                }

                print FORMATTED_FASTPHASE "\n";
            }
            else {
                print FORMATTED_FASTPHASE "$individual_count$snp_count";
                print FORMATTED_FASTPHASE "$line\n";
                @coords = split/\s+/, $1;
            }
        }
    }

    while (my $line = <PED>) {
        chomp $line;
        my ($fam_id, $indiv_id, $paternal_id, $maternal_id, $sex, $phenotype, @genotypes) = split /\s+/, $line;
        print FORMATTED_FASTPHASE "# ID $indiv_id\n";
        $indiv_id = lc($indiv_id);
        my $fastphase_haplotype_line_1;
        my $fastphase_haplotype_line_2;

        foreach my $coord (@coords) {
            my $genotype_loc_next = ($map_data{$coord} * 2 ) + 1;
            my $genotype_loc = $genotype_loc_next - 1;
            $fastphase_haplotype_line_1 .= $genotypes[$genotype_loc];
            $fastphase_haplotype_line_2 .= $genotypes[$genotype_loc_next];
        }

        $fastphase_haplotype_line_1 =~ s/0/?/g;
        $fastphase_haplotype_line_2 =~ s/0/?/g;

        print FORMATTED_FASTPHASE "$fastphase_haplotype_line_1\n";
        print FORMATTED_FASTPHASE "$fastphase_haplotype_line_2\n";
    }

    close PED;
    close PLINK_FASTPHASE;
    close FORMATTED_FASTPHASE;
    return;
}

# To check and load modules.
sub check_and_load_modules {
    my $module_list = shift;

    foreach my $module (@$module_list) {
        my $module_installed = eval ("use $module; 1");

        confess error("Required module not installed: $module")
            if (!$module_installed);
    }
    return;
}

# Check if all options entered by user are valid
sub verify_options {
    my $valid_options = shift;
    my $help = shift;

    if (!$valid_options) {
        pod2usage(-verbose => 0,
                  -msg => "\nSee $0 -h for usage/documentation\n");
    }

    if ($help) {
        pod2usage(-verbose => 99,
                  -sections => "OPTIONS",
                  -msg => "\n");
    }

    return;
}

# File checks
sub verify_input_files {
    my $files = shift;

    foreach my $file (@$files) {

        confess error("MAP or PED file not defined: $!")
            if (!defined $file);

        confess error("File ( $file ) does not exist: $!")
            if (!-e $file);

        confess error("File ($file) is empty: $!")
            if (-s $file == 0);

    }
    return;
}

# Check if output path is mentioned, else create output path where the ped
# and map files exists.
sub check_output_path {
    my $path = shift;

    if (defined $path) {

        print "\nAttempting to create $path ...\n\n"
            if (!$quiet && !-d $path);
        system("mkdir -p $path")
            if (defined ($path) && !-d $path);
    }
    else {
        $path = $ENV{'PWD'};
    }

    $path .= '/'
        if ($path !~ m/\/$/);

    return $path;
}


# Verify intervals and assign to variables,
# i.e. more than 1 interval mentioned at the same time
sub interval {
    my $only_one_interval_at_a_time = 0;
    $which_metric = 0;

    $only_one_interval_at_a_time ++,
    ($upper_interval, $lower_interval) = get_intervals($bp_interval_string),
    $which_metric = 'bp',
    $max_interval = $max_coord
        if (defined ($bp_interval_string) && ( $bp_interval_string ne ''));

    $only_one_interval_at_a_time ++,
    ($upper_interval, $lower_interval) = get_intervals($kb_interval_string),
    $which_metric = 'kb',
    $max_interval = $max_allowed_kb_limit
        if (defined ($kb_interval_string) && ( $kb_interval_string ne ''));

    $only_one_interval_at_a_time ++,
    ($upper_interval, $lower_interval) = get_intervals($mb_interval_string),
    $which_metric = 'mb',
    $max_interval = $max_allowed_mb_limit
        if (defined ($mb_interval_string) && ( $mb_interval_string ne ''));

    return 2
        if ($only_one_interval_at_a_time > 1);

    return 1
        if ($only_one_interval_at_a_time == 1);

    return 0
        if ($only_one_interval_at_a_time == 0);
}

# Clean and get intervals
sub get_intervals {
    my $string = shift;
    $string =~ s/\s+//g;

    my @intervals = split/,/, $string;

    confess error("More than 2 limits have been specified for the interval")
        if (scalar(@intervals) > 2);

    confess error ("Upper interval ( $intervals[0] ) cannot be greater than
                   Lower interval ( $intervals[1] )")
        if ($intervals[0] > $intervals[1]);

    return ($intervals[0], $intervals[1]);
}

# Sub routine to get number of coords in fastPHASE
sub get_no_coords {
    my $file = shift;
    open (FASTPHASE, '<', $file) ||
    confess error("Cannot open $file for reading: $!");

    my $no_of_indivs = <FASTPHASE>;
    my $no_of_snps = <FASTPHASE>;

    while (my $line = <FASTPHASE>) {
        chomp $line;
        if ($line =~ m/^(P\s+)/) {
            $line =~ s/$1//;
            $line =~ s/^\s+//;
            $line =~ s/\s+$//;
            my @no_of_positions = split/\s+/, $line;
            confess error("Number of marker positions in P line is not
                          equal to no of snps in second line")
                if (scalar(@no_of_positions) != $no_of_snps);
            return scalar(@no_of_positions);
        }
    }
}

# wrapper code to run plink and convert plinks fastphase
sub run_plink_and_convert_plinks_fastphase {
    my $out_prefix = $fastphase_in_dir . $filename . 'chr' . $chr;
    my ($plink_from, $plink_to, $sliding_window_size);

    $upper_interval = $user_defined_upper_interval;
    $lower_interval = $user_defined_lower_interval;

    # Just extract the portions specified is stop requested.
    $max_interval = $lower_interval if (defined($stop));

    if (!$no_of_markers) {
        $sliding_window_size = $lower_interval - $upper_interval;

        $plink_from = '--from-kb',
        $plink_to = '--to-kb'
            if ($which_metric eq 'kb');

        $plink_from = '--from-mb',
        $plink_to = '--to-mb'
            if ($which_metric eq 'mb');

        $plink_from = '--from-bp',
        $plink_to = '--to-bp'
            if ($which_metric eq 'bp');

        execute_plink($sliding_window_size,
                      $out_prefix,
                      $plink_from,
                      $plink_to);

        if (defined $overlap && $overlap ne '') {
            confess error("--overlap takes either integer or float values only!")
                if ($overlap !~ m/^(\d+|\d+\.\d+)$/);

            $sliding_window_size = $user_defined_lower_interval - $user_defined_upper_interval;
            $lower_interval = $user_defined_upper_interval + $sliding_window_size + $overlap;
            $upper_interval = $user_defined_lower_interval - $overlap;


            # Just extract the portions specified is stop requested.
            $max_interval = $lower_interval if (defined($stop));

            confess error("Cannot calculate new boundaries due to invalid overlap\n
                          Offset: $overlap, Upper Interval: $upper_interval,
                          Lower Interval: $lower_interval")
                if ($upper_interval > $lower_interval);

            execute_plink($sliding_window_size,
                          $out_prefix,
                          $plink_from,
                          $plink_to);
        }
    }
    elsif ($no_of_markers) {
        $out_prefix .= '_markers';
        my $fastphase_to_be_formatted = $out_prefix . '.recode.phase.inp';

        execute_system_command("plink --noweb --ped $ped --map $map --chr $chr --recode-fastphase --out $out_prefix",
                               "\nCommand:\n--------\nplink --noweb --ped $ped --map $map --chr $chr --recode-fastphase --out $out_prefix\n");

        my $total_no_of_markers = get_no_coords($fastphase_to_be_formatted);

        for (my $i=0; $i<$total_no_of_markers; $i=$i+$no_of_markers) {
            my $marker_end = $i + $no_of_markers;
            my $marker_beg = $i;
            my $marker_beg_in_file_name = $marker_beg + 1;
            my $marker_end_in_file_name = $marker_end;
            $marker_end_in_file_name = $total_no_of_markers
                if ($marker_end > $total_no_of_markers);
            my $fastphase_formatted = $out_prefix . '_' . $marker_beg_in_file_name . 'to' . $marker_end_in_file_name . '.formatted.inp';

            execute_system_command(0,
                                   "\nFormatting fastPHASE file from plink and converting to genotype\n\nNew fastPHASE input is:\n$fastphase_formatted\n");
            recode_fastPHASE_from_plink($fastphase_to_be_formatted, $fastphase_formatted, $marker_end, $marker_beg, $total_no_of_markers);
            execute_system_command(0,
                               "\nRunning fastPHASE ...\n");
            run_fastphase($fastphase_formatted);

        }
    }
}

# Subroutine to run fastphase
sub run_fastphase {
    my $file_in = shift;
    my ($filename, $base, $suffix) = fileparse($file_in, qr/\.[^.]*/);
    my $fastphase_out_prefix = $fastphase_out_dir . $filename;

    if (defined ($fastphase) && $fastphase eq '') {
        # default fastphase
        execute_system_command("fastPHASE -o$fastphase_out_prefix $file_in",
                               "\nCommand:\n--------\nfastPHASE -o$fastphase_out_prefix $file_in\n\n");
    }
    elsif (defined ($fastphase) && $fastphase ne '') {
        execute_system_command("fastPHASE $fastphase -o$fastphase_out_prefix $file_in",
                               "\nCommand:\n--------\nfastPHASE $fastphase -o$fastphase_out_prefix $file_in\n\n");
    }

    return;
}

# Sub routine to execute plink
sub execute_plink {

    my $sliding_window_size = shift;
    my $out_prefix = shift;
    my $plink_from = shift;
    my $plink_to = shift;

    for (my $i=$upper_interval; $i<$max_interval; $i=$i+$sliding_window_size) {

        $lower_interval = $i + $sliding_window_size;
        $lower_interval = $max_interval if ($lower_interval > $max_interval);

        my $outfile = $out_prefix . '_' . $i . 'to' . $lower_interval . '_' . $which_metric;
        my $fastphase_to_be_formatted = $outfile . '.recode.phase.inp';
        my $fastphase_formatted = $outfile . '.formatted.inp';


        next if (-f $fastphase_formatted);

        #print "Command:\n--------\nplink --noweb --cow --ped $ped --map $map --chr $chr $plink_from $i $plink_to $lower_interval --recode-fastphase --out $outfile\n";
        execute_system_command("plink --noweb --cow --ped $ped --map $map --chr $chr $plink_from $i $plink_to $lower_interval --recode-fastphase --out $outfile",
                               "\nCommand:\n--------\nplink --noweb --cow --ped $ped --map $map --chr $chr $plink_from $i $plink_to $lower_interval --recode-fastphase --out $outfile\n");
	# Some times plink may be out of bounds, so skip such events
	next if (!-f $fastphase_to_be_formatted);
	
        execute_system_command(0,
                               "\nFormatting fastPHASE file from plink and converting to genotype\n\nNew fastPHASE input is:\n$fastphase_formatted\n");
        recode_fastPHASE_from_plink($fastphase_to_be_formatted, $fastphase_formatted);
        execute_system_command(0,
                               "\nRunning fastPHASE ...\n");
        run_fastphase($fastphase_formatted);
    }
    return;
}

# Shell msg that differentiates log from error
sub error {
    my $msg = shift;
    print "\nERROR!\n------\n$msg\n\n";
    pod2usage(-verbose => 0,
	      -msg => "\nSee $0 -h for usage/documentation\n");
    return;
}

# Shell msg that differentiates log from warning
sub warning {
    my $msg = shift;
    print "\nWarning!\n------\n$msg\n\n";
    return;
}

# Subroutine to execute system command with verbosity
sub execute_system_command {
    my $command = shift;
    my $print_msg = shift;

    if (!$quiet) {
       print $print_msg if ($print_msg);
       system ("$command") if ($command);
    }
    elsif ($quiet) {
        system ("$command 1>/dev/null 2>/dev/null") if ($command);
    }

    return;
}

__END__

=head1 NAME

recode_plinks_fastPHASE_and_run.pl

=head1 SYNOPSIS

Convert fastPHASE output from plink to genotype data.

=head2 EXAMPLES

    perl recode_plinks_fastPHASE_and_run.pl -h

    perl recode_plinks_fastPHASE_and_run.pl

=head2 DESCRIPTION

plink has an option to generate input file for fastPHASE program from map and
ped files, but it does not retain genotype information, instead replaces it
with 1's and 0's. This code will retain that information when both map and ped
files are supplied. Also, this is a wrapper script to generate fastPHASE input
files at required intervals. From testing, it was found that a haplotype
containing 10 genotype loci is optimum to work with fastPHASE. Alternately, one
can use plink's options such as --from-kb  / --to-kb, --from-mb / --to-mb and
--from-bp / --to-bp to extract snps within those boundaries.

=head2 *** IMPORTANT ***

It is required that the order of markers in map file is in order with the
genotypes column (column 7) of the ped file, as mentioned in plink's
documentation. Also, this code was written to work with plink data for diploid
organisms.

=head1 OPTIONS

=over 5

=item --help or -h (Optional).

    Displays this helpful message.

=item --quite or -q (Optional).

    Supress log message to shell window.

    Default: disabled.

=item --map or -ma (Required).

    Path to MAP file in plink format.

=item --ped or -p (Required).

    Path to PED file in plink format.

=item --all_chrs or -a (Required).

    Run wrapper for all chromosomes present in the map file

=item --chr or -c (Required if --all_chrs not defined).

    plink option: --chr: Select a particular chromosome {N}

=item --fastphase or -fa (Optional).

    Use this script to run fastphase (Will take longer)
    -fa Default fastphase without options.
    -fa '-T10' fastphase run with options.

=item --no-of-markers or -n (Optional).

    Required if none of intervals is specified.
    How many makers should be considered for each chromosome (starting with the
    first marker) to create a fastPHASE input file.
    -n 10  will create fastPHASE input file with  a pair of 10 genotypes per
    line for each diploid individual.

    Default: If -n or -kb or -mb or -bp is not defined, the program
             exits with error.

=item  --bp_interval or -bp (Optional).

    plink options: --from-bp, --to-bp: Select SNPs within this window ..
                                       specified in bases.

    Default: If -n or -bp or -kb or -mb is not defined, all the data is
             converted into a single fastPHASE input file.

=item --kb_interval or -kb (Optional).

    plink options: --from-kb, --to-kb: Select SNPs within this window ..
                                      specified in kilo bases.

    Default: If -n or -bp or -kb or -mb is not defined, all the data is
             converted into a single fastPHASE input file.

=item  --mb_interval or -mb (Optional).

    plink option: --from-mb, --to-mb: Select SNPs within this window ..
                                      specified in mega bases.

    Default: If -n or -bp or -kb or -mb is not defined, all the data is
             converted into a single fastPHASE input file.

=item --overlap or -ov (Optional).

    Produce fastphase input files with the specified overlap.
    ********** IMPORTANT **********
    You are responsible for specifying overlap in same units corresponding to the
    intervals. For ex: if you specified -mb '0,1' and you want an overlap
    of 500kb you must use --overlap 0.5 instead of --overlap 500.

=item --stop or -s (Optional).

    Just extract the specified intervals of the chromosome. If -mb '0, 1' is
    specified, plink and fastPHASE will be run on markers extracted from those
    portions only.

=item --filename or -fi (Optional)

    Define what should be the base name of the output files. If it is not
    defined, all files will have a base name 'plink_fastphase_out'.

=item --out (Optional).

    Path to store fastPHASE input and file(s). If not mentioned, a directory
    will be created in current working path.

    Default: disabled.

=back

=head1 AUTHOR

Kranti Konganti, E<lt>konganti@tamu.eduE<gt>.

=head1 COPYRIGHT

This program is distributed under the Artistic License.

=head1 DATE

Apr-19-2012

=cut
