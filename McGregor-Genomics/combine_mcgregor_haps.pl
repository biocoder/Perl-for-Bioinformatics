#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

my ($CHANGEDBY) = q$LastChangedBy: konganti$ =~ m/.+?\:(.+)/;
my ($LASTCHANGEDDATE) = q$LastChangedDate: 2013-09-16 08:33:47 -0500 (Mon, 16 Sep 2013) $ =~ m/.+?\:(.+)/;
my ($VERSION) = q$LastChangedRevision: 57 $ =~ m/.+?(\d+)/;
my $AUTHORFULLNAME = 'Kranti Konganti';

# Being extra cautious for nice die message instead of Perl's message

check_and_load_modules(['Pod::Usage', 'Getopt::Long', 'File::Basename']);
this_script_info();

# Declare initial global variables

my ($quiet, $breed, $output, $help, $mcgregor_pop, $chr_range, @chrs,
    $start_chr_line, $print_seq, $file_ext, $test_animal_id);
my $parent_utf8->{'S'} = 'S';#"\xb9";
$parent_utf8->{'D'} = 'D';#"\xb2";
$parent_utf8->{''} = 'U';#"\xb3";
my $breed_utf8->{'N'} = 'N';
$breed_utf8->{'A'} = 'A';
$breed_utf8->{''} = '-';

my $BREED_OF_ORIGIN_SCRIPT = 'assign_BreedOfOrigin_from_fastPHASE_for_McGregor_pop.pl';
my $BREED_NELLORE = 'N';
my $BREED_ANGUS = 'A';
my $PARENT_SIRE = 'S';
my $PARENT_DAM = 'D';
my $RECOMBINANT = '';

my $is_valid_option = GetOptions ('help|?' => \$help,
                                  'quiet' => \$quiet,
                                  'breed=s' => \$breed,
                                  'output=s' => \$output,
                                  'mcgregor-pop=s' => \$mcgregor_pop,
				  'chr-range=s' => \$chr_range,
                                  'print-seq' => \$print_seq,
                                  'test-animal-id=s' => \$test_animal_id
                                  );

if (defined $chr_range && $chr_range ne '') {
    if ($chr_range =~ m/^\w+\.\.\w+$/) {
        my ($chr_st, $chr_end) = split/\.\./, $chr_range;
        @chrs = $chr_st..$chr_end;
    }
    elsif ($chr_range =~ m/\d+/) {
        $chrs[0] = $chr_range;
    }
    else {
        error("Cannot recognize chromosome id or chromosome range: $chr_range");
    }
}
else {
    @chrs = 1..30;
}

# Check for the validity of options
verify_options($is_valid_option, $help);
verify_input_files([$mcgregor_pop], ['McGregor Population list ( mcgregor_pop_animal_ids.txt )']);

$output = validate_create_path($output, 'create', 'Output');
$breed = validate_create_path($breed, 'do not create', 'Breed of Origin Directory');

execute_system_command(0,
                       "\nChecking the version requirement for system level commands ...\n");
check_sys_level_cmds(['grep', 'ls', 'sort', 'basename'],
                     ['2.10', '8.13', '8.13', '8.13']);

execute_system_command(0,
                       "\nOutput files will be stored at $output ...\n");

my $animal_info = store_hapguess($mcgregor_pop);

execute_system_command(0,
                       "\nCreating Genotype files ...\n");

$file_ext = '.bo.txt' if (!defined($print_seq));
$file_ext = '.seq.txt' if (defined($print_seq));

foreach my $chr (sort {$a <=> $b} @chrs) {
    my $filename = $output . 'chr' . $chr . $file_ext;
    my $genotype_out_fh = open_file($filename, '>');

    foreach my $animal (keys %$animal_info) {

        next if (defined($test_animal_id) && $test_animal_id ne '' && $animal !~ m/$test_animal_id/i);

        my ($ah1, $ah2) = calc_and_combine_hap_overlap($animal, $chr);
        $ah1 =~ s/S|D|U//g;
        $ah2 =~ s/S|D|U//g;

        print $genotype_out_fh ">$animal | Haplotype 1\n$ah1\n\n>$animal | Haplotype 2\n$ah2\n\n\n\n";
    }
    close $genotype_out_fh;
}

########################## Functions ##########################################

# Magical function that calculates overlap and combines over it.

sub calc_and_combine_hap_overlap {
    my $animal = shift;
    my $chr = shift;
    my $chr_boa_files = get_files($chr);
    my $overlap_files = 0;
    $start_chr_line = 0;

    my ($a1h1, $a1h2, $a2h1, $a2h2, $marker_line_1, $marker_line_2,
        $a1_boa1, $a1_poa1, $a2_boa1, $a2_poa1,
        $a1_boa2, $a1_poa2, $a2_boa2, $a2_poa2,
        $combined_hap1, $combined_hap2,
        $animal_chr_hap1, $animal_chr_hap2,
        @conc_marker_line, $long_a2h1_len,
        $long_a2h1_no_overlap, $long_a2h2_no_overlap) = '';

    foreach my $chr_boa (@$chr_boa_files) {
        if ($chr_boa =~ m/\.bo\.txt$/) {

            #print "$animal\t$chr_boa";

            my $grep_cmd = "grep -P '" . $animal_info->{$animal} . $animal . "\t.*'";
            $grep_cmd =~ s/(.+?\'.+?\t)(\W.+)/$1\\$2/;
            $grep_cmd =~ s/\t/\\t/g;
            my $animal_lines = qx($grep_cmd $chr_boa);

            error("Unable to get animal haplotypes !!\nFailed at command: $grep_cmd $chr_boa\n")
                if (!$animal_lines || $animal_lines eq '');

            $overlap_files++;

            if ($overlap_files == 3) {
                if ($long_a2h1_len != 0) {
                    $long_a2h1_no_overlap = reverse($long_a2h1_no_overlap);
                    $long_a2h2_no_overlap = reverse($long_a2h2_no_overlap);

                    for (1 .. $long_a2h1_len) {
                        chop $long_a2h1_no_overlap;
                        chop $long_a2h2_no_overlap;
                        shift @$marker_line_2;
                    }

                    $long_a2h1_no_overlap = reverse($long_a2h1_no_overlap);
                    $long_a2h2_no_overlap = reverse($long_a2h2_no_overlap);
                }

                $a1h1 .= $long_a2h1_no_overlap;
                $a1h2 .= $long_a2h2_no_overlap;

                push(@conc_marker_line, @$marker_line_2);
                $marker_line_1 = \@conc_marker_line;
                $overlap_files = 2;
            }

            if ($overlap_files == 1) {
                ($a1h1, $a1h2, $a1_boa1, $a1_poa1,
                 $a1_boa2, $a1_poa2) = sort_animal_haps($animal_lines);
                $marker_line_1 = get_markers($chr_boa);
                push(@conc_marker_line, @$marker_line_1);
            }
            elsif ($overlap_files == 2) {
                ($a2h1, $a2h2, $a2_boa1, $a2_poa1,
                 $a2_boa2, $a2_poa2) = sort_animal_haps($animal_lines);
                $marker_line_2 = get_markers($chr_boa);

                my ($a1h1_m, $a2h1_m,
                    $overlap_start_index_1) = get_overlap_info($marker_line_1, $marker_line_2,
                                                               \$a1h1, \$a2h1);
                my ($a1h2_m, $a2h2_m,
                    $overlap_start_index_2) = get_overlap_info($marker_line_1, $marker_line_2,
                                                               \$a1h2, \$a2h2);

                #print scalar(@$marker_line_1), "\t";
                #print length($a1h1), " - ", length($a2h1), "\n";
                #print "\n\n\n@$marker_line_1\n\n\n@$marker_line_2\n\n\n";

                if (!$a1h1_m || !$a2h2_m) {
                    $long_a2h1_len = 0;
                    $long_a2h1_no_overlap = $a2h1;
                    $long_a2h2_no_overlap = $a2h2;

                    if ($start_chr_line < 1) {
                        for (1 .. length($a1h1)) {
                            $animal_chr_hap1 .= $breed_utf8->{$a1_boa1}
                        }
			for (1 .. length($a1h2)) {
			    $animal_chr_hap2 .= $breed_utf8->{$a1_boa2};
			}
			for (1 .. length($a2h1)) {
			    $animal_chr_hap1 .= $breed_utf8->{$a2_boa1};
			}
			for (1 .. length($a2h2)) {
			    $animal_chr_hap2 .= $breed_utf8->{$a2_boa2};
			}
			$start_chr_line++;
		    }
		    else {
			for (1 .. length($a2h1)) {
                            $animal_chr_hap1 .= $breed_utf8->{$a2_boa1};
			}
			for (1 .. length($a2h2)) {
			    $animal_chr_hap2 .= $breed_utf8->{$a2_boa2};
			}
		    }
                    next;
                }

                my ($preserve_a2h1, $preserve_a2h2) = ($a2h1, $a2h2);
                my ($preserve_a2boa1, $preserve_a2boa2,
                    $preserve_a2poa1, $preserve_a2poa2) = ($a2_boa1, $a2_boa2,
                                                           $a2_poa1, $a2_poa2);

                if ($a1h1_m ne $a2h1_m && $a1h2_m ne $a2h2_m) {
                    ($a1h1_m, $a2h1_m,
                    $overlap_start_index_1) = get_overlap_info($marker_line_1, $marker_line_2,
                                                                \$a1h1, \$a2h2);
                    ($a1h2_m, $a2h2_m,
                    $overlap_start_index_2) = get_overlap_info($marker_line_1, $marker_line_2,
                                                                \$a1h2, \$a2h1);
                    $a2h1 = $a2h2;
                    $a2h2 = $preserve_a2h1;
                    $a2_poa1 = $a2_poa2;
                    $a2_poa2 = $preserve_a2poa1;
                    $a2_boa1 = $a2_boa2;
                    $a2_boa2 = $preserve_a2boa1;
                }
                #else {
                    #$a2h1 = $preserve_a2h1;
                    #$a2h2 = $preserve_a2h2;
                    #$a2_poa1 = $preserve_a2poa1;
                    #$a2_poa2 = $preserve_a2poa2;
                    #$a2_boa1 = $preserve_a2boa1;
                    #$a2_boa2 = $preserve_a2boa2;
                #}

                #if ($a1h2_m ne $a2h2_m) {
                    #($a1h2_m, $a2h2_m,
                    # $overlap_start_index_2) = get_overlap_info($marker_line_1, $marker_line_2,
                    #                                            \$a1h2, \$a2h1);
                    #$a2_boa2 = $a2_boa1;
                    #$a2_poa2 = $a2_poa1;
                #}

                $long_a2h1_len = length($a2h1_m);
                $long_a2h1_no_overlap = $a2h1;
                $long_a2h2_no_overlap = $a2h2;

# JUNK CODE. Not necessary, but retaining it for sanity checks.
#
#                if (($a1_boa1 eq $a2_boa1) && ($a1_poa1 eq $a2_poa1)
#                    && ($a1h1_m eq $a2h1_m)) {
#                    $combined_hap1 = combine_genotype($a1h1,
#                                                      $a2h1,
#                                                      $a1_boa1,
#                                                      $a1_poa1,
#                                                      $overlap_start_index_1,
#                                                      "1"
#                                                      );
#                }
#                elsif (($a1_boa1 eq $a2_boa1) && ($a1_poa1 ne $a2_poa1)
#                    && ($a1h1_m eq $a2h1_m)) {
#                    my $other_parent_present = check_other($a1_poa1, $a2_poa1);
#                    $combined_hap1 = combine_genotype($a1h1,
#                                                      $a2h1,
#                                                      $a1_boa1,
#                                                      $other_parent_present,
#                                                      $overlap_start_index_1,
#                                                      "1"
#                                                      );
#                }
#               elsif (($a1_boa1 ne $a2_boa1) && ($a1_poa1 eq $a2_poa1)
#                    && ($a1h1_m eq $a2h1_m)) {
#                    my $other_breed_present = check_other($a1_boa1, $a2_boa1);
#                    $combined_hap1 = combine_genotype($a1h1,
#                                                      $a2h1,
#                                                      $other_breed_present,
#                                                      $a1_poa1,
#                                                      $overlap_start_index_1,
#                                                      "1"
#                                                      );
#                }
#                elsif (($a1_boa1 ne $a2_boa1) && ($a1_poa1 ne $a2_poa1)
#                    && ($a1h1_m eq $a2h1_m)) {
#                    my $other_breed_present = check_other($a1_boa1, $a2_boa1);
#                    my $other_parent_present = check_other($a1_poa1, $a2_poa1);
#
#                    $combined_hap1 = combine_genotype($a1h1,
#                                                      $a2h1,
#                                                      $other_breed_present,
#                                                      $other_parent_present,
#                                                      $overlap_start_index_1,
#                                                      "1"
#                                                      );
#                }
#                #else {
#                    #print "1: elsed\n$chr_boa\t\n$a1h1_m\n$a1h2_m\n\n$a2h1_m\n$a2h2_m\n\n";
#                #}
#
#
#                if (($a1_boa2 eq $a2_boa2) && ($a1_poa2 eq $a2_poa2)
#                    && ($a1h2_m eq $a2h2_m)) {
#                    $combined_hap2 = combine_genotype($a1h2,
#                                                      $a2h2,
#                                                      $a1_boa2,
#                                                      $a1_poa2,
#                                                      $overlap_start_index_2,
#                                                      "2"
#                                                      );
#                }
#                elsif (($a1_boa2 eq $a2_boa2) && ($a1_poa2 ne $a2_poa2)
#                    && ($a1h2_m eq $a2h2_m)) {
#                    my $other_parent_present = check_other($a1_poa2, $a2_poa2);
#                    $combined_hap2 = combine_genotype($a1h2,
#                                                      $a2h2,
#                                                      $a1_boa2,
#                                                      $other_parent_present,
#                                                      $overlap_start_index_2,
#                                                      "2"
#                                                      );
#                }
#                elsif (($a1_boa2 ne $a2_boa2) && ($a1_poa2 eq $a2_poa2)
#                    && ($a1h2_m eq $a2h2_m)) {
#                    my $other_breed_present = check_other($a1_boa2, $a2_boa2);
#                    $combined_hap2 = combine_genotype($a1h2,
#                                                      $a2h2,
#                                                      $other_breed_present,
#                                                      $a1_poa2,
#                                                      $overlap_start_index_2,
#                                                      "2"
#                                                      );
#                }
#                elsif (($a1_boa2 ne $a2_boa2) && ($a1_poa2 ne $a2_poa2)
#                    && ($a1h2_m eq $a2h2_m)) {
#                    my $other_breed_present = check_other($a1_boa2, $a2_boa2);
#                    my $other_parent_present = check_other($a1_poa2, $a2_poa2);
#                    $combined_hap2 = combine_genotype($a1h2,
#                                                      $a2h2,
#                                                      $other_breed_present,
#                                                      $other_parent_present,
#                                                      $overlap_start_index_2,
#                                                      "2"
#                                                      );
#                }
#                $animal_po->{$animal}->{$a2_poa1}++;
#                #else {
#                    #print "2: elsed\n$a1h1_m\n$a1h2_m\n\n$a2h1_m\n$a2h2_m\n\n";
#                #}
#
#                #$overlap_files = 0;
#                #binmode(STDOUT, ":utf8");
#
#                #print "$a1h1, $a1h2\n$a2h1, $a2h2\n\n$a1h1_m, $a1h2_m\n$a2h1_m, $a2h2_m\n\n";
#                #exit;
#
#                if (!$combined_hap1) {
#                   $combined_hap1 =  combine_genotype($a1h1,
#                                                      $a2h1,
#                                                      '',
#                                                      '',
#                                                      $overlap_start_index_1,
#                                                      "1"
#                                                      );
#                   if (defined ($print_seq)) {
#                        $animal_chr_hap1 = $combined_hap1;
#                    }
#                    else  {
#                        #$animal_chr_hap1 = $combined_hap1;
#                        #print "else1 - h1: $animal_chr_hap1 - $combined_hap1\n";
#                    }
#                }
#                if (!$combined_hap2) {
#                   $combined_hap2 = combine_genotype($a1h2,
#                                                     $a2h2,
#                                                     '',
#                                                     '',
#                                                     $overlap_start_index_2,
#                                                     "2"
#                                                     );
#                   if (defined ($print_seq)) {
#                        $animal_chr_hap2 = $combined_hap2;
#                    }
#                    else  {
#                        #$animal_chr_hap2 = $combined_hap2;
#                        #print "else1 - h1: $animal_chr_hap2 - $combined_hap2\n";
#                    }
#                }
#                #elsif ($combined_hap1 =~ m/S$/ && $combined_hap2 =~ m/D$/) {
#                #    $animal_chr_hap1 = $combined_hap1;
#                #    $animal_chr_hap2 = $combined_hap2;
#                #}
#                #elsif ($combined_hap1 =~ m/D$/ && $combined_hap2 =~ m/S$/) {
#                #    $animal_chr_hap1 = $combined_hap2;
#                #    $animal_chr_hap2 = $combined_hap1;
#                #}
#                #elsif (($combined_hap1 =~ m/S$/ && $combined_hap2 =~ m/U$/)
#                #       || ($combined_hap1 =~ m/U$/ && $combined_hap2 =~ m/D$/)) {
#                #    $animal_chr_hap1 = $combined_hap1;
#                #    $animal_chr_hap2 = $combined_hap2;
#                #}
#                #elsif (($combined_hap1 =~ m/U$/ && $combined_hap2 =~ m/S$/)
#                ##       || ($combined_hap1 =~ m/D$/ && $combined_hap2 =~ m/U$/)) {
#                 #   $animal_chr_hap1 = $combined_hap2;
#                 #   $animal_chr_hap2 = $combined_hap1;
#                #}
#                #elsif ($combined_hap1 =~ m/U$/ && $combined_hap2 =~ m/U$/
#                #       && $combined_hap1 ne '' && $combined_hap2 ne '') {
#                #    $animal_chr_hap1 = $combined_hap1;
#                #    $animal_chr_hap2 = $combined_hap2;
#                #}
#                #else {
#                #    $combined_hap1 = combine_genotype($a1h1,
#                #                                      $a2h1,
#                #                                      '',
#                #                                      '',
#                #                                      $overlap_start_index_1
#                #                                      );
#                #    $combined_hap2 = combine_genotype($a1h2,
#                #                                      $a2h2,
#                #                                      '',
#                #                                      '',
#                #                                      $overlap_start_index_2
#                #                                      );
#                #    $animal_chr_hap1 = $combined_hap1;
#                #    $animal_chr_hap2 = $combined_hap2;
#                #}
#                else {
#                    if (defined ($print_seq)) {
#                        $animal_chr_hap1 = $combined_hap1;
#                        $animal_chr_hap2 = $combined_hap2;
#                    }
#                    else  {
#                        #$animal_chr_hap1 = $combined_hap1;
#                        #$animal_chr_hap2 = $combined_hap2;
#                        #print "else2 - h1: $animal_chr_hap1 - $combined_hap1\n";
#                        #print "else2 - h2: $animal_chr_hap2 - $combined_hap2\n";
#                    }
#                }
#                #else {
#                    #print "\n$chr_boa$animal\n\n$a1h1 - $a1_boa1|$a1_poa1\t$a1h2 - $a1_boa2|$a1_poa2\n\n$a2h1 - $a2_boa1|$a2_poa1\t$a2h2 - #$a2_boa2|$a1_poa2\n\n1:$combined_hap1\n2:$combined_hap2\n\n";
#                    #exit;
#               # }
#                #$animal_chr_hap1 .= $combined_hap1;
#                #$animal_chr_hap2 .= $combined_hap2;
#                #print "\n\n";
#                #print "$combined_hap1\n";
#                #print "$combined_hap2\n";
#                #exit;
#                my $long_poa1_place_holder = $combined_hap1;
#                #print $combined_hap2, "\n";
#                $long_poa1 = chop $long_poa1_place_holder;
#                $long_boa1 = chop $long_poa1_place_holder;
#                my $long_poa2_place_holder = $combined_hap2;
#                $long_poa2 = chop $long_poa2_place_holder;
#                $long_boa2 = chop $long_poa2_place_holder;
#                #print "\n-----\n$combined_hap1\n$combined_hap2\n-----\n";

                #if ($breed_utf8->{$a2_boa1} eq '-') {
                    #print "@$marker_line_2\n\n\n";
                #}

                if ($start_chr_line < 1) {
                    for (1 .. length($a1h1)) {
                        $animal_chr_hap1 .= $breed_utf8->{$a1_boa1};
                    }
                    for (1 .. length($a1h2)) {
                        $animal_chr_hap2 .= $breed_utf8->{$a1_boa2};
                    }
                    for ($overlap_start_index_1 .. length($a2h1)) {
                        $animal_chr_hap1 .= $breed_utf8->{$a2_boa1};
                    }
                    for ($overlap_start_index_2 .. length($a2h2)) {
                        $animal_chr_hap2 .= $breed_utf8->{$a2_boa2};
                    }
                    $start_chr_line++;
                }
                else {
                    for ($overlap_start_index_1 .. length($a2h1)) {
                        $animal_chr_hap1 .= $breed_utf8->{$a2_boa1};
                    }
                    for ($overlap_start_index_2 .. length($a2h2)) {
                        $animal_chr_hap2 .= $breed_utf8->{$a2_boa2};
                    }
                }
            }
        }
        else {
            confess error("This file ( $chr_boa ) does not contain Breed of Origin / Parent of Origin Information");
        }
    }

    return ($a1h1, $a1h2) if (defined ($print_seq));
    return ($animal_chr_hap1, $animal_chr_hap2) if (!defined ($print_seq));
}

# JUNK FUNCTIONS. Retaining it for sanity checks.
#
#sub check_other {
#    my $check_1 = shift;
#    my $check_2 = shift;
#
#    my $one_known;
#    if ((!$check_1 || $check_1 eq '') && ($check_2 || $check_2 ne '')) {
#        $one_known = $check_2;
#    }
#    elsif ((!$check_2 || $check_2 eq '') && ($check_1 || $check_1 ne '')) {
#        $one_known = $check_1;
#    }
#    else {
#        $one_known = '';
#    }
#    return $one_known;
#}
#
#sub combine_genotype {
#    my $hap_1 = shift;
#    my $hap_2 = shift;
#    my $hap_length_1 = length($hap_1);
#    my $hap_length_2 = length($hap_2);
#    my $hap = '';
#    my $boa = shift;
#    my $poa = shift;
#    my $overlap_start = shift;
#    my $ishap1ORhap2 = shift;
#
    # Hap 1
#    if ($print_seq) {
#        $hap = $hap_1;
#    }
#    else {
#        for (1 .. $hap_length_1) {
#            $hap .= $breed_utf8->{$boa};
#        }
#    }
#
    # Hap 2
#    if ($print_seq) {
#        my @hap_nucs = split//, $hap_2;
#        for ($overlap_start .. $hap_length_2) {
#            $_--;
#            $hap .= $hap_nucs[$_];
#        }
#    }
#    else {
#        for ($overlap_start .. $hap_length_2) {
#            $hap .= $breed_utf8->{$boa};
#        }
#    }
#
#    $hap .= $breed_utf8->{$boa};
#    $hap .= $parent_utf8->{$poa};
#
#    return $hap;
#}

# Subroutine to get overlap indexes for coordinate line.

sub get_overlap_info {
    my $marker_line_1 = shift;
    my $marker_line_2 = shift;
    my $hap_1 = shift;
    my $hap_2 = shift;

    $$hap_1 =~ s/\s+//g;
    $$hap_2 =~ s/\s+//g;

    my %markers_line_1;

    for (0 .. $#$marker_line_1) {
        $markers_line_1{@$marker_line_1[$_]} = $_;
    }

    foreach my $marker (sort @$marker_line_2) {
        if (exists $markers_line_1{$marker}) {
            my $marker_line_2_extract_size = scalar(@$marker_line_1) - $markers_line_1{$marker};
            my $marker_line_1_substr_end = substr($$hap_1, $markers_line_1{$marker});
            my $marker_line_2_substr_beg = substr($$hap_2, 0, $marker_line_2_extract_size);
            my $overlap_start_index = $marker_line_2_extract_size + 1;
            return ($marker_line_2_substr_beg, $marker_line_1_substr_end, $overlap_start_index);
        }
    }
    return (0, 0, 0);
}

# Uses Unix's grep tool to get marker line

sub get_markers {
    my $file = shift;
    my $grep_cmd = 'grep -ohP "Animal_ID\t[\d\s]*\tHapBreedOrigin"';
    my $marker_loc_line = `$grep_cmd $file`;
    chomp $marker_loc_line;
    $marker_loc_line =~ s/Animal_ID\t//;
    $marker_loc_line =~ s/\tHapBreedOrigin//;
    my @marker_locs = split/\s+/, $marker_loc_line;
    return (\@marker_locs);
}

# Sorts animal haps by hap index number generated by assign_BreedOfOrigin_from_fastPHASE_for_McGregor_pop

sub sort_animal_haps {
    my $animal_hap_lines = shift;
    my ($animal_hap_1, $animal_hap_2,
        $boa_hap_1, $poa_hap_1,
        $boa_hap_2, $poa_hap_2);
    chomp $animal_hap_lines;

    my @animal_haps = split/[\n\r]/, $animal_hap_lines;

    foreach my $animal_hap (@animal_haps) {
        if ($animal_hap =~ m/HAP0/) {
            #print $animal_hap, "\n";
            #exit;
            (my $sire, my $dam, my $generation, my $sex, my $animal, $animal_hap_1, $boa_hap_1, $poa_hap_1) = split/\t/, $animal_hap;
            $animal_hap_1 =~ s/HAP0\://;
        }
        else {
            #print $animal_hap, "\n";
            (my $sire, my $dam, my $generation, my $sex, my $animal, $animal_hap_2, $boa_hap_2, $poa_hap_2) = split/\t/, $animal_hap;
            $animal_hap_2 =~ s/HAP1\://;
        }
    }
    return ($animal_hap_1, $animal_hap_2, $boa_hap_1, $poa_hap_1, $boa_hap_2, $poa_hap_2);
}

# Uses Unix's sort tool to sort files generated by assign_BreedOfOrigin_from_fastPHASE_for_McGregor_pop
# by chromosome number and overlap information.

sub get_files {
    my $chr = shift;
    my $boa_file_prefix = $breed . 'mcgregor_chr' . $chr . '_*';
    my $sort_k_param = '-k1.' . length($boa_file_prefix) . 'n';
    my @boa_files = `ls $boa_file_prefix | sort $sort_k_param`;

    confess error("Cannot find breed of origin files. Unable to calculate overlap !\nDid you run $BREED_OF_ORIGIN_SCRIPT script for chr$chr to assign Breed Of Origin and Parent of Origin for each animal haplotype ?\n")
        if ($#boa_files <= 0);

    return \@boa_files;
}

# Stores happlotype info and returns a hashref

sub store_hapguess {
    my $file = shift;
    my $fh = open_file($file, '<');
    my $hapguess;

    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ m/(.+?\t.+?\t.+?\t.+?\t)(.+)/) {
            $hapguess->{$2} = $1;
        }
    }

    close $fh;
    return $hapguess;
}

# Does what the function's name says

sub strip_leading_and_trailing_spaces {
    my $line = shift;
    $line =~ s/^\s+//;
    $line =~ s/\s+$//;
    $line = uc($line); # upper case all nucleotides
    return $line;
}

# To check and load modules.

sub check_and_load_modules {
    my $module_list = shift;
    my $is_module_loadable = 1;
    my $req_modules = '';

    foreach my $module (@$module_list) {
        my $module_installed = eval ("use $module; 1");
	$is_module_loadable = 0,
	$req_modules .= "$module, "
	    if (!$module_installed);
    }

    $req_modules =~ s/\,\s+$//;

    confess error("Required module(s) not installed. Following modules and its dependencies must be installed at system level:\n$req_modules\n")
	if (!$is_module_loadable);

    return;
}

# Check if all options entered by user are valid

sub verify_options {
    my $valid_options = shift;
    my $help = shift;

    if (!$valid_options) {
        pod2usage(-verbose => 2,
                  -msg => "\nSee $0 -h for options.\n");
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
    my $what_files = shift;
    my $file_no = 0;

    foreach my $file (@$files) {

        confess error("@$what_files[$file_no] file not specified: $!")
            if (!defined $file);

        confess error("@$what_files[$file_no] file ( $file ) does not exist: $!")
            if (!-e $file);

        confess error("@$what_files[$file_no] file ($file) is empty: $!")
            if (-s $file == 0);

        # Removing DOS Carriage returns
        `perl -i -p -e 's/\r/\n/g' $file`;
        `perl -i -p -e 's/^\n//' $file`;

        $file_no++;
    }
    return;
}

# Check if output path is mentioned, else create output path where the ped
# and map files exists.

sub validate_create_path {
    my $path = shift;
    my $create_dir = shift;
    my $msg = shift;

    confess error ("$msg path not defined or entered!")
        if (!defined $path || $path eq '');

    if ($create_dir eq 'create') {
        if (defined $path) {
            execute_system_command("mkdir -p $path",
                                   "\nAttempting to create $path ...\n\n",
                                   )
                if (defined ($path) && !-d $path);
        }
        else {
            $path = $ENV{'PWD'};
        }
    }

    confess error("Path ( $path ) does not exist: $!")
        if (!-d $path);

    $path .= '/'
        if ($path !~ m/\/$/);

    return $path;
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

# Subroutine to check the least required version of system command.

sub check_sys_level_cmds {
    my $cmds = shift;
    my $versions = shift;
    my ($curr_version_unf, $req_version);
    for (my $i=0; $i<scalar(@$cmds); $i++) {
        my $cmd_version_out = `@$cmds[$i] --version`;
        $req_version = @$versions[$i];
        $req_version =~ s/\.//g;
        if ($cmd_version_out =~ m/^@$cmds[$i].*?(\d+[\.\d]*)/) {
            my $curr_version = $curr_version_unf = $1;
            $curr_version =~ s/\.//g;
            if ($curr_version < $req_version) {
                confess error("\nAt least Version @$versions[$i] required for system level command: @$cmds[$i]\nCurrent Version: $curr_version_unf\n");
            }
        }
    }
    return;
}

# Shell msg that differentiates log from error

sub error {
    my $msg = shift;
    print "\nERROR!\n------\n$msg\n\n";
    exit;
}

# Shell msg for warning

sub warning {
    my $msg = shift;
    print "\nWARNING!\n--------\n$msg\n\n";
    return;
}

# Subroutine to open files and return file handle

sub open_file {
    my $file = shift;
    my $mode = shift;
    open (my $file_handle, "$mode", $file) ||
        confess error("Cannot open ( $file ) in mode ( $mode ): $!");
    return $file_handle;
 }

# Subroutine to print SCRIPT Version
# Print this script's info

sub this_script_info {
    print "\n", '@', '-' x 80, '@', "\n";
    print "  Program Name       :  " , basename($0), "\n";
    print "  Version            :  $VERSION\n" if ($VERSION);
    print "  Author             :  $AUTHORFULLNAME\n" if ($AUTHORFULLNAME);
    print "  Last Changed By    : $CHANGEDBY\n" if ($CHANGEDBY);
    print "  Last Changed Date  : $LASTCHANGEDDATE\n";
    print '@', '-' x 80, '@', "\n\n";
    return;
}

__END__

=head1 NAME

combine_mcgregor_haps.pl

=head1 SYNOPSIS

This script will attempt to combine overlapping regions and tries to determine
the breed of origin over the strech of markers. If the breed of origin is Nellore,
it marks the location as 'N' and as 'A' for Angus. It leaves '-' if it is
unable to determine the breed of origin for the overlapped region. This may be
attributed due to possible genotyping errors. This script runs on the output
files from assign_BreedOfOrigin_from_fastPHASE_for_McGregor_pop.pl.

Examples:

    perl combine_mcgregor_haps.pl -h

    generate_mcgregor_genotype_imp.pl -m mcgregor_pop_animal_ids.txt -o <output dir to store results from this script> -b <assign_BreedOfOrigin_from_fastPHASE_for_McGregor_pop.pl output dir>

=head1 DESCRIPTION

After the wrapper script assign_BreedOfOrigin_from_fastPHASE_for_McGregor_pop.pl
has been run on the output from recode_plinks_fastPHASE_and_run.pl, for every
haplotype guess from fastPHASE an output file with possible breed of origin and
parent of origin information is assigned for each animal haplotype. This ouput
is used as input to this script to produce chromosome files for each animal, each
animal has 2 haplotypes, one predicted as coming from Sire and another from Dam.
Each haplotype shows a strech of N's or A's or -'s depending upon predicted
breed of origin information.

=head1 OPTIONS

combine_mcgregor_haps.pl takes the following arguments:

=over 4

=item -h or --help (Optional)

  Displays this helpful message.

=item -q or --quiet (Optional)

  Providing this option suppresses the log messages to the shell.
  Default: disabled

=item -m or --mcgregor-pop (Required)

  Path to McGregor population file containing, Sire, Dam and Sex information for
  each animal.
  ( mcgregor_pop_animal_ids.txt ).

=item -b or --breed (Required)

  Path to result directory of  assign_BreedOfOrigin_from_fastPHASE_for_McGregor_pop.pl
  script.

=item -o or --output (Required)

  Path to output directory .i.e where should the results from this script be
  stored?

=item -c or --chr-range (Optional)

  If mentioned, (ex: 27..30), chromosome files for that range will be created.

=back

=head1 AUTHOR

Kranti Konganti, E<lt>konganti@tamu.eduE<gt>.

=head1 COPYRIGHT

This program is distributed under the Artistic License.

=head1 DATE

Oct-05-2012

=cut
