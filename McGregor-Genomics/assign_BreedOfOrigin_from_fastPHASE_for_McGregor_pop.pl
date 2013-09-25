#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

my ($LASTCHANGEDDATE) = q$LastChangedDate: 2013-09-16 08:33:47 -0500 (Mon, 16 Sep 2013) $ =~ m/.+?\:(.+)/;
my ($VERSION) = q$LastChangedRevision: 57 $ =~ m/.+?(\d+)/;
my ($CHANGEDBY) = q$LastChangedBy: konganti $ =~ m/.+?\:(.+)/;
my $AUTHOR = 'Kranti Konganti';

# Being extra cautious for nice die message instead of Perl's message

check_and_load_modules(['Pod::Usage', 'Getopt::Long', 'File::Basename']);

# Print script info

this_script_info();

# Declare initial global variables

my ($quiet, $fastphase, $output, $help, $generation_file, $map, $family, $pedigree,
    $boa, $poa, $allow_mismatches, $is_fa, $fastphase_in_dir, $fastphase_out_dir,
    @fastphase_out_files, $hapguess_fc);

my $RECODE_SCRIPT = 'recode_plinks_fastPHASE_and_run.pl';
my $BREED_NELLORE = 'N';
my $BREED_ANGUS = 'A';
my $PARENT_SIRE = 'S';
my $PARENT_DAM = 'D';
my $RECOMBINANT = '';

my $is_valid_option = GetOptions ('help|?' => \$help,
                                  'quiet|?' => \$quiet,
                                  'fastphase=s' => \$fastphase,
                                  'output=s' => \$output,
                                  'generation=s' => \$generation_file,
                                  'family=s' => \$family,
                                  'map=s' => \$map,
                                  'pedigree=s' => \$pedigree,
                                  'allow-mismatches=i' => \$allow_mismatches,
				  'is-fa' => \$is_fa
                                  );

# Check for the validity of options
verify_options($is_valid_option, $help);
verify_input_files([$generation_file, $map, $family, $pedigree],
                   ['Generation ( animal_family_generation_sorted.list )',
                    'Map ( mcgregor_merged_umd3.map in plink format )',
                    'Family ( mcgregor_merged_umd3.fam in plink format )',
                    'Pedigree ( dam_sire_pedigree_list.v2.txt )']);

$output = validate_create_path($output, 'create', 'Output');
$fastphase = validate_create_path($fastphase, 'do not create', 'Fastphase results');

execute_system_command(0,
                       "\nOutput files will be stored at $output ...\n");

# Check for fastPhase directories produced by recode_plinks_fastPHASE_and_run
if (!defined ($is_fa)) {
    $fastphase_in_dir = $fastphase . 'fastPHASE_corrected_input/';
    $fastphase_out_dir = $fastphase . 'fastPHASE_output/';
}
else {
    $fastphase_out_dir = $fastphase;
}

execute_system_command(0,
                       "\nChecking for $fastphase_in_dir ...\n")
    if (!defined ($is_fa));
confess error("Directory ( $fastphase_in_dir ) does not exist.\n
              Did you ran $RECODE_SCRIPT?")
    if (!defined ($is_fa) && !-d $fastphase_in_dir);

execute_system_command(0,
                       "\nChecking for $fastphase_out_dir ...\n")
    if (!defined ($is_fa));
confess error("Directory ( $fastphase_out_dir ) does not exist.\n
              Did you ran $RECODE_SCRIPT?")
    if (!defined ($is_fa) && !-d $fastphase_out_dir);

# Store pedigree relation and generation relation information
my $generation_rel = store_and_get_generation();
my $family_rel = store_and_get_family([$pedigree, $family]);

# Look for fastPHASE input and output files and parse results
opendir (FASTPHASE_OUT, $fastphase_out_dir) ||
    confess error("Cannot open ( $fastphase_out_dir ): $!");

if (!defined ($is_fa)) {
    @fastphase_out_files = grep {/hapguess\_switch\.out$/} readdir FASTPHASE_OUT;
}
else {
    @fastphase_out_files = grep {/\.seq\.txt$/} readdir FASTPHASE_OUT;
}
close FASTPHASE_OUT;

$|++ if (!defined ($quiet));
if (!defined ($is_fa)) {
    $hapguess_fc = `ls $fastphase_out_dir*_hapguess_switch.out | wc -l`;
}
else {
    $hapguess_fc = `ls $fastphase_out_dir*.seq.txt | wc -l`;
}
chomp $hapguess_fc;
my $parsed_fc = 0;

foreach my $fastphase_out_file (@fastphase_out_files) {
    my ($fastphase_in_file, $marker_info_from_map_file,
        $breedOfOrigin_filename, $breedOfOrigin_fh, $haplotype,
        $marker_loc_string);

     $boa = {};
     $poa = {};

    if (!defined ($is_fa) && $fastphase_out_file =~ m/^(.+?)(\_)(chr.+?)(formatted)(\_hapguess\_switch\.)out$/) {
        $fastphase_in_file = $fastphase_in_dir . $1 . $2 . $3 . $4 . '.inp';
        $fastphase_out_file = $fastphase_out_dir . $fastphase_out_file;
        $breedOfOrigin_filename = $output . $1 . $2 . $3 . 'bo.txt';
        $breedOfOrigin_fh = open_file($breedOfOrigin_filename, '>');

        execute_system_command(0,
                               "\nParsing $map file ...\n");
        $marker_info_from_map_file = get_marker_info($3);
    }
    elsif (defined ($is_fa) && $fastphase_out_file =~ m/^(.+?)\.seq\.txt$/) {
	$fastphase_out_file = $fastphase_out_dir . $fastphase_out_file;
	$breedOfOrigin_filename = $output . $1 . '.bo.txt';
	$breedOfOrigin_fh = open_file($breedOfOrigin_filename, '>');
    }

    execute_system_command(0,
                           "\nParsing $fastphase_out_file ...\n");

    if ($fastphase_out_file eq '.' || $fastphase_out_file eq '..') {
        next;
    }

    execute_system_command(0,
                               "\nParsing $fastphase_in_file file ...\n")
	if (!defined ($is_fa));

    my $hapguess = store_hapguess($fastphase_out_file);

    $haplotype = store_boa_for_founders($hapguess, $haplotype);
    $haplotype = store_boa_poa_for_animals_with_parent_info($hapguess, $haplotype,
                                                            '^F1');
    $haplotype = store_boa_for_f2_and_later($hapguess, $haplotype);
    $haplotype = store_boa_for_hs($hapguess, $haplotype);

    # Print header for file for later parsing of marker locations and to
    # build fasta with N's and A's to view breakpoints.
    my $marker_locs;
    if (!defined ($is_fa)) {
	$marker_locs = store_marker_info_from_fastphase_input($fastphase_in_file);
	foreach my $marker_loc (@$marker_locs) {
	    $marker_loc_string .= "$marker_loc ";
	}
	chop $marker_loc_string; # Removing last space
    }

    print $breedOfOrigin_fh "#Sire_ID\tDam_ID\tGeneration\tSex(1=M|2=F)\t",
    "Animal_ID\t$marker_loc_string\tHapBreedOrigin\tHapParentOrigin\n"
	if (!defined ($is_fa));

    foreach my $animal_id (keys %$hapguess) {
        my ($sire, $dam, $sex, $generation) = get_animal_info($hapguess, $animal_id);
	my $wait = '';
	my $hap_count = 0;

        foreach my $hap (keys %{$haplotype->{$animal_id}}) {
	    $hap_count++;
	    #if (!exists $boa->{$animal_id}->{$hap} || !exists $poa->{$animal_id}->{$hap}) {
	    #print $boa->{$animal_id}->{$hap}, "\t$animal_id\n";#}
            my ($printable_hap_id) = ($hap =~ m/(HAP0|HAP1)/);
	    my $printable_hap = $hap;
            $printable_hap =~ s/HAP(0|1)\://;

	    if ($hap =~ m/HAP1/ && $hap_count == 1) {
		$wait = "$sire\t$dam\t$generation\t$sex\t$animal_id\t$hap\t$boa->{$animal_id}->{$hap}\t$poa->{$animal_id}->{$hap}\n";
		next;
	    }

	    if (!defined ($is_fa) && $hap_count == 2 && $wait ne '') {
		print $breedOfOrigin_fh "$sire\t$dam\t$generation\t$sex\t",
		    "$animal_id\t$hap\t$boa->{$animal_id}->{$hap}\t$poa->{$animal_id}->{$hap}\n",
		    $wait;
	    }
	    else {
		print $breedOfOrigin_fh "$sire\t$dam\t$generation\t$sex\t",
            "$animal_id\t$hap\t$boa->{$animal_id}->{$hap}\t$poa->{$animal_id}->{$hap}\n";
	    }

	    print $breedOfOrigin_fh ">$animal_id | $printable_hap_id\n$printable_hap\t",
	    "$boa->{$animal_id}->{$hap}\t$poa->{$animal_id}->{$hap}\n\n"
		if (defined ($is_fa));
        }
	print $breedOfOrigin_fh "\n\n" if (defined ($is_fa));
    }

    close $breedOfOrigin_fh;
    $parsed_fc++;
    print "\r$parsed_fc / $hapguess_fc files parsed";
    #print "$breedOfOrigin_filename\n";
    #exit;
}

print "\n\nDone!\n\n";





####################### Functions ###########################

# This function assign breed of origin and parent of origin information
# for McGregor Half sibs (HS) generation.
sub store_boa_for_hs {
    my $hapguess = shift;
    my $haplotype = shift;

    foreach my $animal_id (keys %$hapguess) {
        my ($sire, $dam, $sex, $generation) = get_animal_info($hapguess, $animal_id);
        my $animal_haps = get_animal_haplotypes($hapguess, $animal_id);
        my $sire_haps = get_animal_haplotypes($hapguess, $sire);
        my $dam_haps = get_animal_haplotypes($hapguess, $dam);
        my $is_match_sire;

        if ($generation =~ m/^HS/i && exists $boa->{$sire} && !exists $boa->{$dam}) {
            for (my $i=0; $i < scalar(@$animal_haps); $i++) {
                foreach my $hap (keys %{$boa->{$sire}}) {
                    my $hap_dna = $hap;
                    $hap_dna =~ s/HAP(0|1)\://;
		    if ($allow_mismatches) {
			$is_match_sire = mismatch_accept($hap_dna, @$animal_haps[$i]);
		    }
		    else {
			$is_match_sire = string_compare($hap_dna, @$animal_haps[$i]);
		    }
                    if ($is_match_sire) {
                        $haplotype->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = 1;
                        $boa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = $boa->{$sire}->{$hap};
                        $poa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = $PARENT_SIRE;

                        my $animal_other_hap_index = get_other_index($animal_haps, $i);

                        $haplotype->{$animal_id}->{"HAP$animal_other_hap_index:@$animal_haps[$animal_other_hap_index]"} = 1;
                        $boa->{$animal_id}->{"HAP$animal_other_hap_index:@$animal_haps[$animal_other_hap_index]"} = '';
                        $poa->{$animal_id}->{"HAP$animal_other_hap_index:@$animal_haps[$animal_other_hap_index]"} = '';
                        last;
                    }
                }
		last if ($is_match_sire);
            }
        }
        elsif ($generation =~ m/^HS/i && exists $boa->{$sire} && exists $boa->{$dam}) {
            confess error("Half-sib ( $animal_id ) has a Dam, which has breed of origin information
                          Sire: $sire, Dam: $dam. This is not currently handled.");
        }

	if ($generation =~ m/^HS/i && !$is_match_sire) {
	    for (my $i=0; $i < scalar(@$animal_haps); $i++) {
		$haplotype->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = 1;
		$boa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = '';
		$poa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = '';
	    }
	}
    }
    return $haplotype;
}


# This function stores and returns hash reference for breed of origin,
# parent of origin and haplotype info for F2 generation of McGregor genomics
sub store_boa_for_f2_and_later {
    my $hapguess = shift;
    my $haplotype = shift;

    foreach my $animal_id (keys %$hapguess) {
        my ($sire, $dam, $sex, $generation) = get_animal_info($hapguess, $animal_id);
        if ($generation =~ m/^NILL/i && $sire =~ m/^NILL/ && $dam =~/^NILL/) {
            print warning("Animal: $animal_id does not have neither Sire, Dam or Generation information.\nFiles checked: $generation_file, $pedigree, $family")
		if (!$quiet);
		my $animal_haps = get_animal_haplotypes($hapguess, $animal_id);
		for (my $i=0; $i < scalar(@$animal_haps); $i++) {
		    $haplotype->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = 1;
		    $boa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = '';
		    $poa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = '';
		}
        }

        if ($generation =~ m/^F2/i && exists $boa->{$sire} && exists $boa->{$dam}) {
            $haplotype = search_parent_haps_for_f2_and_later($hapguess,
                                                             $animal_id,
                                                             $haplotype);
        }
        elsif ($generation =~ m/^NILL/i && $sire !~ m/^NILL/ && $dam !~/^NILL/) {
            if (!exists $boa->{$dam}) {
		my $dam_haps = get_animal_haplotypes($hapguess, $dam);
		assign_empty_boa_poa_for_both_haps($dam_haps, $dam, $boa);
                #$boa->{$dam}->{"HAP0:UNAVAILABLE"} = '';
                #$boa->{$dam}->{"HAP1:UNAVAILABLE"} = '';
            }
            if (!exists $boa->{$sire}) {
		my $sire_haps = get_animal_haplotypes($hapguess, $sire);
		assign_empty_boa_poa_for_both_haps($sire_haps, $sire, $boa);
            }
	    if (!exists $poa->{$dam}) {
		my $dam_haps = get_animal_haplotypes($hapguess, $dam);
		assign_empty_boa_poa_for_both_haps($dam_haps, $dam, $poa);
            }
            if (!exists $poa->{$sire}) {
		my $sire_haps = get_animal_haplotypes($hapguess, $sire);
		assign_empty_boa_poa_for_both_haps($sire_haps, $sire, $poa);
            }
            $haplotype = search_parent_haps_for_f2_and_later($hapguess,
                                                             $animal_id,
                                                             $haplotype);
        }
    }
    return $haplotype;
}

# This function searches for exact match of animal's haplotype in parent's
# haplotypes and assigns breed and sire/dam information for cycle 1 F2
sub search_parent_haps_for_f2_and_later {
    my $hapguess = shift;
    my $animal_id = shift;
    my $haplotype = shift;

    my ($sire, $dam, $sex, $generation) = get_animal_info($hapguess, $animal_id);
    my $animal_haps = get_animal_haplotypes($hapguess, $animal_id);
    my $sire_haps = get_animal_haplotypes($hapguess, $sire);
    my $dam_haps = get_animal_haplotypes($hapguess, $dam);
    my ($is_match_sire, $is_match_dam, $is_match_sire_index, $is_match_dam_index);


    # Search in sire haps. If a match is found, search for match in dam with
    # animal's other haplotype string;

    for (my $i=0; $i < scalar(@$animal_haps); $i++) {
        ($is_match_sire, $is_match_sire_index) = split/\|/, string_compare($sire_haps, @$animal_haps[$i]);
        ($is_match_dam, $is_match_dam_index) = split/\|/, string_compare($dam_haps, @$animal_haps[$i]);

        if ($is_match_sire && $is_match_dam) {
            next if ($i < $#$animal_haps);

	    my $animal_hap_other_index = get_other_index($animal_haps, $i);
	    my $sire_hap_other_index = get_other_index($sire_haps, $is_match_sire_index);
	    my $dam_hap_other_index = get_other_index($dam_haps, $is_match_dam_index);

	    if ((@$animal_haps[$i] ne @$animal_haps[$animal_hap_other_index])
		&& ($boa->{$sire}->{"HAP$is_match_sire_index:$is_match_sire"} eq $boa->{$dam}->{"HAP$is_match_dam_index:$is_match_dam"})
		&& ($boa->{$sire}->{"HAP$sire_hap_other_index:@$sire_haps[$sire_hap_other_index]"} eq $boa->{$dam}->{"HAP$dam_hap_other_index:@$dam_haps[$dam_hap_other_index]"})
		&& (@$sire_haps[$is_match_sire_index] eq @$dam_haps[$is_match_dam_index])
		&& (@$sire_haps[$sire_hap_other_index] eq @$dam_haps[$dam_hap_other_index])
		&& (@$animal_haps[$i] eq @$sire_haps[$is_match_sire_index])
		&& (@$animal_haps[$i] eq @$dam_haps[$is_match_dam_index])
		&& (@$animal_haps[$animal_hap_other_index] eq @$sire_haps[$sire_hap_other_index])
		&& (@$animal_haps[$animal_hap_other_index] eq @$dam_haps[$dam_hap_other_index])) {

		$haplotype->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = 1;
		$boa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = $boa->{$dam}->{"HAP$is_match_sire_index:$is_match_sire"};
		$poa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = '';

		$haplotype->{$animal_id}->{"HAP$animal_hap_other_index:@$animal_haps[$animal_hap_other_index]"} = 1;
		$boa->{$animal_id}->{"HAP$animal_hap_other_index:@$animal_haps[$animal_hap_other_index]"} = $boa->{$sire}->{"HAP$sire_hap_other_index:@$sire_haps[$sire_hap_other_index]"};
		$poa->{$animal_id}->{"HAP$animal_hap_other_index:@$animal_haps[$animal_hap_other_index]"} = '';
	    }
	    elsif ((@$animal_haps[$i] eq @$animal_haps[$animal_hap_other_index])
		   && ($boa->{$sire}->{"HAP$is_match_sire_index:$is_match_sire"} eq $boa->{$dam}->{"HAP$is_match_dam_index:$is_match_dam"})
		   && ($boa->{$sire}->{"HAP$sire_hap_other_index:@$sire_haps[$sire_hap_other_index]"} eq $boa->{$dam}->{"HAP$dam_hap_other_index:@$dam_haps[$dam_hap_other_index]"})
		   && (@$sire_haps[$is_match_sire_index] eq @$dam_haps[$is_match_dam_index])
		   && (@$animal_haps[$animal_hap_other_index] eq @$sire_haps[$is_match_sire_index])
		   && (@$animal_haps[$animal_hap_other_index] eq @$dam_haps[$is_match_dam_index])) {

		$haplotype->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = 1;
		$boa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = $boa->{$sire}->{"HAP$is_match_sire_index:$is_match_sire"};
		$poa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = '';

		$haplotype->{$animal_id}->{"HAP$animal_hap_other_index:@$animal_haps[$animal_hap_other_index]"} = 1;
		$boa->{$animal_id}->{"HAP$animal_hap_other_index:@$animal_haps[$animal_hap_other_index]"} = $boa->{$sire}->{"HAP$is_match_sire_index:$is_match_sire"};
		$poa->{$animal_id}->{"HAP$animal_hap_other_index:@$animal_haps[$animal_hap_other_index]"} = '';
	    }
	    elsif ((@$animal_haps[$i] eq @$animal_haps[$animal_hap_other_index])
		       && ($boa->{$sire}->{"HAP$is_match_sire_index:$is_match_sire"} ne $boa->{$dam}->{"HAP$is_match_dam_index:$is_match_dam"})
		       && ($boa->{$sire}->{"HAP$sire_hap_other_index:@$sire_haps[$sire_hap_other_index]"} ne $boa->{$dam}->{"HAP$dam_hap_other_index:@$dam_haps[$dam_hap_other_index]"})
		       && (@$sire_haps[$is_match_sire_index] eq @$dam_haps[$is_match_dam_index])
		       && (@$sire_haps[$sire_hap_other_index] ne @$dam_haps[$dam_hap_other_index])
		       && (@$animal_haps[$animal_hap_other_index] eq @$sire_haps[$is_match_sire_index])
		       && (@$animal_haps[$animal_hap_other_index] eq @$dam_haps[$is_match_dam_index])) {

		    $haplotype->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = 1;
		    $boa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = '';
		    $poa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = $PARENT_DAM;

		    $haplotype->{$animal_id}->{"HAP$animal_hap_other_index:@$animal_haps[$animal_hap_other_index]"} = 1;
		    $boa->{$animal_id}->{"HAP$animal_hap_other_index:@$animal_haps[$animal_hap_other_index]"} = '';
		    $poa->{$animal_id}->{"HAP$animal_hap_other_index:@$animal_haps[$animal_hap_other_index]"} = $PARENT_SIRE;
		}
	    else {
		$haplotype->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = 1;
		$boa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = '';
		$poa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = '';

		$haplotype->{$animal_id}->{"HAP$animal_hap_other_index:@$animal_haps[$animal_hap_other_index]"} = 1;
		$boa->{$animal_id}->{"HAP$animal_hap_other_index:@$animal_haps[$animal_hap_other_index]"} = '';
		$poa->{$animal_id}->{"HAP$animal_hap_other_index:@$animal_haps[$animal_hap_other_index]"} = '';
	    }
        }
        elsif ($is_match_sire && !$is_match_dam) {

            $haplotype->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = 1;
            $boa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = $boa->{$sire}->{"HAP$is_match_sire_index:$is_match_sire"};
            $poa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = $PARENT_SIRE;

	    my $not_matched_in_dam_index = get_other_index($animal_haps, $i);

	    if ($allow_mismatches) {
		my $snp_allowed_match = mismatch_accept($dam_haps, @$animal_haps[$not_matched_in_dam_index]);
		assign_boa_poa($snp_allowed_match,
			       $animal_haps,
			       $animal_id,
			       $dam,
			       $not_matched_in_dam_index,
			       $haplotype,
			       "",
			       $PARENT_DAM
			       );
	    }
	    elsif (!$allow_mismatches) {
		($is_match_dam, $is_match_dam_index) = split/\|/, string_compare($dam_haps, @$animal_haps[$not_matched_in_dam_index]);
		assign_boa_poa($is_match_dam,
			       $animal_haps,
			       $animal_id,
			       $dam,
			       $not_matched_in_dam_index,
			       $haplotype,
			       "",
			       $PARENT_DAM
			       );
	    }
	    last;
        }
	elsif ($is_match_dam && !$is_match_sire) {

	    $haplotype->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = 1;
            $boa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = $boa->{$dam}->{"HAP$is_match_dam_index:$is_match_dam"};
            $poa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = $PARENT_DAM;

	    my $not_matched_in_sire_index = get_other_index($animal_haps, $i);

	    if ($allow_mismatches) {
		my $snp_allowed_match = mismatch_accept($sire_haps, @$animal_haps[$not_matched_in_sire_index]);
		assign_boa_poa($snp_allowed_match,
				$animal_haps,
				$animal_id,
				$sire,
				$not_matched_in_sire_index,
				$haplotype,
				"",
				$PARENT_SIRE
				);
	    }
	    elsif (!$allow_mismatches) {
		($is_match_sire, $is_match_sire_index) = split/\|/, string_compare($sire_haps, @$animal_haps[$not_matched_in_sire_index]);
		assign_boa_poa($is_match_sire,
			       $animal_haps,
			       $animal_id,
			       $sire,
			       $not_matched_in_sire_index,
			       $haplotype,
			       "",
			       $PARENT_SIRE
			       );
	    }
	    last;
	}
	elsif (!$is_match_sire && !$is_match_dam) {
	    ($is_match_sire, $is_match_sire_index) = string_compare($sire_haps, @$animal_haps[0]) if (!$allow_mismatches);
	    $is_match_sire = mismatch_accept($sire_haps, @$animal_haps[0]) if ($allow_mismatches);
	    if ($is_match_sire) {
		assign_boa_poa($is_match_sire,
			       $animal_haps,
			       $animal_id,
			       $sire,
			       "0",
			       $haplotype,
			       "",
			       $PARENT_SIRE
			       );
		($is_match_dam, $is_match_dam_index) = string_compare($dam_haps, @$animal_haps[1]) if (!$allow_mismatches);
		$is_match_dam = mismatch_accept($dam_haps, @$animal_haps[1]) if ($allow_mismatches);
		assign_boa_poa($is_match_dam,
			       $animal_haps,
			       $animal_id,
			       $dam,
			       "1",
			       $haplotype,
			       "",
			       $PARENT_DAM
			       );

	    }
	    else {
		($is_match_dam, $is_match_dam_index) = string_compare($dam_haps, @$animal_haps[0]) if (!$allow_mismatches);
		$is_match_dam = mismatch_accept($dam_haps, @$animal_haps[0]) if ($allow_mismatches);
		assign_boa_poa($is_match_dam,
			       $animal_haps,
			       $animal_id,
			       $dam,
			       "0",
			       $haplotype,
			       "",
			       $PARENT_DAM
			       );
		($is_match_sire, $is_match_sire_index) = string_compare($sire_haps, @$animal_haps[1]) if (!$allow_mismatches);
		$is_match_sire = mismatch_accept($sire_haps, @$animal_haps[1]) if ($allow_mismatches);
		assign_boa_poa($is_match_sire,
			       $animal_haps,
			       $animal_id,
			       $sire,
			       "1",
			       $haplotype,
			       "",
			       $PARENT_SIRE
			       );
	    }
	    last;
	}
    }
    return $haplotype;
}

# This function automates the assignment of breed of origin and parent of
# origin for hash ref
sub assign_boa_poa {
    my $parent_dna = shift;
    my $animal_haps = shift;
    my $animal_id = shift;
    my $parent_id = shift;
    my $index = shift;
    my $haplotype = shift;
    my $nill = shift;
    my $which_parent = shift;

    if ($parent_dna && exists $boa->{$parent_id}->{"HAP0:$parent_dna"}) {
	$haplotype->{$animal_id}->{"HAP$index:@$animal_haps[$index]"} = 1;
	$boa->{$animal_id}->{"HAP$index:@$animal_haps[$index]"} = $boa->{$parent_id}->{"HAP0:$parent_dna"};
	$poa->{$animal_id}->{"HAP$index:@$animal_haps[$index]"} = $which_parent;
    }
    elsif ($parent_dna && exists $boa->{$parent_id}->{"HAP1:$parent_dna"}) {
	$haplotype->{$animal_id}->{"HAP$index:@$animal_haps[$index]"} = 1;
	$boa->{$animal_id}->{"HAP$index:@$animal_haps[$index]"} = $boa->{$parent_id}->{"HAP1:$parent_dna"};
	$poa->{$animal_id}->{"HAP$index:@$animal_haps[$index]"} = $which_parent;
    }
    else {
	$haplotype->{$animal_id}->{"HAP$index:@$animal_haps[$index]"} = 1;
	$boa->{$animal_id}->{"HAP$index:@$animal_haps[$index]"} = $nill;
	$poa->{$animal_id}->{"HAP$index:@$animal_haps[$index]"} = $nill;
	return 0;
    }
    return 1;
}

# This function is an attempt to reduce code rendundancy.
sub assign_empty_boa_poa_for_both_haps {
    my $haps = shift;
    my $indiv = shift;
    my $hash_ref = shift;

    for(my $i=0; $i < scalar(@$haps); $i++) {
	$hash_ref->{$indiv}->{"HAP$i:@$haps[$i]"} = '';
    }
    return;
}

# This function stores and returns the hash reference for breed of origin,
# parent of origin and haplotype info for founders of McGregor genomics
# population
sub store_boa_for_founders {
    my $hapguess = shift;
    my $haplotype = shift;

    # Try to assign breed of origin for the founders if they have no parent
    # information or only 1 parent information with dna.
    foreach my $animal_id (keys %$hapguess) {

        my ($sire, $dam, $sex, $generation) = get_animal_info($hapguess, $animal_id);
        my $animal_haps = get_animal_haplotypes($hapguess, $animal_id);
        my $sire_haps = get_animal_haplotypes($hapguess, $sire);
        my $dam_haps = get_animal_haplotypes($hapguess, $dam);

        if ($generation =~ m/^GP|GGP/i) {
            if (!$sire && !$dam && $sex == 1) {
                for (my $i=0; $i < scalar(@$animal_haps); $i++) {
                    $haplotype->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = 1;
                    $boa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = $BREED_NELLORE;
                    $poa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = '';
                }
            }
            elsif (!$sire && !$dam && $sex == 2) {
               for (my $i=0; $i < scalar(@$animal_haps); $i++) {
                    $haplotype->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = 1;
                    $boa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = $BREED_ANGUS;
                    $poa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = '';
                }
            }

            if ($sire && $dam && !exists $hapguess->{$dam} && !exists $hapguess->{$sire}) {
                my $animal_boa = '';
                $animal_boa = $BREED_NELLORE if ($sex == 1);
                $animal_boa = $BREED_ANGUS if ($sex == 2);

                for (my $i=0; $i < scalar(@$animal_haps); $i++) {
                    $haplotype->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = 1;
                    $boa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = $animal_boa;#$BREED_NELLORE;
                    $poa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = '';
                }

                for (my $i=0; $i < scalar(@$animal_haps); $i++) {
                    my $hap_index_equal_to_sire = string_compare($sire_haps, @$animal_haps[$i]);
                    if ($hap_index_equal_to_sire) {
                        $poa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = $PARENT_SIRE;

                        my $j = get_other_index($animal_haps, $i);

                        $poa->{$animal_id}->{"HAP$j:@$animal_haps[$j]"} = $PARENT_DAM;
                        last;
                    }
                }
            }
            elsif ($sire && (!$dam || !exists $hapguess->{$dam})) {

                for (my $i=0; $i < scalar(@$animal_haps); $i++) {
                    $haplotype->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = 1;
                    $boa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = $BREED_NELLORE;
                    $poa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = '';
                }

                for (my $i=0; $i < scalar(@$animal_haps); $i++) {
                    my $hap_index_equal_to_sire = string_compare($sire_haps, @$animal_haps[$i]);
                    if ($hap_index_equal_to_sire) {
                        $poa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = $PARENT_SIRE;

                        my $j = get_other_index($animal_haps, $i);

                        $poa->{$animal_id}->{"HAP$j:@$animal_haps[$j]"} = $PARENT_DAM;
                        last;
                    }
                }
            }
            elsif ($dam && (!$sire || !exists $hapguess->{$sire})) {

                for (my $i=0; $i < scalar(@$animal_haps); $i++) {
                    $haplotype->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = 1;
                    $boa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = $BREED_ANGUS;
                    $poa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = '';
                }

                for (my $i=0; $i < scalar(@$animal_haps); $i++) {
                    my $hap_index_equal_to_dam = string_compare($dam_haps, @$animal_haps[$i]);
                    if ($hap_index_equal_to_dam) {
                        $poa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = $PARENT_DAM;

                        my $j = get_other_index($animal_haps, $i);

                        $poa->{$animal_id}->{"HAP$j:@$animal_haps[$j]"} = $PARENT_SIRE;
                        last;
                    }
                }
            }
        }
    }

    $haplotype = store_boa_poa_for_animals_with_parent_info($hapguess, $haplotype,
                                                            '^GP|GGP');
    return $haplotype;
}

# This function stores breed of origin and parent of origin for animals
# which have DNA, Parent and breed of origin information and returns a hashref
# For F1 indivs of McGregor genomics population, since only N sire x A dam cross
# was made to obtain F1, then all of them are same breed of origin
# i.e NA (Nellore-Angus). We only need to determine parent of origin.

sub store_boa_poa_for_animals_with_parent_info {
    my $hapguess = shift;
    my $haplotype = shift;
    my $generation_to_match = shift;

    foreach my $animal_id (keys %$hapguess) {

        my ($sire, $dam, $sex, $generation) = get_animal_info($hapguess, $animal_id);

        if ($generation =~ m/$generation_to_match/i && exists $boa->{$sire} && exists $boa->{$dam}) {
            my @sire_s_boa = values %{$boa->{$sire}};
            my @dam_s_boa = values %{$boa->{$dam}};
            $haplotype = search_parent_haps_for_founders_and_f1($hapguess,
                                                                $animal_id,
                                                                $haplotype,
                                                                \@sire_s_boa,
                                                                \@dam_s_boa);
        }
        elsif ($generation =~ m/^F1/i && exists $boa->{$sire} && !exists $boa->{$dam}
               && !exists $haplotype->{$animal_id}) {
            my @sire_s_boa = values %{$boa->{$sire}};
            my @dam_s_boa = ($BREED_ANGUS, $BREED_ANGUS);
            $haplotype = search_parent_haps_for_founders_and_f1($hapguess,
                                                                $animal_id,
                                                                $haplotype,
                                                                \@sire_s_boa,
                                                                \@dam_s_boa);
        }
        elsif ($generation =~ m/^F1/i && !exists $boa->{$sire} && exists $boa->{$dam}
               && !exists $haplotype->{$animal_id}) {
            my @dam_s_boa = values %{$boa->{$dam}};
            my @sire_s_boa = ($BREED_NELLORE, $BREED_NELLORE);
            $haplotype = search_parent_haps_for_founders_and_f1($hapguess,
                                                                $animal_id,
                                                                $haplotype,
                                                                \@sire_s_boa,
                                                                \@dam_s_boa);
        }
        elsif ($generation =~ m/^F1/i && !exists $boa->{$sire} && !exists $boa->{$dam}) {
            confess error("\nCannot find parent information,  either Sire or Dam for F1 individual $animal_id\n");
        }
    }
    return $haplotype;
}

# This function searches for exact match of animal's haplotype in parent's
# haplotypes and assigns breed and sire/dam information
sub search_parent_haps_for_founders_and_f1 {
    my $hapguess = shift;
    my $animal_id = shift;
    my $haplotype = shift;
    my $sire_s_boa = shift;
    my $dam_s_boa = shift;
    my $no_match_for_both_haps = 0;

    my ($sire, $dam, $sex, $generation) = get_animal_info($hapguess, $animal_id);
    my $animal_haps = get_animal_haplotypes($hapguess, $animal_id);
    my $sire_haps = get_animal_haplotypes($hapguess, $sire);
    my $dam_haps = get_animal_haplotypes($hapguess, $dam);
    my ($is_match_sire, $is_match_dam);


    for (my $i=0; $i < scalar(@$animal_haps); $i++) {
	if ($allow_mismatches) {
	    $is_match_sire = mismatch_accept($sire_haps, @$animal_haps[$i]);
	    $is_match_dam = mismatch_accept($dam_haps, @$animal_haps[$i]);
	}
	else {
	    $is_match_sire = string_compare($sire_haps, @$animal_haps[$i]);
	    $is_match_dam = string_compare($dam_haps, @$animal_haps[$i]);
	}

        if ($is_match_sire) {
            # Since founders breed of origin is either pure Nellore or
            # Angus, we are using hard coded N or A
            $haplotype->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = 1;
            $boa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = @$sire_s_boa[0];
            $poa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = $PARENT_SIRE;

            my $j = get_other_index($animal_haps, $i);

            $haplotype->{$animal_id}->{"HAP$j:@$animal_haps[$j]"} = 1;
            $boa->{$animal_id}->{"HAP$j:@$animal_haps[$j]"} = @$dam_s_boa[0];
            $poa->{$animal_id}->{"HAP$j:@$animal_haps[$j]"} = $PARENT_DAM;
            last;
        }
        elsif ($is_match_dam) {
            # Since founders breed of origin is either pure Nellore or
            # Angus, we are using hard coded N or A
            $haplotype->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = 1;
            $boa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = @$dam_s_boa[0];
            $poa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = $PARENT_DAM;

            my $j = get_other_index($animal_haps, $i);

            $haplotype->{$animal_id}->{"HAP$j:@$animal_haps[$j]"} = 1;
            $boa->{$animal_id}->{"HAP$j:@$animal_haps[$j]"} = @$sire_s_boa[0];
            $poa->{$animal_id}->{"HAP$j:@$animal_haps[$j]"} = $PARENT_SIRE;
            last;
        }
        elsif (!$is_match_dam && !$is_match_sire) {
            $no_match_for_both_haps++;
        }
    }

    # If we cannot find a match for the animal's haplotype in both
    # father and mother, then leave blank.
    if ($no_match_for_both_haps == 2) {
        for (my $i=0; $i < scalar(@$animal_haps); $i++) {
            $haplotype->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = 1;
            $boa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = '';
            $poa->{$animal_id}->{"HAP$i:@$animal_haps[$i]"} = '';
        }
    }
    return $haplotype;
}

# This function compares strings of argument 1 to argument 2
# it only does exact matches
sub string_compare {
    my $parent_haps = shift;
    my $animal_hap = shift;

    confess error("You cannot define --allow_mismatches 0
                  This script will search for exact matches by default")
        if (defined($allow_mismatches) && $allow_mismatches == 0);

    if (ref($parent_haps) eq 'ARRAY') {
	my ($exact_match) = grep { $_ eq $animal_hap } @$parent_haps;
	my ($exact_match_index) = grep { @$parent_haps[$_] eq $animal_hap } 0..$#$parent_haps;
	return "$exact_match|$exact_match_index"
	    if ($exact_match && ($exact_match_index == 0 || $exact_match_index == 1));
    }
    elsif ($parent_haps eq $animal_hap) {
        return 1;
    }
    return 0;
}

# This function will try to allow mismatches in unassigned haplotypes
sub mismatch_accept {
    my $parent_haps = shift;
    my $animal_hap = shift;

    if (ref($parent_haps) eq 'ARRAY') {
        foreach my $parent_hap (@$parent_haps) {
            next if ($parent_hap eq  '' || $parent_hap eq 'UNAVAILABLE');
            return $parent_hap if (get_no_of_mismatches($parent_hap, $animal_hap));
        }
    }
    else {
        return 0 if ($parent_haps eq '');
        return $parent_haps if (get_no_of_mismatches($parent_haps, $animal_hap));
    }
    return 0;
}

# This function compares animal hap to parent's hap to count number of
# nucleotide mismatches;
sub get_no_of_mismatches {
    my $parent_hap_string = shift;
    my $animal_hap_string = shift;

    $parent_hap_string =~ s/\s+//g;
    $animal_hap_string =~ s/\s+//g;
    my @parent_bases = split//, $parent_hap_string;
    my @bases = split//, $animal_hap_string;
    my $no_of_mismatches = 0;

    #my @indexes;
    for(my $i=0; $i<scalar(@bases); $i++) {
        if ($bases[$i] ne $parent_bases[$i]) {
            #push(@indexes, $i);
	    $no_of_mismatches++;
	    if ($no_of_mismatches > $allow_mismatches) {
		#$RECOMBINANT = '';
		return 0;
	    }
        }
    }

    #if (is_recombinant(\@indexes)) {
#	print "@indexes\n";
#	$RECOMBINANT = 'Rec';
 #   }
  #  else {
#	$RECOMBINANT = '';
 #   }

    return 1;
}

sub is_recombinant {
    my $indexes = shift;
    my $first_val = $indexes->[0];

    #return 0 if (scalar(@$indexes) != $allow_mismatches);
    for my $index (0 .. $#$indexes) {
        return unless $indexes->[$index] == $first_val + $index;
    }
    return 1;
}

# This function returns the remaining index value of the other
# haplotype
sub get_other_index {
    my $animal_haps = shift;
    my $i = shift;
    my $j = scalar(@$animal_haps) - $i - 1;
    $j = scalar(@$animal_haps) - 1 if ($j < 0);
    return $j;
}

# This function loops through each haplotype for each animal
# Split and return animal info
sub get_animal_info {
    my $hapguess = shift;
    my $animal_id = shift;
    my ($sire, $dam, $sex, $generation);

    if (!exists $family_rel->{$animal_id}) {
        $sire = $dam = $sex = $generation = 'NILL';
    }
    else {
         ($sire, $dam, $sex) = split/\|/, $family_rel->{$animal_id};
         if (!exists $generation_rel->{$animal_id}) {
            $generation = 'NILL';
        }
        else {
            $generation = $generation_rel->{$animal_id};
        }
    }
    return ($sire, $dam, $sex, $generation);
}

# Return haplotypes
sub get_animal_haplotypes {
    my $hapguess = shift;
    my $animal_id = shift;
    my ($hap_line_1, $hap_line_2);

    if (!exists $hapguess->{$animal_id}) {
        $hap_line_1 = $hap_line_2 = '';
    }
    else {
        ($hap_line_1, $hap_line_2) = split/\|/, $hapguess->{$animal_id};
    }
    return ([$hap_line_1, $hap_line_2]);
}

# Store family information into memory and return reference to hash
sub store_and_get_family {
    my $store;
    my $files = shift;

    foreach my $file (@$files) {
        my $family_fh = open_file($file, '<');

        while (my $line = <$family_fh>) {
            chomp $line;
            next if ($line =~ m/#/);
            my @fam_attributes = split/\t/, $line;

            if (scalar(@fam_attributes) == 0) {
                confess error("There seems to be an empty line in family file
                              ( $file )");
            }
            elsif (scalar(@fam_attributes) == 6) { # plink family file
                # The generation file is very inconsistent when compared to map
                # and ped files given to me. I have to dig down deep / debug
                # so many times and I found that some animal ids have '-' and in
                # map file the ids have '/'. For ex: animal id 403-5 is present in
                # generation file where as the same animal id is represented as
                # 403/5 in ped file. This drives the code nuts in not being
                # able to assign sire and dam information.
                $fam_attributes[1] =~ s/\W/\-/g;
                $fam_attributes[2] =~ s/\W/\-/g;
                $fam_attributes[3] =~ s/\W/\-/g;

                $fam_attributes[1] = uc($fam_attributes[1]); # animal id
                $fam_attributes[2] = uc($fam_attributes[2]); # sire id
                $fam_attributes[3] = uc($fam_attributes[3]); # dam id
                $store->{$fam_attributes[1]} = "$fam_attributes[2]|$fam_attributes[3]|$fam_attributes[4]"
                    if (!exists $store->{$fam_attributes[1]});
            }
            elsif (scalar(@fam_attributes) == 4) { # dam_sire_pedigree.list.v2 file
                $fam_attributes[0] =~ s/\W/\-/g;
                $fam_attributes[1] =~ s/\W/\-/g;
                $fam_attributes[2] =~ s/\W/\-/g;

                $fam_attributes[0] = uc($fam_attributes[0]); # animal id
                $fam_attributes[1] = uc($fam_attributes[1]); # sire id
                $fam_attributes[2] = uc($fam_attributes[2]); # dam id
                $store->{$fam_attributes[0]} = "$fam_attributes[1]|$fam_attributes[2]|$fam_attributes[3]"
                    if (!exists $store->{$fam_attributes[0]});
            }
            else {
                confess error("Could not get all columns from line $line
                              in family file ( $file )");
            }
        }
        close $family_fh;
    }
    return $store;
}

# Store generation information of each animal and return reference to hash
sub store_and_get_generation {
    my $store;
    my $generation_fh = open_file($generation_file, '<');

    while (my $line = <$generation_fh>) {
        chomp $line;
        if ($line =~ m/^(\S+)\,(\S+)\,(\S+)$/) {
            my $animal_id_alias = uc($1);
            my $animal_id = uc($2);
            my $generation = uc($3);

            # The generation file is very inconsistent when compared to map
            # and ped files given to me. I have to dig down deep / debug
            # so many times and I found that some animal ids have '-' and in
            # map file the ids have '/'. For ex: animal id 403-5 is present in
            # generation file where as the same animal id is represented as
            # 403/5 in ped file. This drives the code nuts in not being
            # able to assign sire and dam information.
            $animal_id_alias =~ s/\W/\-/g;
            $animal_id =~ s/\W/\-/g;
            $store->{$animal_id_alias} = $generation;
            $store->{$animal_id} = $generation;
        }
        else {
            confess error("Could not get all the values from the line\n
                          Is the line well formatted?\n$line");
        }
    }
    close $generation_fh;
    return $store;
}

# Store marker info (P line) from fastPHASE input
sub store_marker_info_from_fastphase_input {
    my $file = shift;
    my $found_marker_string = 0;
    my @marker_locs;

    my $fastphase_inp_fh = open_file($file, '<');

    while (my $line = <$fastphase_inp_fh>) {
        chomp $line;
        if ($line =~ m/^P\s+(\d+.+)/) {
            @marker_locs = split/\s+/, $1;
            $found_marker_string = 1;
            last;
        }
    }

    confess error("Cannot find marker positions in ( $file ) i.e P string of
                  fastPHASE input file.")
        if (!$found_marker_string);

    close $fastphase_inp_fh;
    return \@marker_locs;
}

# Store hapguess information from fastPHASE
sub store_hapguess {
    my $file = shift;
    my ($animal_id, $hapguess);

    if (defined ($is_fa)) {
	check_and_load_modules(['Bio::SeqIO']);
	my $fa_file_obj = Bio::SeqIO->new('-file' => $file,
					  '-format' => 'fasta');
	while (my $seq_in_hap1 = $fa_file_obj->next_seq) {
	    my $seq_id = $seq_in_hap1->id;
	    $seq_id = uc($seq_id);
	    my $seq_in_hap2 = $fa_file_obj->next_seq;
	    $hapguess->{$seq_id} = $seq_in_hap1->seq . '|' . $seq_in_hap2->seq;
	}
    }
    else {
	my $fastphase_res_fh = open_file($file, '<');
	while (my $line = <$fastphase_res_fh>) {
	    chomp $line;
	    next if ($line =~ m/^\s+/);

	    if ($line =~ m/^\#\s+ID\s+(.+)/) {
		$animal_id = uc($1);

		# The generation file is very inconsistent when compared to map
		# and ped files given to me. I have to dig down deep / debug
		# so many times and I found that some animal ids have '-' and in
		# map file the ids have '/'. For ex: animal id 403-5 is present in
		# generation file where as the same animal id is represented as
		# 403/5 in ped file. This drives the code nuts in not being
		# able to assign sire and dam information.
		$animal_id =~ s/\W/\-/g;

		my $hap_line_1 = <$fastphase_res_fh>;
		$hap_line_1 = strip_leading_and_trailing_spaces($hap_line_1);

		my $hap_line_2 = <$fastphase_res_fh>;
		$hap_line_2 = strip_leading_and_trailing_spaces($hap_line_2);

		$hapguess->{$animal_id} = "$hap_line_1|$hap_line_2";
	    }
	}
	close $fastphase_res_fh;
    }
    return $hapguess;
}

# Store marker info from map file and return reference to hash
sub get_marker_info {
    my $chr = shift;
    my $store;
    my $map_fh = open_file ($map, '<');

    while (my $line = <$map_fh>) {
        chomp $line;
        if ($line =~ m/^\w+\s+\S+\s+\S+\s+\d+$/) {
            my ($chr_no, $snp_id, $gen_distance, $chr_pos) = split /\s+/, $line;

            $chr = lc($chr);
            $chr_no = lc($chr_no);

            $chr_no = 'chr'.$chr_no if ($chr_no !~ m/^chr/);
            $store->{$chr_pos} = $snp_id if ($chr eq $chr_no);
        }
        else {
            confess error("Could not get all the values from the line\n
                          Is the line file well formatted?\n$line");
        }
    }
    close $map_fh;
    return $store;
}

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

    foreach my $module (@$module_list) {
        my $module_installed = eval ("use $module; 1");

        confess error("Required module not installed: $module\n
                      The module and its dependencies must be installed at system level"),
            if (!$module_installed);
    }
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

# Shell msg that differentiates log from error
sub error {
    my $msg = shift;
    print "\nERROR!\n------\n$msg\n\n";
    return;
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

# Print this script's info
sub this_script_info {
    print "\n", '@', '-' x 80, '@', "\n";
    print "  Program Name       :  " , basename($0), "\n";
    print "  Version            :  $VERSION\n";
    print "  Author             :  $AUTHOR\n";
    print "  Last Changed By    : $CHANGEDBY\n";
    print "  Last Changed Date  : $LASTCHANGEDDATE\n";
    print '@', '-' x 80, '@', "\n\n";
    return;
}

__END__

=head1 NAME

assign_BreedOfOrigin_from_fastPHASE_for_McGregor_pop.pl

=head1 SYNOPSIS

This script will attempt to assign breed of origin and parent of origin for
output files from fastPHASE which were run on output files from the
recode_plinks_fastPHASE_and_run.pl script.

Examples:

    perl assign_BreedOfOrigin_from_fastPHASE_for_McGregor_pop.pl -h

    assign_BreedOfOrigin_from_fastPHASE_for_McGregor_pop.pl -m ../mcgregor_merged_umd3.map -g ../animal_family_generation_sorted.list.txt -fam ../mcgregor_merged_umd3.fam -fas . -p ../dam_sire_pedigree.list.v2.txt -o breed -q

=head1 DESCRIPTION

After the wrapper script recode_plinks_fastPHASE_and_run.pl has been run to
extract specified chunks of markers, fastPHASE is run on those input files
and the output of the fastPHASE is stored in a directory as mentioned to the
recode_plinks_fastPHASE_and_run.pl script. This script will take that fastPHASE
output directory and will attempt to assign breed of origin and parent of origin
information to each of the output files in the directory mentioned to this script
with --fastphase option. This code works for only McGregor population, since
assiging breed of origin has been broken down depending upon establishing
breed of origin for Founders, F1, HS, F2 and later.

=head1 OPTIONS

assign_BreedOfOrigin_from_fastPHASE_for_McGregor_pop.pl takes the following arguments:

=over 4

=item -h or --help (Optional)

  Displays this helpful message.

=item -q or --quiet (Optional)

  Providing this option suppresses the log messages to the shell.
  Default: disabled

=item -m or --map (Required)

  Path to MAP file in plink format for McGregor population
  ( mcgregor_merged_umd3.map ).

=item -g or --generation (Required)

  Path to generation file describing, sire - dam - generation information for
  each animal of McGregor population ( animal_family_generation_sorted.list ).

=item -fam or --family (Required)

  Path to Family file in plink format for McGregor population
  ( mcgregor_merged_umd3.fam ).

=item -p or --pedigree (Required)

  Path to pedigree file of McGregor population
  ( dam_sire_pedigree_list.v2.txt ).

=item -fas or --fastphase (Required)

  Path to fastPHASE output directory which was a result of running the
  wrapper script recode_plinks_fastPHASE_and_run.pl.

=item -is or --is-fa (Optional)

  If the phased files are in fasta format, single fasta per haplotype, then this
  option will adjust the program accordingly.

=item -a or --allow-mismatches (Optional)

  -a 3 allows 3 mismatches when comparing haplotype of an animal to its
  parents' haplotype.
  Default: disabled (which means exact match search).

=item -o or --output (Required)

  Name of the output directory.

=back

=head1 AUTHOR

Kranti Konganti, E<lt>konganti@tamu.eduE<gt>.

=head1 COPYRIGHT

This program is distributed under the Artistic License.

=head1 DATE

May-05-2012

=cut
