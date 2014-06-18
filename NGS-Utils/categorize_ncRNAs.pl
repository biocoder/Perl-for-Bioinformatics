#!/usr/bin/env perl

use strict;
use Carp;
use Getopt::Long qw(:config pass_through);
use IO::Routine;
use Set::IntervalTree;
use Parallel::ForkManager;

my ($LASTCHANGEDBY) = q$LastChangedBy: konganti $ =~ m/.+?\:(.+)/;
my ($LASTCHANGEDDATE) = q$LastChangedDate: 2013-10-09 12:46:11 -0500 (Wed, 09 Oct 2013) $ =~ m/.+?\:(.+)/;
my ($VERSION) = q$LastChangedRevision: 64 $ =~ m/.+?(\d+)/;
my $AUTHORFULLNAME = 'Kranti Konganti';

my ($help, $quiet, $cuffcmp, $genePred, $out, $sample_names,
    $fpkm_cutoff, $cov_cutoff, $refGenePred, $length, $categorize,
    $min_exons, $overlap, $novel, $extract_pat, $no_tmp,
    $antisense_only, $disp_anti_option, $gtf_bin, $num_cpu,
    $linc_rna_prox);

my ($p_file_names_gtf, $p_file_names_txt) = [];
my $ncRNA_class = {};

my $is_valid_option = GetOptions('help|?'              => \$help,
				 'quiet'               => \$quiet,
				 'cuffcmp=s'           => \$cuffcmp,
				 'annotation=s'        => \$refGenePred,
				 'out=s'               => \$out,
				 'sample-names=s'      => \$sample_names,
				 'fpkm-cutoff=f'       => \$fpkm_cutoff,
				 'cov-cutoff=f'        => \$cov_cutoff, 
				 'genePred'            => \$genePred,
				 'categorize'          => \$categorize,
				 'length=i'            => \$length,
				 'min-exons=i'         => \$min_exons,
				 'overlap=f'           => \$overlap,
				 'include-novel'       => \$novel,
				 'clean-tmp'           => \$no_tmp,
				 'bin-gtfToGenePred=s' => \$gtf_bin,
				 'antisense-only'      => \$antisense_only,
				 'cpu=i'               => \$num_cpu,
				 'linc-rna-prox=i'     => \$linc_rna_prox);

my $io = IO::Routine->new($help, $quiet);
my $s_time = $io->start_timer;

$io->verify_options([$is_valid_option, $sample_names, 
		     $refGenePred, $out, $cuffcmp]);

$io->this_script_info($io->file_basename($0),
		      $VERSION,
		      $AUTHORFULLNAME,
		      $LASTCHANGEDBY,
		      $LASTCHANGEDDATE);

$io->c_time('Analysis started...');
$io->c_time('Verifying options...');

# Define Defaults
$fpkm_cutoff = 0.0 if (!defined $fpkm_cutoff || $fpkm_cutoff eq '');
$cov_cutoff = 0.0 if (!defined $cov_cutoff || $cov_cutoff eq '');
$length = 200 if (!defined $length || $length eq '');
$overlap = 0 if (!defined $overlap || $overlap eq '');
$min_exons = 2 if (!defined $min_exons || $min_exons eq '');

$io->c_time('Validating output path...');
my $output = $io->validate_create_path($out, 'create', 
				       'Output directory');

$io->c_time('Checking for required GNU core utils...');
$io->check_sys_level_cmds(['grep', 'sed'], 
			  ['2.6.3', '4.2.1']);

$io->c_time('Checking for annotation file in Gene Prediction format...');
$io->verify_files([$refGenePred], ['Gene Prediction']);
my $refGenePred_fh = $io->open_file('<', $refGenePred);

$io->c_time('Checking for Cuffcompare tracking and Cufflinks assembled transcript files...');
my $cuffcmp_fh = $io->open_file('<', $cuffcmp);

$io->error('Cufflinks assembled transcript files not provided.')
    if ($#ARGV < 0);

$io->c_time('Tainting sample names...');
$sample_names = $io->strip_leading_and_trailing_spaces($sample_names);
$sample_names =~ s/\s+/\_/g;
$sample_names =~ s/\,$//;
my @lables = split/\,/, $sample_names;

for (0 .. $#ARGV) {
    $io->verify_files([$ARGV[$_]],
                      ["Cufflinks assembled transcript"]);
    push @{$p_file_names_gtf}, $output . $io->file_basename($ARGV[$_]) . '.' . $lables[$_] . '.putative_ncRNAs.gtf';
    push @{$p_file_names_txt}, $output . $io->file_basename($ARGV[$_]) . '.' . $lables[$_] . '.putative_ncRNAs.txt';
    unlink $p_file_names_gtf->[$_] if (-e $p_file_names_gtf->[$_] && !defined($genePred) && !defined($categorize));
}

$disp_anti_option = 'True' if (defined($antisense_only));
$disp_anti_option = 'False' if (!defined($antisense_only));

if (defined($categorize)) {    
    $io->execute_system_command(0,
				"Using options:\n--------------\n" .
				"Minimum transcript length              : $length\n" .
				"Minimum exon overlap percentage        : $overlap\n" .
				"Minimum number of exons per transcript : $min_exons\n" .
				"Extract only Antisense exon overlaps   : $disp_anti_option");
    class_ncRNAs();
}
elsif (defined($genePred)) {
    $io->execute_system_command(0,
                                "Using options:\n--------------\n" .
                                "Minimum transcript length              : $length\n" .
                                "Minimum exon overlap percentage        : $overlap\n" .
				"Minimum number of exons per transcript : $min_exons\n" . 
				"Extract only Antisense exon overlaps   : $disp_anti_option");
    get_genePred();
    class_ncRNAs();
}
else {
    $io->execute_system_command(0,
				"Using options:\n--------------\n" .
				"FPKM cutoff                            : $fpkm_cutoff\n" .
				"Coverage cutoff                        : $cov_cutoff\n".
				"Minimum transcript length              : $length\n" .
				"Minimum exon overlap percentage        : $overlap\n" .
				"Minimum number of exons per transcript : $min_exons\n" . 
				"Extract only Antisense exon overlaps   : $disp_anti_option");
    get_putative_ncRNAs();
    get_genePred();
    class_ncRNAs();
}

if (defined($no_tmp)) {
    $io->c_time('Removing intermediate files...');
    for (0 .. $#ARGV) {
	unlink $p_file_names_gtf->[$_];
	unlink $p_file_names_txt->[$_];
    }
}

$io->c_time('categorize_ncRNAs Finished!');
$io->end_timer($s_time);
exit;



############################# Functions ###############################

# Get putative list of ncRNAs from Cuffcompare tracking file corresponding
# to cufflinks assembled transcript fragments.
sub get_putative_ncRNAs {
    $io->c_time('Getting putative list of ncRNAs in GTF format...');
    if (defined $novel) {
	$extract_pat = 'j|i|o|u|x';
    }
    else {
	$extract_pat = 'i|o|u|x';
    }
    
    while (my $line = <$cuffcmp_fh>) {
	chomp $line;
	$line = $io->strip_leading_and_trailing_spaces($line);
	my ($t_id, $loc_id, $loc_name, $class_code, @cols) = split/\t/, $line;
	if ($class_code =~ m/$extract_pat/i) {
	    for (0 .. $#cols) {
		my ($q_loc_id, $q_t_id, $discard) = split/\|/, $cols[$_];
		next if (!$q_t_id || $q_t_id eq '');
		$q_t_id =~ s/\./\\./g;
		my $t_lines = $io->execute_get_sys_cmd_output("grep -iP \'$q_t_id\' $ARGV[$_]", 0);
		if ($t_lines =~ m/.+?[FR]PKM.+?\"(.+?)\".+?cov.+?\"(.+?)\".+/i) {
		    $io->execute_system_command("grep -iP \'$q_t_id\' $ARGV[$_] | sed -e \'s\/\$\/ class_code \"$class_code\"\;\/' >> $p_file_names_gtf->[$_]") if ($1 >= $fpkm_cutoff && $2 >= $cov_cutoff);
		}
		else {
		    $io->execute_system_command("grep -iP \'$q_t_id\' $ARGV[$_] | sed -e \'s\/\$\/ class_code \"$class_code\"\;\/' >> $p_file_names_gtf->[$_]");
		}
	    }
	}
    }
    return;
}

# Convert GTF to gene prediction format.
sub get_genePred {
    my $cpu = '';
    $io->c_time('Converting putative ncRNAs list to Gene Prediction format using gtfToGenePred tool');

    my $exe_gtfToGenePred = 'gtfToGenePred';
    $exe_gtfToGenePred = $gtf_bin if (defined($gtf_bin) && $gtf_bin ne '');

    my $check_for_gtfToGenePred = $io->execute_get_sys_cmd_output($exe_gtfToGenePred, 0);
    
    $io->error('Cannot find gtfToGenePred tool in your path') 
	if ($check_for_gtfToGenePred !~ m/.*?gtfToGenePred.*?convert a GTF file to a genePred/i);

    if (defined $num_cpu) {
        $cpu = Parallel::ForkManager->new($num_cpu);
	$cpu->set_max_procs($num_cpu);
    }
    
    for (0 .. $#$p_file_names_gtf) {
	$cpu->start and next if (defined $num_cpu);
	$io->verify_files([$p_file_names_gtf->[$_]], ['GTF']);
	$io->execute_system_command("$exe_gtfToGenePred -genePredExt -geneNameAsName2 $p_file_names_gtf->[$_] $p_file_names_txt->[$_] 2> /dev/null",
	    "$exe_gtfToGenePred -genePredExt -geneNameAsName2 $p_file_names_gtf->[$_] $p_file_names_txt->[$_] 2> /dev/null");
	$cpu->finish if (defined $num_cpu);
    }
    $cpu->wait_all_children;
    return;
}

# Categorize ncRNAs
sub class_ncRNAs {
    my $refAnnot = store_coords($refGenePred);
    my $cpu = '';

    if (defined $num_cpu) {
	$cpu = Parallel::ForkManager->new($num_cpu);
	$cpu->set_max_procs($num_cpu);
    }
  

    for (0 .. $#ARGV) {

	$cpu->start and next if (defined $num_cpu);

        my $p_gtf = $p_file_names_gtf->[$_];
        my $p_ncRNAs = store_coords($p_file_names_txt->[$_]);
        my $c_ncRNAs = $output . $io->file_basename($ARGV[$_]) . '.' . $lables[$_] . '.putative.class.ncRNAs.gtf';
	my $u_ncRNAs = $output . $io->file_basename($ARGV[$_]) . '.' . $lables[$_] . '.putative.noClass.ncRNAs.gtf';
        unlink $c_ncRNAs if (-e $c_ncRNAs);
	unlink $u_ncRNAs if (-e $u_ncRNAs);
	
	my ($num_ex_ov, $num_incs, $num_concs, $num_poncs, $num_lincs, $noclass, $num_noSense,
	    $discard) = 0;

	$io->c_time('Categorizing ncRNAs (Exonic overlaps) [' . $io->file_basename($p_file_names_gtf->[$_], 'suffix') . ' ]...');
	($num_ex_ov, $num_noSense) = calc_overlaps('exonic', $p_gtf, $p_ncRNAs, $c_ncRNAs, $refAnnot, $u_ncRNAs);

	$io->c_time('Categorizing ncRNAs (Intronic overlaps - incs) [' . $io->file_basename($p_file_names_gtf->[$_], 'suffix') . ' ]...');
	($num_incs, $discard) = calc_overlaps('Inc', $p_gtf, $p_ncRNAs, $c_ncRNAs, $refAnnot);

	$io->c_time('Categorizing ncRNAs (Intronic overlaps - concs) [' . $io->file_basename($p_file_names_gtf->[$_], 'suffix') . ' ]...');
	($num_concs, $discard) = calc_overlaps('Conc', $p_gtf, $p_ncRNAs, $c_ncRNAs, $refAnnot);

	$io->c_time('Categorizing ncRNAs (Intronic overlaps - poncs) [' . $io->file_basename($p_file_names_gtf->[$_], 'suffix') . ' ]...');
        ($num_poncs, $discard) = calc_overlaps('Ponc', $p_gtf, $p_ncRNAs, $c_ncRNAs, $refAnnot);

	$io->c_time('Categorizing ncRNAs (lincRNA) [' . $io->file_basename($p_file_names_gtf->[$_], 'suffix') . ' ]...');
        ($num_lincs, $noclass) = calc_lincRNAs($p_gtf, $p_ncRNAs, $c_ncRNAs, $refAnnot, $u_ncRNAs);
	
	$io->c_time("\n\nncRNA Summary [" . $io->file_basename($p_file_names_gtf->[$_], 'suffix') . " ] :\n" . 
		    "----------------------------------------------------------------------\n" .
		    "LincRNAs: $num_lincs\n" . 
		    "Intronic overlaps - concs: $num_concs\n" .
                    "Intronic overlaps - poncs: $num_poncs\n" . 
		    "Intronic overlaps - incs: $num_incs\n" . 
		    "Exonic overlaps: $num_ex_ov\n" . 
		    "Total Categorized: " . ($num_lincs + 
					     $num_concs + 
					     $num_incs + 
					     $num_ex_ov +
					     $num_poncs) . 
		    "\nUncategorized: " . 
		    ($noclass + $num_noSense) . "\n" );
	
	$cpu->finish if (defined $num_cpu);
    }

    $cpu->wait_all_children if (defined $num_cpu);
    return;
}

# Calculate LincRNAs
sub calc_lincRNAs {
    my $p_gtf = shift;
    my $p_ncRNAs = shift;
    my $c_ncRNAs = shift;
    my $refAnnot = shift;
    my $u_ncRNAs = shift;
    my $found = 0;
    my $num_noClass = 0;
    my $ov_found_chk_flag = 0;

    foreach my $chr (keys %{$refAnnot}) {
	my $ref_gene_tree = Set::IntervalTree->new();
	   
	foreach my $ref_gene (values @{$refAnnot->{$chr}}) {
	    my ($ref_strand,
		$ref_tr_start,
		$ref_tr_end,
		$ref_exons,
		$ref_exon_starts,
		$ref_exon_ends,
		$ref_tr_id) = get_parts($ref_gene);
	    
	    if (defined $linc_rna_prox && $linc_rna_prox > 0) {
		$ref_tr_start += $linc_rna_prox;
		$ref_tr_end += $linc_rna_prox;
		$ov_found_chk_flag = 1;
	    }

	    $ref_gene_tree->insert($ref_tr_id, $ref_tr_start, $ref_tr_end);
	} 

	
	foreach my $transfrag (values @{$p_ncRNAs->{$chr}}) {
	    my ($nc_strand,
		$nc_tr_start,
		$nc_tr_end,
		$nc_exons,
		$nc_exon_starts,
		$nc_exon_ends,
		$nc_tr_id) = get_parts($transfrag);
	 
	    my $unique_key = "$nc_tr_id$nc_strand$nc_tr_start$nc_tr_end$nc_exons" .
                @$nc_exon_starts . @$nc_exon_ends;
	    my $ncRNA_length = $nc_tr_end - $nc_tr_start;
	    my $ncRNA_length = $nc_tr_end - $nc_tr_start;
	    my $ov_found = $ref_gene_tree->fetch($nc_tr_start, $nc_tr_end);
	    
	    if (scalar(@$ov_found) == $ov_found_chk_flag &&
		!exists $ncRNA_class->{$unique_key} &&
		$ncRNA_length >= $length) {
		$found++;
		$ncRNA_class->{$unique_key} = 1;
		$io->execute_system_command("grep -iP \'$nc_tr_id\' $p_gtf | sed -e \'s\/\$\/ transcript_length \"$ncRNA_length\"\; ncRNA_type \"LincRNA\";\/\' >> $c_ncRNAs", 0);
	    }
	    elsif (!exists $ncRNA_class->{$unique_key}) {
		$num_noClass++;
		$ncRNA_class->{$unique_key} = 1;
		$io->execute_system_command("grep -iP \'$nc_tr_id\' $p_gtf | sed -e \'s\/\$\/ transcript_length \"$ncRNA_length\"\; ncRNA_type \"No Class \(\?\)\"\;\/\' >> $u_ncRNAs", 0);
	    }
	}
    }
    return ($found, $num_noClass);
}

# Calcualte overlaps of putative ncRNAs with reference information.
sub calc_overlaps {
    my $mode = shift;
    my $p_gtf = shift;
    my $p_ncRNAs = shift;
    my $c_ncRNAs = shift;
    my $refAnnot = shift;
    my $u_ncRNAs = shift;
    my $found = 0;
    my $num_noSense = 0;

    # Calculate overlap.
    foreach my $nc_chr (keys %{$p_ncRNAs}) {

	my $nc_int_tree = Set::IntervalTree->new() if ($mode =~ m/^ponc|conc$/i);
	    
	for (my $ncRNA_line=0; $ncRNA_line <= $#{$p_ncRNAs->{$nc_chr}}; $ncRNA_line++) {
	    	  
	    my ($nc_strand,
		$nc_tr_start,
		$nc_tr_end,
		$nc_exons,
		$nc_exon_starts,
		$nc_exon_ends, 
		$nc_tr_id) = get_parts($p_ncRNAs->{$nc_chr}->[$ncRNA_line]);
	    
	    my $unique_key = "$nc_tr_id$nc_strand$nc_tr_start$nc_tr_end$nc_exons" .
		@$nc_exon_starts . @$nc_exon_ends;

	    my $ncRNA_length = $nc_tr_end - $nc_tr_start;
	    $nc_int_tree->insert($nc_tr_id, $nc_tr_start, $nc_tr_end) if ($mode =~ m/^ponc|conc$/i);

	    foreach my $ref_gene (values @{$refAnnot->{$nc_chr}}) {
		
		my ($ref_strand, 
		    $ref_tr_start, 
		    $ref_tr_end, 
		    $ref_exons, 
		    $ref_exon_starts, 
		    $ref_exon_ends,
		    $ref_tr_id) = get_parts($ref_gene);
		
		my $is_strand_Antisense = is_Antisense($ref_strand, $nc_strand);

		# Only lncRNAs, adjustable by user with -length and -min-exons parameters.
		
		if ($ncRNA_length >= $length &&
		    $nc_exons >= $min_exons) { 

		    # Complete overlap reference gene with ncRNA intron.
		    if ($mode =~ m/^conc$/i &&
                        $nc_tr_start < $ref_tr_start &&
                        $nc_tr_end > $ref_tr_end) {

			my $ov_tr_found = $nc_int_tree->fetch($ref_tr_start ,$ref_tr_end);
			my $is_ncRNA_Conc = 0;
			$is_ncRNA_Conc = 1 if (scalar(@$ov_tr_found) >= 1);
			my $is_ncRNA_Conc = is_intronicOverlap($ref_tr_start, $ref_tr_end, $nc_exon_starts, $nc_exon_ends);
			
			if ($is_ncRNA_Conc &&
			    $is_strand_Antisense &&
			    !exists $ncRNA_class->{$unique_key}) {
			    $ncRNA_class->{$unique_key} = "ncRNA_type \"Antisense intronic overlap (Conc) with $ref_tr_id\";";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line--;
			    $found++;
			    last;
			}
			elsif ($is_ncRNA_Conc &&
			       !$is_strand_Antisense &&
			       !exists $ncRNA_class->{$unique_key}) {
			    $ncRNA_class->{$unique_key} = "ncRNA_type \"Sense intronic overlap (Conc) with $ref_tr_id\";";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
                            $ncRNA_line--;
			    $found++;
			    last;
			}
			elsif (!$is_ncRNA_Conc &&
			    !exists $ncRNA_class->{$unique_key}) {
			    $ncRNA_class->{$unique_key} = "ncRNA_type \"Intronic overlap (Conc)\";";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line--;
			    $found++;
			    last;
			}
		    }
		    		    
		    # Complete overlap of ncRNA within reference intron	
		    if ($mode =~ m/^inc$/i &&
			$nc_tr_start >= $ref_tr_start &&
                        $nc_tr_end <= $ref_tr_end) {
			
			my $is_ncRNA_Inc = is_intronicOverlap($nc_tr_start, $nc_tr_end, $ref_exon_starts, $ref_exon_ends);
			
			if ($is_ncRNA_Inc &&
			    $is_strand_Antisense &&
			    !exists $ncRNA_class->{$unique_key}) {
			    $ncRNA_class->{$unique_key} = "ncRNA_type \"Antisense intronic overlap (Inc) with $ref_tr_id\";";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line--;
			    $found++;
			    last;
			}
			elsif ($is_ncRNA_Inc &&
			       !$is_strand_Antisense &&
			       !exists $ncRNA_class->{$unique_key}) {
			    $ncRNA_class->{$unique_key} = "ncRNA_type \"Sense intronic overlap (Inc) with $ref_tr_id\";";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line--;
			    $found++;
			    last;
			}
			elsif ($is_ncRNA_Inc &&
			       !exists $ncRNA_class->{$unique_key}) {
			    $ncRNA_class->{$unique_key} = "ncRNA_type \"Intronic overlap (Inc) with $ref_tr_id\";";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line--;
			    $found++;
			    last;
			}
		    }

		    # Antisense and Sense partial reference intronic overlap.
		    if ($mode =~ m/^ponc$/i &&
			(( ($nc_tr_start > $ref_tr_start) &&
			   ($nc_tr_end > $ref_tr_end) 
			 ) || 
			 ( ($nc_tr_start < $ref_tr_start) &&
			   ($nc_tr_end < $ref_tr_end)
			 ))
			) {
			 
			my $found_intron_ov = $nc_int_tree->fetch($ref_tr_start, $ref_tr_end);
			my $is_ncRNA_exonicOverlap = is_exonicOverlap($ref_exon_starts, $ref_exon_ends, $nc_exon_starts, $nc_exon_ends);
			my $retain_overlap = $overlap;
			$overlap = 0;
			
			if (!$is_ncRNA_exonicOverlap &&
			    $is_strand_Antisense &&
			    !exists $ncRNA_class->{$unique_key} &&
			    scalar(@$found_intron_ov >= 1)) {
			    $ncRNA_class->{$unique_key} = "ncRNA_type \"Antisense partial intronic overlap (Ponc) with $ref_tr_id\";";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line--;
			    $found++;
			    last;
			}
			elsif (!$is_ncRNA_exonicOverlap &&
			       !$is_strand_Antisense &&
			       !exists $ncRNA_class->{$unique_key} &&
			       scalar(@$found_intron_ov >= 1)) {
			    $ncRNA_class->{$unique_key} = "ncRNA_type \"Sense partial intronic overlap (Ponc) with $ref_tr_id\";";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line--;
			    $found++;
			    last;
			}
			elsif (!$is_ncRNA_exonicOverlap &&
			       !exists $ncRNA_class->{$unique_key} &&
			       scalar(@$found_intron_ov >= 1)) {
			    $ncRNA_class->{$unique_key} = "ncRNA_type \"Partial intronic overlap (Ponc) with $ref_tr_id\";";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line--;
			    $found++;
			    last;
			}
			$overlap = $retain_overlap;
		    }
		    
		    # Antisense and Sense reference exonic overlap.
		    if ($mode =~ m/^exonic$/i) {
			
			my $is_ncRNA_exonicOverlap = is_exonicOverlap($ref_exon_starts, $ref_exon_ends, $nc_exon_starts, $nc_exon_ends);

			if ($is_ncRNA_exonicOverlap &&
			    $is_strand_Antisense &&
			    !exists $ncRNA_class->{$unique_key}) {
			    $ncRNA_class->{$unique_key} = "ncRNA_type \"Antisense exonic overlap with $ref_tr_id\";";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line--;
			    $found++;
			    last;
			}
			elsif ($is_ncRNA_exonicOverlap &&
			       !$is_strand_Antisense &&
			       !exists $ncRNA_class->{$unique_key}) {
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line--; 
			    $num_noSense++,
			    $io->execute_system_command("grep -iP \'$nc_tr_id\' $p_gtf | sed -e \'s\/\$\/ transcript_length \"$ncRNA_length\"\; ncRNA_type \"No Class \(\?\)\"\;\/\' >> $u_ncRNAs", 0),
			    last if (defined($antisense_only) && $antisense_only);
			    $found++;
			    $ncRNA_class->{$unique_key} = "ncRNA_type \"Sense exonic overlap with $ref_tr_id\";";
			    last;
			}
			elsif ($is_ncRNA_exonicOverlap &&
			       !exists $ncRNA_class->{$unique_key}) {
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line--;
			    $num_noSense++,
			    $io->execute_system_command("grep -iP \'$nc_tr_id\' $p_gtf | sed -e \'s\/\$\/ transcript_length \"$ncRNA_length\"\; ncRNA_type \"No Class \(\?\)\"\;\/\' >> $u_ncRNAs", 0),
			    $ncRNA_class->{$unique_key} = "ncRNA_type \"No Class \(\?\)\"\;",
			    last if (defined($antisense_only) && $antisense_only);
			    $found++;
			    $ncRNA_class->{$unique_key} = "ncRNA_type \"Exonic overlap with $ref_tr_id\";";
			    last;
			}
		    }
		}
	    }
	    $io->execute_system_command("grep -iP \'$nc_tr_id\' $p_gtf | sed -e \'s\/\$\/ transcript_length \"$ncRNA_length\"\; $ncRNA_class->{$unique_key}\/\' >> $c_ncRNAs", 0)
		if ($ncRNA_class->{$unique_key} ne '');
	}
    }
    return ($found, $num_noSense);
}

# Store genome coordinates.
sub store_coords {
    my $f = shift;
    my $fh = $io->open_file('<', $f);
    my $store = {};

    $io->c_time('Reading information from gene prediction format file [ ' . 
		$io->file_basename($f, 'suffix') . ' ] ...');
    while (my $line = <$fh>) {
	chomp $line;
	$line = $io->strip_leading_and_trailing_spaces($line);
	my ($t_id, $chr, $strand, $tr_start, $tr_end,
	    $cds_start, $cds_end, $num_exons, $exon_starts, $exon_ends, @rem) = split/\t/, $line;
	$exon_starts =~ s/\,$//;
	$exon_ends =~ s/\,$//;

	$chr = 'chr$chr' if ($line =~ m/^ens/i && $strand =~ m/^\+|\-|\.$/ && $num_exons =~ m/\d+/);

	$io->error('Supplied file [ ' . $io->file_basename($f, 'suffix') . ' ] does not seem to be in gene prediction format...' .
		   "\n\nError occured on line:\n\n$line\n")
	    if ($chr !~ m/^(chr|ens|uc)/i || $strand !~ m/^\+|\-|\.$/ || $num_exons !~ m/\d+/);
	
	push @{$store->{lc($chr)}}, "$strand|$tr_start|$tr_end|$cds_start|$cds_end|$num_exons|$exon_starts|$exon_ends|$t_id";
    }
    return $store;
}

# Split and return columns.
sub get_parts {
    my @line_parts = split /\|/, shift;
    
    return ($line_parts[0],
	    $line_parts[1],
	    $line_parts[2],
	    $line_parts[5],
	    [split /\,/, $line_parts[6]],
	    [split /\,/, $line_parts[7]],
	    $line_parts[8]);
}

# Calculate exonic overlap.
sub is_exonicOverlap {
    my ($s_ex_st, $s_ex_end, $c_ex_st, $c_ex_end) = @_;
    
    for (0 .. $#$c_ex_st) {

	my $nc_ex_coord = $_;
	my $c_ex_len = $c_ex_end->[$nc_ex_coord] - $c_ex_st->[$nc_ex_coord];
	
	for (0 .. $#$s_ex_st) {

	    my $s_ex_len = $s_ex_end->[$_] - $s_ex_st->[$_];
	    my $ex_ov_per = ($c_ex_len / $s_ex_len) * 100;
	    
	    if ($c_ex_st->[$nc_ex_coord] <= $s_ex_st->[$_] && 
		$c_ex_st->[$nc_ex_coord] <= $s_ex_end->[$_] && 
		$c_ex_end->[$nc_ex_coord] >= $s_ex_st->[$_] && 
		$c_ex_end->[$nc_ex_coord] <= $s_ex_end->[$_]) {
		return 1 if ( defined($overlap) && ($ex_ov_per >= $overlap) );
	    }
	    elsif ($c_ex_st->[$nc_ex_coord] >= $s_ex_st->[$_] &&
		   $c_ex_st->[$nc_ex_coord] <= $s_ex_end->[$_] &&
		   $c_ex_end->[$nc_ex_coord] >= $s_ex_st->[$_] &&
		   $c_ex_end->[$nc_ex_coord] <= $s_ex_end->[$_]) {
		return 1 if ( defined($overlap) && ($ex_ov_per >= $overlap) );
	    }
	    elsif ($c_ex_st->[$nc_ex_coord] >= $s_ex_st->[$_] &&
		   $c_ex_st->[$nc_ex_coord] <= $s_ex_end->[$_] &&
		   $c_ex_end->[$nc_ex_coord] >= $s_ex_st->[$_] &&
		   $c_ex_end->[$nc_ex_coord] >= $s_ex_end->[$_]) {
		return 1 if ( defined($overlap) && ($ex_ov_per >= $overlap) );
	    }
	    elsif ($c_ex_st->[$nc_ex_coord] <= $s_ex_st->[$_] &&
		   $c_ex_st->[$nc_ex_coord] <= $s_ex_end->[$_] &&
		   $c_ex_end->[$nc_ex_coord] >= $s_ex_st->[$_] && 
		   $c_ex_end->[$nc_ex_coord] >= $s_ex_end->[$_]) {
		return 1 if ( defined($overlap) && ($ex_ov_per >= $overlap) );
	    }
	}
    }
    return 0;
}

# Calculate ncRNA intronic overlap with reference gene information, both Concs and Poncs.
sub is_intronicOverlap {
    my ($tr_start, $tr_end, 
	$ex_st, $ex_end) = @_;
    
    # Return false if it is only single exon transcript.
    return 0 if ($#$ex_st == 0);

    for (1 .. $#$ex_st) {
	my $prev_exon_end = $_ - 1;
	return 1 if ($tr_end < $ex_st->[$_] &&
		     $tr_start > $ex_end->[$prev_exon_end]);
	
	#return 1 if ( ($tr_end - $tr_start) < ($ex_st->[$_] - $ex_end->[$prev_exon_end]) )
    }
    return 0;
}

# Calculate strandedness.
sub is_Antisense {
    my ($ref_strand, $comp_strand) = @_;

    if ($ref_strand ne $comp_strand &&
	$ref_strand =~ m/^\+|\-$/ &&
	$comp_strand =~ m/^\+|\-$/) {
	return 1;
    } 
    return 0;
}

__END__

=head1 NAME

categorize_ncRNAs.pl

=head1 SYNOPSIS

This script will categorize the Cufflinks' assembled transcripts into ncRNA classes as mentioned 
in the paper: http://genome.cshlp.org/content/22/3/577.full. 

=head2 DOCUMENTATION

    perldoc categorize_ncRNAs.pl

=head3 EXAMPLES:

    perl categorize_ncRNAs.pl -min-exons 1 -fpkm 2 -length 200 -annotation refSeq_ucscKnown_ensemble.txt -sample-names 'M0,M1,M2' -out lncRNA -cuffcmp cuffcmp.tracking m0/transcripts.gtf m1/transcripts.gtf m2/cufflinks/transcripts.gtf

=head1 DESCRIPTION

Cufflinks includes a program called Cuffcompare, which will compare the transcripts generated 
across all the samples with provided reference annotation of choice and generates class code 
for each of the assembled transcripts, while tracking the transcripts at the respective loci in 
each sample. This script will take in the Cuffcompare tracking file and the corresponding assembled 
transcripts for each sample in GTF format and produces a putative list of novel ncRNAs and categorizes
them as mentioned in the paper (http://genome.cshlp.org/content/22/3/577.full). Best results can be 
obtained when all the known annotation resources are combined. For organisms like Human and Mouse, 
combining RefSeq, UCSC Known genes and Ensembl genes can help filter a lot of known protein-coding 
genes and known ncRNAs. When Cufflinks' assembled transript files are supplied in GTF format, they
are first converted to Gene Prediction format using gtfToGenePred tool from UCSC, which can be downloaded
from http://hgdownload.cse.ucsc.edu/admin/exe. The tool is flexible and can resume in parts as mentioned
in the OPTIONS' section.

=head3 OUTPUT:

The final output is written to *.putative.class.gtf files and depends upon number of Cufflinks assembled transcript files
supplied.

=head2 KNOWN ISSUES

The annotation should be in the Gene Prediction format. When Gene Prediction format files for RefSeq
are downloaded from UCSC, the first column is generally a "bin" id unique to SQL table from UCSC. This
column *MUST* be removed prior to using it as an annotation file. This can be done using the "cut"
command.

Ex: cut -f 2- refGene.txt > refGene.nobin.txt

=head1 OPTIONS

categorize_ncRNAs.pl takes the following arguments:

=over 4

=item -h or --help (Optional)

    Displays this helpful message.

=item -q or --quiet (Optional)

    Turn off logging.

=item -cuff or --cuffcmp (Required)

    Path to Cuffcompare tracking file.
    Ex: cuffcmp.tracking

=item -annot or --annotation (Required)

    Path to annotation file in Gene Prediction format.
    Ex: refGene.txt

=item --out (Required)

    Path to output directory.

=item -sample or --sample-names (Required)

    Sample names in order of supplied transcripts files.
    Ex: --sample-names 'Sample1,Sample2,Sample3';

=item -fpkm or --fpkm-cutoff (Optional)

    Default: disabled

    Extract transcript features whose FPKM / RPKM value is at least this much.
    This can be a floating point value.

=item -cov or --cov-cutoff (Optional)

    Default: disabled

    Extract transript features whose coverage is at least this much. This
    can be a floating point value.

=item -inc or --include-novel (Optional)

    Default: disabled

    By default, transcript features belonging to class codes "x", "o", "i" and "u"
    are extracted. Providing this option also extracts Cufflinks transcripts 
    classified as novel ("j") isoforms.

=item -gene or -genePred (Optional)

    Default: disabled

    The script first extracts transcripts that belong to the Cufflinks'
    class codes "x", "o", "i" and "u" or class codes "x", "o", "i", "u"
    and "j" and generates the respective *putative_ncRNAs.gtf files.
    If for any reason the script fails on moving forward, it can be 
    asked to skip the extraction step and resume from converting the 
    *putative_ncRNAs.gtf files to Gene Prediction format and then 
    the categorization step with this option.

=item -cat or --categorize (Optional)

    Default: disabled

    Providing this option skips the extraction and convertion steps and
    continues the pipeline from categorizing putative ncRNAs.

=item -len or --length (Optional)

    Default: 200

    Extract transcripts whose length is at least this much.


=item -min or --min-exons (Optional)

    Default: 2

    Extract transcripts which contain at least this many number of exons.

=item -ov or --overlap (Optional)

    Default: disabled

    When calculating exonic overlaps with reference exon boundaries, consider
    it as an exonic overlap if the Cufflinks assembled transcripts' exon overlaps
    with reference exon by at least this much percentage. This can be floating 
    point value.

=item -clean or --clean-tmp (Optional)

    Default: disabled

    Remove intermediate files. Specifically, *putative_ncRNAs.gtf and
    *putative_ncRNAs.txt files are removed.

=item -anti or --antisense (Optional)

    Default: disabled
    
    When reporting exonic overlaps with reference exons, report only Antisense
    exonic overlaps.

=item --linc-rna-prox (Optional)

    Default: disabled

    When reporting Long intergenic ncRNAs, report only those which are within this
    many number of bases upstream or downstream of the reference gene.

=item -cpu or --cpu (Optional)

    Default: 1
    
    Use this many number of CPUs to run in parallel.

=back

=head1 AUTHOR

Kranti Konganti, E<lt>konganti@tamu.eduE<gt>.

=head1 COPYRIGHT

This program is distributed under the Artistic License.

=head1 DATE

Oct-09-2013

=cut
