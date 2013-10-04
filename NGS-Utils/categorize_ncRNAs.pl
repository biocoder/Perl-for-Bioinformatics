#!/usr/bin/env perl

use strict;
use Carp;
use Getopt::Long qw(:config pass_through);
use IO::Routine;
use File::Basename;
use Set::IntervalTree;

my ($LASTCHANGEDBY) = q$LastChangedBy: konganti $ =~ m/.+?\:(.+)/;
my ($LASTCHANGEDDATE) = q$LastChangedDate: 2013-10-04 09:08:44 -0500 (Fri, 04 Oct 2013) $ =~ m/.+?\:(.+)/;
my ($VERSION) = q$LastChangedRevision: 62 $ =~ m/.+?(\d+)/;
my $AUTHORFULLNAME = 'Kranti Konganti';

my $io = IO::Routine->new();
my $s_time = $io->start_timer();

my ($help, $quiet, $cuffcmp, $genePred, $out, $sample_names,
    $fpkm_cutoff, $cov_cutoff, $refGenePred, $length, $categorize,
    $min_exons, $overlap, $novel, $extract_pat);
my ($p_file_names_gtf, $p_file_names_txt) = [];
my $ncRNA_class = {};

my $is_valid_option = GetOptions('help|?'         => \$help,
				 'quiet'          => \$quiet,
				 'cuffcmp=s'      => \$cuffcmp,
				 'annotation=s'   => \$refGenePred,
				 'out=s'          => \$out,
				 'sample-names=s' => \$sample_names,
				 'fpkm-cutoff=f'  => \$fpkm_cutoff,
				 'cov-cutoff=f'   => \$cov_cutoff, 
				 'genePred'       => \$genePred,
				 'categorize'     => \$categorize,
				 'length=i'       => \$length,
				 'min-exons=i'    => \$min_exons,
				 'overlap=i'      => \$overlap,
				 'novel'          => \$novel);


$io->verify_options([$is_valid_option, $sample_names, 
		     $refGenePred, $out],
                    $help);

$io->this_script_info($0, $VERSION, $AUTHORFULLNAME, $LASTCHANGEDBY, $LASTCHANGEDDATE, $quiet);

$io->c_time('Analysis started...', $quiet);
$io->c_time('Verifying options...', $quiet);

# Define defaults
$fpkm_cutoff = 0.0 if (!defined $fpkm_cutoff || $fpkm_cutoff eq '');
$cov_cutoff = 0.0 if (!defined $fpkm_cutoff || $fpkm_cutoff eq '');
$length = 200 if (!defined $length || $length eq '');
$overlap = 0 if (!defined $overlap || $overlap eq '');
$min_exons = 2 if (!defined $min_exons || $min_exons eq '');

$io->c_time('Validating output path...', $quiet);
my $output = $io->validate_create_path($out, 'create', 
				       'Output directory');

$io->c_time('Checking for required GNU core utils...', $quiet);
$io->check_sys_level_cmds(['grep', 'sed'], 
			  ['2.6.3', '4.2.1']);

$io->c_time('Checking for annotation file in Gene Prediction format...');
$io->verify_files([$refGenePred], ['Gene Prediction']);
my $refGenePred_fh = $io->open_file('<', $refGenePred);

$io->c_time('Checking for Cuffcompare tracking and Cufflinks assembled transcript files...');
my $file_basename = sub {
    my @file_p_attrs = fileparse(shift, qr/\.[^.]*/);
    return "$file_p_attrs[0]$file_p_attrs[2]" if (shift eq 'suffix');
    return $file_p_attrs[0]
};
my $cuffcmp_fh = $io->open_file('<', $cuffcmp);

$io->error('Cufflinks assembled transcript files not provided.')
    if ($#ARGV < 0);

$io->c_time('Tainting sample names...', $quiet);
$sample_names = $io->strip_leading_and_trailing_spaces($sample_names);
$sample_names =~ s/\s+/\_/g;
my @lables = split/\,/, $sample_names;

for (0 .. $#ARGV) {
    $io->verify_files([$ARGV[$_]],
                      ["Cufflinks assembled transcript"]);
    push @{$p_file_names_gtf}, $output . $file_basename->($ARGV[$_]) . '.' . $lables[$_] . '.putative_ncRNAs.gtf';
    push @{$p_file_names_txt}, $output . $file_basename->($ARGV[$_]) . '.' . $lables[$_] . '.putative_ncRNAs.txt';
    unlink $p_file_names_gtf->[$_] if (-e $p_file_names_gtf->[$_] && !defined($genePred) && !defined($categorize));
}

if (defined($categorize)) {
    class_ncRNAs();
}
elsif (defined($genePred)) {
    get_genePred();
    class_ncRNAs();
}
else {
    get_putative_ncRNAs();
    get_genePred();
    class_ncRNAs();
}

$io->c_time('Done!', $quiet);
print $io->end_timer($s_time, $quiet), "\n";
exit;



############################# Functions ###############################



# Get putative list of ncRNAs from Cuffcompare tracking file corresponding
# to cufflinks assembled transcript fragments.
sub get_putative_ncRNAs {
    $io->c_time('Getting putative list of ncRNAs in GTF format...', $quiet);
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
		if ($t_lines =~ m/.+?FPKM.+?\"(.+?)\".+?cov.+?\"(.+?)\".+/i) {
		    $io->execute_system_command("grep -iP \'$q_t_id\' $ARGV[$_] | sed -e \'s\/\$\/ Class code \"$class_code\"\;\/' >> $p_file_names_gtf->[$_]")if ($1 >= $fpkm_cutoff && $2 >= $cov_cutoff);
		}	
	    }
	}
    }
    return;
}

# Convert GTF to gene prediction format.
sub get_genePred {
    $io->c_time('Converting putative ncRNAs list to Gene Prediction format using gtfToGenePred tool', $quiet);
    my $check_for_gtfToGenePred = $io->execute_get_sys_cmd_output('gtfToGenePred', 0);
    
    $io->error('Cannot find gtfToGenePred tool in your path') 
	if ($check_for_gtfToGenePred !~ m/.*?gtfToGenePred.*?convert a GTF file to a genePred/i);
    
    for (0 .. $#$p_file_names_gtf) {
	$io->verify_files([$p_file_names_gtf->[$_]], ['GTF']);
	$io->execute_system_command("gtfToGenePred -genePredExt -geneNameAsName2 $p_file_names_gtf->[$_] $p_file_names_txt->[$_]",
	    "gtfToGenePred -genePredExt -geneNameAsName2 $p_file_names_gtf->[$_] $p_file_names_txt->[$_]");
    }
    return;
}

# Categorize ncRNAs
sub class_ncRNAs {
    my $refAnnot = store_coords($refGenePred);

    for (0 .. $#ARGV) {
        my $p_gtf = $p_file_names_gtf->[$_];
        my $p_ncRNAs = store_coords($p_file_names_txt->[$_]);
        my $c_ncRNAs = $output . $file_basename->($ARGV[$_]) . '.' . $lables[$_] . '.class.ncRNAs.gtf';
        unlink $c_ncRNAs if (-e $c_ncRNAs);
	
	my ($num_ex_ov, $num_incs, $num_concs, $num_poncs, $num_lincs, $noclass) = 0;

	$io->c_time('Categorizing ncRNAs (Exonic overlaps)...', $quiet);
	my $num_ex_ov = calc_overlaps('exonic', $p_gtf, $p_ncRNAs, $c_ncRNAs, $refAnnot);

	$io->c_time('Categorizing ncRNAs (Intronic overlaps - incs)...', $quiet);
	$num_incs = calc_overlaps('Inc', $p_gtf, $p_ncRNAs, $c_ncRNAs, $refAnnot);

	$io->c_time('Categorizing ncRNAs (Intronic overlaps - concs)...', $quiet);
	$num_concs = calc_overlaps('Conc', $p_gtf, $p_ncRNAs, $c_ncRNAs, $refAnnot);

	$io->c_time('Categorizing ncRNAs (Intronic overlaps - poncs)...', $quiet);
        $num_poncs = calc_overlaps('Ponc', $p_gtf, $p_ncRNAs, $c_ncRNAs, $refAnnot);

	$io->c_time('Categorizing ncRNAs (lincRNA)...', $quiet);
        ($num_lincs, $noclass) = calc_lincRNAs($p_gtf, $p_ncRNAs, $c_ncRNAs, $refAnnot);
	
	$io->c_time("\n\nncRNA Summary:\n" . 
		    "------------------\n" .
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
		    "\nUncategorized: $noclass\n" );
    }
	return;
}

# Calculate LincRNAs
sub calc_lincRNAs {
    my $p_gtf = shift;
    my $p_ncRNAs = shift;
    my $c_ncRNAs = shift;
    my $refAnnot = shift;
    my $found = 0;
    my $num_noClass = 0;

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
	    
	    my $ov_found = $ref_gene_tree->fetch($nc_tr_start, $nc_tr_end);
	    
	    if (scalar(@$ov_found) == 0) {
		$found++;
		    $io->execute_system_command("grep -iP \'$nc_tr_id\' $p_gtf | sed -e \'s\/\$\/ LincRNA;\/\' >> $c_ncRNAs", 0)
	    }
	    else {
		$num_noClass++;
		$io->execute_system_command("grep -iP \'$nc_tr_id\' $p_gtf | sed -e \'s\/\$\/ No Class \(\?\);\/\' >> $c_ncRNAs", 0)
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
    my $found = 0;

    # Calculate overlap.
    foreach my $nc_chr (keys %{$p_ncRNAs}) {

	my $nc_int_tree = Set::IntervalTree->new() if ($mode =~ m/^ponc|conc$/i);
	    
	for (my $ncRNA_line=0; $ncRNA_line <= $#{$p_ncRNAs->{$nc_chr}}; $ncRNA_line++) {
	    my $ncRNA = ${$p_ncRNAs->{$nc_chr}}[$ncRNA_line];
	  
	    my ($nc_strand,
		$nc_tr_start,
		$nc_tr_end,
		$nc_exons,
		$nc_exon_starts,
		$nc_exon_ends, 
		$nc_tr_id) = get_parts($ncRNA);
	    
	    my $unique_key = "$nc_tr_id$nc_strand$nc_tr_start$nc_tr_end$nc_exons" .
		@$nc_exon_starts . @$nc_exon_ends;

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
		if (($nc_tr_end - $nc_tr_start) >= $length &&
		    $nc_exons >= $min_exons) { 

		    # Complete overlap reference gene with ncRNA intron.
		    if ($mode =~ m/^conc$/i &&
                        $nc_tr_start < $ref_tr_start &&
                        $nc_tr_end > $ref_tr_end) {

			my $ov_tr_found = $nc_int_tree->fetch($ref_tr_start ,$ref_tr_end);

			my $is_ncRNA_Conc = 0;
			$is_ncRNA_Conc = 1 if (scalar(@$ov_tr_found) >= 1);
						
			#my $is_ncRNA_Conc = is_intronicOverlap($ref_tr_start, $ref_tr_end, $nc_exon_starts, $nc_exon_ends);
			
			if ($is_ncRNA_Conc &&
			    $is_strand_Antisense &&
			    !exists $ncRNA_class->{$unique_key}) {
			    $ncRNA_class->{$unique_key} = "Sense intronic overlap (Conc) with $ref_tr_id;";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line-- if ($ncRNA_line > 0);
                            $ncRNA_line = 0 if ($ncRNA_line = 0);
			    $found++;
			}
			elsif ($is_ncRNA_Conc &&
			       !$is_strand_Antisense &&
			       !exists $ncRNA_class->{$unique_key}) {
			    $ncRNA_class->{$unique_key} = "Sense intronic overlap (Conc) with $ref_tr_id;";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
                            $ncRNA_line-- if ($ncRNA_line > 0);
                            $ncRNA_line = 0 if ($ncRNA_line = 0);
			    $found++;
			}
			elsif (!$is_ncRNA_Conc &&
			    !exists $ncRNA_class->{$unique_key}) {
			    $ncRNA_class->{$unique_key} = "Sense intronic overlap (Conc);";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line-- if ($ncRNA_line > 0);
                            $ncRNA_line = 0 if ($ncRNA_line = 0);
			    $found++;
			}
		    }
		    
		    # Complete overlap of ncRNA within reference intron	
		    if ($mode =~ m/^inc$/i &&
			$nc_tr_start >= $ref_tr_start &&
                        $nc_tr_end <= $ref_tr_end) {
			
			my $retain_overlap = $overlap;
			$overlap = 0;
			my $is_ncRNA_Inc = is_intronicOverlap($nc_tr_start, $nc_tr_end, $ref_exon_starts, $ref_exon_ends);
			
			if ($is_ncRNA_Inc &&
			    $is_strand_Antisense &&
			    !exists $ncRNA_class->{$unique_key}) {
			    $ncRNA_class->{$unique_key} = "Antisense intronic overlap (Inc) with $ref_tr_id;";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line-- if ($ncRNA_line > 0);
			    $ncRNA_line = 0 if ($ncRNA_line = 0);
			    $found++;
			}
			elsif ($is_ncRNA_Inc &&
			       !$is_strand_Antisense &&
			       !exists $ncRNA_class->{$unique_key}) {
			    $ncRNA_class->{$unique_key} = "Sense intronic overlap (Inc) with $ref_tr_id;";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line-- if ($ncRNA_line > 0);
			    $ncRNA_line = 0 if ($ncRNA_line = 0);
			    $found++;
			}
			elsif ($is_ncRNA_Inc &&
			       !exists $ncRNA_class->{$unique_key}) {
			    $ncRNA_class->{$unique_key} = "Intronic overlap (Inc) with $ref_tr_id;";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line-- if ($ncRNA_line > 0);
			    $ncRNA_line = 0 if ($ncRNA_line = 0);
			    $found++;
			}
			
			$overlap = $retain_overlap;
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
			    $ncRNA_class->{$unique_key} = "Antisense partial intronic overlap (Ponc) with $ref_tr_id;";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line-- if ($ncRNA_line > 0);
			    $ncRNA_line = 0 if ($ncRNA_line = 0);
			    $found++;
			}
			elsif (!$is_ncRNA_exonicOverlap &&
			       !$is_strand_Antisense &&
			       !exists $ncRNA_class->{$unique_key} &&
			       scalar(@$found_intron_ov >= 1)) {
			    $ncRNA_class->{$unique_key} = "Sense partial intronic overlap (Ponc) with $ref_tr_id;";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line-- if ($ncRNA_line > 0);
			    $ncRNA_line = 0 if ($ncRNA_line = 0);
			    $found++;
			}
			elsif (!$is_ncRNA_exonicOverlap &&
			       !exists $ncRNA_class->{$unique_key} &&
			       scalar(@$found_intron_ov >= 1)) {
			    $ncRNA_class->{$unique_key} = "Partial intronic overlap (Ponc) with $ref_tr_id;";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line-- if ($ncRNA_line > 0);
			    $ncRNA_line = 0 if ($ncRNA_line = 0);
			    $found++;
			}
			$overlap = $retain_overlap;
		    }
		    
		    # Antisense and Sense reference exonic overlap.
		    if ($mode =~ m/^exonic$/i) {
			
			my $is_ncRNA_exonicOverlap = is_exonicOverlap($ref_exon_starts, $ref_exon_ends, $nc_exon_starts, $nc_exon_ends);

			if ($is_ncRNA_exonicOverlap &&
			    $is_strand_Antisense &&
			    !exists $ncRNA_class->{$unique_key}) {
			    $ncRNA_class->{$unique_key} = "Antisense exonic overlap with $ref_tr_id;";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line-- if ($ncRNA_line > 0);
			    $ncRNA_line = 0 if ($ncRNA_line == 0);
			    $found++;
			}
			elsif ($mode =~ m/^exonic$/i &&
			       $is_ncRNA_exonicOverlap &&
			       !$is_strand_Antisense &&
			       !exists $ncRNA_class->{$unique_key}) {
			    $ncRNA_class->{$unique_key} = "Sense exonic overlap with $ref_tr_id;";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line-- if ($ncRNA_line > 0);
			    $ncRNA_line = 0 if ($ncRNA_line = 0);
			    $found++;
			}
			elsif ($mode =~ m/^exonic$/i &&
			       $is_ncRNA_exonicOverlap &&
			       !exists $ncRNA_class->{$unique_key}) {
			    $ncRNA_class->{$unique_key} = "Exonic overlap with $ref_tr_id;";
			    splice(@{$p_ncRNAs->{$nc_chr}}, $ncRNA_line, 1);
			    $ncRNA_line-- if ($ncRNA_line > 0);
			    $ncRNA_line = 0 if ($ncRNA_line = 0);
			    $found++;
			}
		    }
		}
	    }
	    $io->execute_system_command("grep -iP \'$nc_tr_id\' $p_gtf | sed -e \'s\/\$\/ $ncRNA_class->{$unique_key}\/\' >> $c_ncRNAs", 0)
		if ($ncRNA_class->{$unique_key} ne '');
	}
    }
    return $found;
}

# Store genome coordinates.
sub store_coords {
    my $f = shift;
    my $fh = $io->open_file('<', $f);
    my $store = {};

    $io->c_time('Reading information from gene prediction format file [ ' . 
		$file_basename->($f, 'suffix') . ' ] ...', $quiet);
    while (my $line = <$fh>) {
	chomp $line;
	$line = $io->strip_leading_and_trailing_spaces($line);
	my ($t_id, $chr, $strand, $tr_start, $tr_end,
	    $cds_start, $cds_end, $num_exons, $exon_starts, $exon_ends, @rem) = split/\t/, $line;
	$exon_starts =~ s/\,$//;
	$exon_ends =~ s/\,$//;
	$io->error('Supplied file [ ' . $file_basename->($f, 'suffix') . ' ] does not seem to be in gene prediction format...' .
		   "\nError occured on line:\n\n$line\n")
	    if ($chr !~ m/^chr/i || $strand !~ m/^\+|\-|\.$/ || $num_exons !~ m/\d+/);
	
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
	
	for (0 .. $#$s_ex_st) {
	    if ($c_ex_st->[$nc_ex_coord] <= $s_ex_st->[$_] && 
		$c_ex_st->[$nc_ex_coord] <= $s_ex_end->[$_] && 
		$c_ex_end->[$nc_ex_coord] >= $s_ex_st->[$_] && 
		$c_ex_end->[$nc_ex_coord] <= $s_ex_end->[$_]) {
		return 1 if ( defined($overlap) &&
			      ( ( ( ($c_ex_end->[$_] - $s_ex_st->[$_]) / ($s_ex_end->[$_] - $s_ex_st->[$_]) ) * 100 ) >= $overlap) );
	    }
	    elsif ($c_ex_st->[$nc_ex_coord] >= $s_ex_st->[$_] &&
		   $c_ex_st->[$nc_ex_coord] <= $s_ex_end->[$_] &&
		   $c_ex_end->[$nc_ex_coord] >= $s_ex_st->[$_] &&
		   $c_ex_end->[$nc_ex_coord] <= $s_ex_end->[$_]) {
		return 1 if ( defined($overlap) &&
			      ( ( ( ($c_ex_end->[$_] - $c_ex_st->[$_]) / ($s_ex_end->[$_] - $s_ex_st->[$_]) ) * 100 ) >= $overlap) );
	    }
	    elsif ($c_ex_st->[$nc_ex_coord] >= $s_ex_st->[$_] &&
		   $c_ex_st->[$nc_ex_coord] <= $s_ex_end->[$_] &&
		   $c_ex_end->[$nc_ex_coord] >= $s_ex_st->[$_] &&
		   $c_ex_end->[$nc_ex_coord] >= $s_ex_end->[$_]) {
		return 1 if ( defined($overlap) &&
			      ( ( ( ($s_ex_end->[$_] - $c_ex_st->[$_]) / ($s_ex_end->[$_] - $s_ex_st->[$_]) ) * 100 ) >= $overlap) );
	    }
	    elsif ($c_ex_st->[$nc_ex_coord] <= $s_ex_st->[$_] &&
		   $c_ex_st->[$nc_ex_coord] <= $s_ex_end->[$_] &&
		   $c_ex_end->[$nc_ex_coord] >= $s_ex_st->[$_] && 
		   $c_ex_end->[$nc_ex_coord] >= $s_ex_end->[$_]) {
		return 1 if ( defined($overlap) &&
			      ( ( ( ($s_ex_end->[$_] - $s_ex_st->[$_]) / ($s_ex_end->[$_] - $s_ex_st->[$_]) ) * 100 ) >= $overlap) );
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
    return 0 if ($#$ex_st < 1);

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
