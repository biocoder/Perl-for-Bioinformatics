#!/usr/bin/env perl

use strict;
use Carp;
use Getopt::Long qw(:config pass_through);
use IO::Routine;
use File::Basename;

my $io = IO::Routine->new();
my $s_time = $io->start_timer();

my ($help, $quiet, $cuffcmp, $genePred, $out, $sample_names,
    $fpkm_cutoff, $cov_cutoff, $refGenePred, $length, $categorize,
    $min_exons, $overlap);
my ($p_file_names_gtf, $p_file_names_txt) = [];

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
				 'overlap=i'      => \$overlap);


$io->verify_options([$is_valid_option, $sample_names, 
		     $refGenePred, $out],
                    $help);

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


# Get putative list of ncRNAs from Cuffcompare tracking file corresponding
# to cufflinks assembled transcript fragments.
sub get_putative_ncRNAs {
    $io->c_time('Getting putative list of ncRNAs in GTF format...', $quiet);
    while (my $line = <$cuffcmp_fh>) {
	chomp $line;
	$line = $io->strip_leading_and_trailing_spaces($line);
	my ($t_id, $loc_id, $loc_name, $class_code, @cols) = split/\t/, $line;
	for (0 .. $#cols) {
	    if ($class_code =~ m/j|i|o|u|x/i) {
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

# Categorize ncRNAs.
sub class_ncRNAs {
    $io->c_time('Categorizing ncRNAs...', $quiet);
    my $refAnnot = store_coords($refGenePred);
    for (0 .. $#ARGV) {
	my $p_gtf = $p_file_names_gtf->[$_];
	my $p_ncRNAs = store_coords($p_file_names_txt->[$_]);
	my $c_ncRNAs = $output . $file_basename->($ARGV[$_]) . '.' . $lables[$_] . '.class.ncRNAs.txt';
	unlink $c_ncRNAs if (-e $c_ncRNAs);
    

	# Calculate overlap.
	foreach my $nc_chr (keys %{$p_ncRNAs}) {
	    foreach my $ncRNA (values @{$p_ncRNAs->{$nc_chr}}) {
		my ($nc_strand,
		    $nc_tr_start,
		    $nc_tr_end,
		    $nc_exons,
		    $nc_exon_starts,
		    $nc_exon_ends, 
		    $nc_tr_id) = get_parts($ncRNA); 
		
		foreach my $ref_gene (values @{$refAnnot->{$nc_chr}}) {
		    
		    my ($ref_strand, 
			$ref_tr_start, 
			$ref_tr_end, 
			$ref_exons, 
			$ref_exon_starts, 
			$ref_exon_ends,
			$ref_tr_id) = get_parts($ref_gene);
		    
		    # Only lncRNAs, adjustable by user with -length and -min-exons parameters
		    if (($nc_tr_end - $nc_tr_start) >= $length && $nc_exons >= 2) { 
			if ($ref_strand ne $nc_strand &&
			    $ref_strand =~ m/^\+|\-$/ &&
			    $nc_strand =~ m/^\+|\-$/) {  # Antisense ncRNA
			    if (is_exonicOverlap($ref_exon_starts, $ref_exon_ends, $nc_exon_starts, $nc_exon_ends)) {
				$io->execute_system_command("grep -iP \'$nc_tr_id\' $p_gtf | sed -e \'s\/\$\/ Antisense exonic overlap; $ref_tr_id\/\' >> $c_ncRNAs", 0);

			    }
			}
		    }
		}
	    }
	}
    }
    return;
}

# Store genome coordinates.
sub store_coords {
    my $f = shift;
    my $fh = $io->open_file('<', $f);
    my $store = {};

    $io->c_time('Reading information from gene prediction format file [ ' . 
		$file_basename->($f, 'suffix') . ' ]', $quiet);
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

# Convert GTF to gene prediction format.
sub get_genePred {
    $io->c_time('Converting putative ncRNAs list to Gene Prediction format using gtfToGenePred tool', $quiet);
    my $check_for_gtfToGenePred = $io->execute_get_sys_cmd_output('gtfToGenePred', 0);
    
    $io->error('Cannot find gtfToGenePred tool in your path') 
	if ($check_for_gtfToGenePred !~ m/.*?gtfToGenePred.*?convert a GTF file to a genePred/i);
    
    for (0 .. $#$p_file_names_gtf) {
	$io->verify_files([$_], ['GTF']);
	$io->execute_system_command("gtfToGenePred -genePredExt -geneNameAsName2 $p_file_names_gtf->[$_] $p_file_names_txt->[$_]",
	    "gtfToGenePred -genePredExt -geneNameAsName2 $p_file_names_gtf->[$_] $p_file_names_txt->[$_]");
    }
    return;
}

# Split and return columns
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

# Calculate exonic overlap
sub is_exonicOverlap {
    my ($s_ex_st, $s_ex_end, $c_ex_st, $c_ex_end) = @_;
    for (0 .. $#$s_ex_st) {	
	if ($c_ex_st->[$_] <= $s_ex_st->[$_] && $c_ex_st->[$_] <= $s_ex_end->[$_] && $c_ex_end->[$_] >= $s_ex_st->[$_] && $c_ex_end->[$_] <= $s_ex_end->[$_]) {
	    return 1 if ( defined($overlap) && ( ( ( ($c_ex_end->[$_] - $s_ex_st->[$_]) / ($s_ex_end->[$_] - $s_ex_st->[$_]) ) * 100 ) >= $overlap) );
	}
	elsif ($c_ex_st->[$_] >= $s_ex_st->[$_] && $c_ex_st->[$_] <= $s_ex_end->[$_] && $c_ex_end->[$_] >= $s_ex_st->[$_] && $c_ex_end->[$_] <= $s_ex_end->[$_]) {
	    return 1 if ( defined($overlap) && ( ( ( ($c_ex_end->[$_] - $c_ex_st->[$_]) / ($s_ex_end->[$_] - $s_ex_st->[$_]) ) * 100 ) >= $overlap) );
	}
	elsif ($c_ex_st->[$_] >= $s_ex_st->[$_] && $c_ex_st->[$_] <= $s_ex_end->[$_] && $c_ex_end->[$_] >= $s_ex_st->[$_] && $c_ex_end->[$_] >= $s_ex_end->[$_]) {
	   return 1 if ( defined($overlap) && ( ( ( ($s_ex_end->[$_] - $c_ex_st->[$_]) / ($s_ex_end->[$_] - $s_ex_st->[$_]) ) * 100 ) >= $overlap) );
	}
	elsif ($c_ex_st->[$_] <= $s_ex_st->[$_] && $c_ex_st->[$_] <= $s_ex_end->[$_] && $c_ex_end->[$_] >= $s_ex_st->[$_] && $c_ex_end->[$_] >= $s_ex_end->[$_]) {
	    return 1 if ( defined($overlap) && ( ( ( ($s_ex_end->[$_] - $s_ex_st->[$_]) / ($s_ex_end->[$_] - $s_ex_st->[$_]) ) * 100 ) >= $overlap) );
	}
    }
    return 0;
}
