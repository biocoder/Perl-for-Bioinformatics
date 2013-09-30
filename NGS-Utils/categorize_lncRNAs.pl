#!/usr/bin/env perl

use strict;
use Carp;
use Getopt::Long qw(:config pass_through);
use IO::Routine;
use File::Basename;

my $io = IO::Routine->new();
my $s_time = $io->start_timer();

my ($help, $quiet, $cuffcmp, $genePred, $out, $bed, $sample_names,
    $fpkm_cutoff, $cov_cutoff, $overwrite);
my ($p_file_names_gtf, $p_file_names_txt) = [];

my $is_valid_option = GetOptions('help|?'         => \$help,
				 'quiet'          => \$quiet,
				 'cuffcmp=s'      => \$cuffcmp,
				 'bed=s'          => \$bed,
				 'out=s'          => \$out,
				 'sample-names=s' => \$sample_names,
				 'fpkm-cutoff=f'  => \$fpkm_cutoff,
				 'cov_cutoff=f'   => \$cov_cutoff, 
				 'overwrite'      => \$overwrite);

$io->c_time('Analysis started...', $quiet);

$io->c_time('Verifying options...', $quiet);
$io->verify_options([$is_valid_option, $sample_names],
                    $help);

$fpkm_cutoff = 0.0 if (!defined $fpkm_cutoff || $fpkm_cutoff eq '');
$cov_cutoff = 0.0 if (!defined $fpkm_cutoff || $fpkm_cutoff eq '');

$io->c_time('Validating output path...', $quiet);
my $output = $io->validate_create_path($out, 'create', 
				       'Output directory');

$io->c_time('Checking for required GNU core utils...', $quiet);
$io->check_sys_level_cmds(['grep'], 
			  ['2.6.3']);

$io->c_time('Checking for Cuffcompare tracking and Cufflinks assembled transcript files...');

my $cuffcmp_filename = sub {my @file_p_attrs = fileparse(shift, qr/\.[^.]*/); return $file_p_attrs[0]};
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
    push @{$p_file_names_gtf}, $output . $cuffcmp_filename->($ARGV[$_]) . '.' . $lables[$_] . '.putative_ncRNAs.gtf';
    push @{$p_file_names_txt}, $output . $cuffcmp_filename->($ARGV[$_]) . '.' . $lables[$_] . '.putative_ncRNAs.txt';
    unlink $p_file_names_gtf->[$_] if (-e $p_file_names_gtf->[$_] && defined($overwrite));
}

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
	    my $t_lines = $io->execute_get_sys_cmd_output("grep -iP \'$q_t_id\' $ARGV[$_]",0);
	    if ($t_lines =~ m/.+?FPKM.+?\"(.+?)\".+?cov.+?\"(.+?)\".+/i) {
		$io->execute_system_command("grep -iP \'$q_t_id\' $ARGV[$_] >> $p_file_names_gtf->[$_]")if ($1 >= $fpkm_cutoff && $2 >= $cov_cutoff);
	    }
	}
    }
}

$io->c_time('Converting putative ncRNAs list to Gene Prediction format using gtfToGenePred tool', $quiet);
my $check_for_gtfToGenePred = $io->execute_get_sys_cmd_output('gtfToGenePred', 0);

$io->error('Cannot find gtfToGenePred tool in your path') 
    if ($check_for_gtfToGenePred !~ m/.*?gtfToGenePred.*?convert a GTF file to a genePred/i);

for (0 .. $#$p_file_names_gtf) {
    $io->execute_system_command("gtfToGenePred -genePredExt -geneNameAsName2 $p_file_names_gtf->[0] $p_file_names_txt->[0]", '');
}

print $io->end_timer($s_time, $quiet), "\n";

