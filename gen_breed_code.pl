#!/usr/bin/perl

use strict;
use warnings;
use Carp;
use IO::Routine;
use Getopt::Long;

my $io = IO::Routine->new();
$io->check_and_load_modules(['Pod::Usage', 'Getopt::Long', 'File::Basename',
                             'Bio::SeqIO']);

my ($help, $quiet, $bo_dir, $output, $bo_code, $pos);

my $is_valid_option = GetOptions('help|?' => \$help,
                                 'quiet' => \$quiet,
                                 'bo-dir=s' => \$bo_dir,
                                 'output=s' => \$output,
                                 'pos=i' => \$pos #Skip nth nucleotide
                                 );


$io->verify_options($is_valid_option);
$output = $io->validate_create_path($output, 'create', 'Output directory');
$bo_dir = $io->validate_create_path($bo_dir, 'do not create', 'Breed of Origin directory');

opendir(BO_DIR, $bo_dir) || $io->error("Cannot open $bo_dir for reading files.\n");
my @files = readdir BO_DIR;
close BO_DIR;

$pos--;

$bo_code->{'NA'} = 1;
$bo_code->{'AN'} = 1;
$bo_code->{'AA'} = 0;
$bo_code->{'NN'} = 2;
$bo_code->{'A-'} = 4;
$bo_code->{'N-'} = 3;
$bo_code->{'-A'} = 4;
$bo_code->{'-N'} = 3;
$bo_code->{'--'} = '-';

foreach my $file (@files) {
    #print $file, "\n";
    next if ($file eq '.' || $file eq '..' || $file !~ m/\.txt$/);
    my ($write_filename, $write_path, $write_suff) = File::Basename::fileparse($file, qr/\.[^.]*/);
    $write_filename = $output . $write_filename . '.code.txt';
    my $w_fh = $io->open_file($write_filename, '>');

    $file = $bo_dir . $file;
    my $seq_in_obj = Bio::SeqIO->new(-file => $file,
                                     -format => 'fasta');
    while (my $seq_in = $seq_in_obj->next_seq) {
        my $seq_in_2 = $seq_in_obj->next_seq;
        gen_code($seq_in->id, $seq_in->seq, $seq_in_2->id, $seq_in_2->seq, $w_fh);
    }

}
########### Functions #################
sub gen_code {
    my $h1_id = shift;
    my $h1_seq = shift;
    my $h2_id = shift;
    my $h2_seq = shift;
    my $fh = shift;

    my @nucs_1 = split//, $h1_seq;
    my @nucs_2 = split//, $h2_seq;

    print $fh ">$h1_id\n";

    for (0 .. $#nucs_1) {
        next if ($_ == $pos);
        #print "$nucs_1[$_]$nucs_2[$_]\n";
        print $fh $bo_code->{"$nucs_1[$_]$nucs_2[$_]"};

    }
    print $fh "\n\n";
    return;
}