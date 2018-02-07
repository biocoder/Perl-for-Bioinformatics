#!/usr/bin/env perl

=head1 NAME

split_fasta_seqs.pl

=head1 SYNOPSIS

split directory containing multifasta files into files with specified
number of sequences of specified chunk size

Examples:

    perl split_fasta_seqs.pl --help

    perl split_fasta_seqs.pl --fasta_dir_or_fasta_file /path/to/dir --destination_dir /path/to/result/chunk/files --seqlen 500 --no-of-seqs 1000

=head1 DESCRIPTION

This script will take the input directory containing multifasta files and split them into small files with specified sequence size (chunk size),
limiting each file to specified number of sequences

=head1 ARGUMENTS

split_fasta_seqs.pl takes the following arguments:

=over 4

=item -h or --help  (Optional)

  Displays the usage message.

=item -f or --fasta_dir_or_fasta_file (Required)

  Path to directory containing fasta files to be processed or
  path to fasta file.

=item -d or --destination_dir (Optional)

  Path to directory where the split files will be stored. If not mentioned,
  new directory will be created at the same path as the fasta file or directory
  containing fasta files.

=item -s or --seqlen (Optional)

  How long should each sequence should be?
  If it is not mentioned, total sequence length will be used.

=item -n or --no-of-seqs (Required)

  How many sequences should each file contain?

=item -r or --remove-dir (Optional)

  @ ****************************************************************** @
  |    !!! KNOW WHAT YOU ARE DOING WHEN MENTIONING THIS OPTION  !!!    |
  |        ===================================================         |
  |      Remove source fasta directory after script has finished.      |
  @ ****************************************************************** @

=back

=head1 AUTHOR

Kranti Konganti, E<lt>konganti@tamu.eduE<gt>.

=head1 COPYRIGHT

This program is distributed under the Artistic License 2.0.

=head1 DATE

Feb-14-2012

=cut

my ($CHANGEDBY) = q$LastChangedBy: konganti $ =~ m/.+?\:(.+)/;
my ($LASTCHANGEDDATE) = q$LastChangedDate: 2012-12-01 10:58:31 -0600 (Sat, 01 Dec 2012) $ =~ m/.+?\:(.+)/;
my ($VERSION) = q$LastChangedRevision: 43 $ =~ m/.+?(\d+)/;
my $AUTHORFULLNAME = 'Kranti Konganti';

use strict;
use warnings;
use Carp;

my $use_module_list = ['Bio::SeqIO', 'Getopt::Long', 'Pod::Usage', 
		       'File::Basename'];

# Check if modules are installed else exit with generous message rather than perl's die message
check_if_modules_installed($use_module_list);
this_script_info();

# Printing new line to make error messages separated from shell prompt
print "\n";

my ($help, $fasta_dir_or_fasta_file, $destination_dir, $seqlen, $no_of_seqs,
    $remove_dir);

my $is_valid_option = GetOptions('help|?' => \$help,
                                 'fasta_dir_or_fasta_file|fasta_file=s' => \$fasta_dir_or_fasta_file,
                                 'destination_dir=s' => \$destination_dir,
                                 'seqlen=i' => \$seqlen,
                                 'no-of-seqs=i' => \$no_of_seqs,
                                 'remove-dir' => \$remove_dir);

if (!$is_valid_option) {
    print "\nSee $0 -h for usage\n\n";
    exit;
}

if ($help) {
    pod2usage(-exitval => 1,
              -msg => "\n");
}

if (!$fasta_dir_or_fasta_file) {
    usage ("--fasta_dir_or_fasta_file");
}

if (!$seqlen) {
    $seqlen = 'use_total';
}

if (!$no_of_seqs) {
    usage ("--no_of_seqs");
}

if (-f $fasta_dir_or_fasta_file) {
    my ($file, $path, $suffix) = fileparse($fasta_dir_or_fasta_file, qr/\.[^.]*/);

    # if path is current working directory return pwd result since fileparse does not want to
    # do that.
    if ($path eq './') {
        $path = `pwd`;
        chomp $path;
        $path .= '/';
        $fasta_dir_or_fasta_file = $path . $fasta_dir_or_fasta_file;
    }

    my $tmp_fasta_source_dir = $path . 'tmp/fasta_file_to_be_splitted';
    my $tmp_link_to_orig_fasta_file = $tmp_fasta_source_dir . '/' . $file . $suffix;

    print "\nCreating temporary directory $tmp_fasta_source_dir\n",
    `mkdir -p $tmp_fasta_source_dir`
        if (!-d $tmp_fasta_source_dir);

    if (!-e $tmp_link_to_orig_fasta_file) {
        print "\nCreating link to original fasta file\n";
        `ln -s $fasta_dir_or_fasta_file $tmp_fasta_source_dir/`;
    }

    $fasta_dir_or_fasta_file = $tmp_fasta_source_dir;
}

if (!$destination_dir) {
    my ($file, $path, $suffix) = fileparse($fasta_dir_or_fasta_file, qr/\.[^.]*/);
    my $new_fasta_destination_dir = $path . 'splitted_fasta_files';
    $destination_dir = $new_fasta_destination_dir;

    print "\n--destination_dir not defined\n\nSplitted files will be stored to $destination_dir\n\n";
}

# Add trailing / to destination and source directories if they are not present
if ($fasta_dir_or_fasta_file !~ m/\/$/) {
    $fasta_dir_or_fasta_file .= '/';
}

if ($destination_dir !~ m/\/$/) {
    $destination_dir .= '/';
}

if (!-d $destination_dir) {
    `mkdir -p $destination_dir`;
}

if (is_dir_empty($fasta_dir_or_fasta_file)) {
    die "\n$fasta_dir_or_fasta_file does not contain any files.\n\n";
}

opendir (FASTA_DIR, $fasta_dir_or_fasta_file) || die "\nCannot open $fasta_dir_or_fasta_file: $!\n\n";
my @fasta_files = readdir FASTA_DIR;

foreach my $fasta_file (@fasta_files) {
    if ($fasta_file eq '.' || $fasta_file eq '..' || -d $fasta_file) {
        next;
    }

    my $path_to_fasta_file = $fasta_dir_or_fasta_file . $fasta_file;
    my ($fasta_file_name, $fasta_file_path, $suffix) = fileparse($path_to_fasta_file, qr/\.[^.]*/);

    print "Splitting $path_to_fasta_file ...\n\n";

    # Read fasta
    my $seq_file_in = Bio::SeqIO->new(-file => $path_to_fasta_file,
                                      -format => 'fasta');

    # If the fasta sequence size of a single sequence in multifasta file is
    # less than the chunk size, then just store it upto the no_of_seqs and
    # write the multifasta file

    my $seq_file_no = 0;
    my @seqs_upto_mentioned_no_of_seqs;

    while (my $seq_in = $seq_file_in->next_seq) {
        my $seq_length = $seq_in->length;

        if ($seqlen eq 'use_total' || $seq_length < $seqlen) {
            push (@seqs_upto_mentioned_no_of_seqs, $seq_in);

            if (((scalar(@seqs_upto_mentioned_no_of_seqs) % $no_of_seqs) == 0) && (scalar(@seqs_upto_mentioned_no_of_seqs) > 0)) {
                $seq_file_no++;
                my $new_fasta_file = $destination_dir . $fasta_file_name . '_' . $seq_file_no . $suffix;
                write_multifasta($fasta_file_name, $new_fasta_file, \@seqs_upto_mentioned_no_of_seqs);
                @seqs_upto_mentioned_no_of_seqs = ();
            }

        }
        elsif ($seq_length > $seqlen) {
            my $sub_seq_id_no = 1;

            for (my $i=1; $i<$seq_length; $i=$i+$seqlen) {
                my $sub_seq_id = $seq_in->id . '{split_seq_' . $sub_seq_id_no . '}';
                my $sub_seq_end = (($i + $seqlen) > $seq_length) ? ($seq_length + 1) : ($i + $seqlen);
                my $sub_seq;

                if ($sub_seq_end == $seq_length + 1) {
                    my $j = $i - $seqlen;
                    $sub_seq = $seq_in->trunc($j, $sub_seq_end - 1);
                }
                else {
                    $sub_seq = $seq_in->trunc($i, $sub_seq_end - 1);
                }

                $sub_seq->id($sub_seq_id);
                $sub_seq->desc($seq_in->desc);
                push (@seqs_upto_mentioned_no_of_seqs, $sub_seq);

                if (((scalar(@seqs_upto_mentioned_no_of_seqs) % $no_of_seqs) == 0) && (scalar(@seqs_upto_mentioned_no_of_seqs) > 0)) {
                    $seq_file_no++;
                    my $new_fasta_file = $destination_dir . $fasta_file_name . '_' . $seq_file_no . $suffix;
                    write_multifasta($fasta_file_name, $new_fasta_file, \@seqs_upto_mentioned_no_of_seqs);
                    @seqs_upto_mentioned_no_of_seqs = ();
                }

                $sub_seq_id_no++;
            }
        }
    }

    # Flush out remaining seqs
    $seq_file_no++;
    my $last_new_fasta_file_name = $destination_dir . $fasta_file_name . '_' . $seq_file_no . $suffix;
    write_multifasta($fasta_file_name, "$last_new_fasta_file_name", \@seqs_upto_mentioned_no_of_seqs);

    $seq_file_in->close;
    close FASTA_DIR;
}

if ($remove_dir) {
    print "Removing $fasta_dir_or_fasta_file ...\n";
    `rm -rf $fasta_dir_or_fasta_file`
}
print "Done!\n\n";

sub write_multifasta {
    my ($fasta_file_name, $new_fasta_file, $seqs) = @_;

    my $query_fasta = Bio::SeqIO->new (-file => ">$new_fasta_file",
                                       -format => "fasta");
    $query_fasta->write_seq(@$seqs);
    $query_fasta->close;
}

sub check_if_modules_installed {
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

sub usage {
    my $msg = shift;
    print "\n$msg not defined. This option is required\n";
    print "\nSee $0 -h for usage\n\n";
    exit;
}

sub is_dir_empty {
    my $dirname = shift;

    opendir (my $dh, $dirname) || die "\n$dirname not a directory: $!\n\n";
    my $there_are_files = 0;

    while ( defined (my $file = readdir $dh) ) {
        my $path_to_file = $dirname . $file;

        if ($file eq '.' || $file eq '..') {
            next;
        }
        elsif (-f $path_to_file) { # if atleast 1 file is present
            $there_are_files++;
            last;
        }
    }

    close $dh;

    if ($there_are_files) {
        return 0;
    }
    else {
        return 1;
    }
}

# Shell msg that differentiates log from error
sub error {
    my $msg = shift;
    print "\nERROR!\n------\n$msg\n\n";
    pod2usage(-exitval => 2);
}

# Shell msg for warning                                                                                                                                              
sub warning {
    my $msg = shift;
    print "\nWARNING!\n--------\n$msg\n\n";
    return;
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
