#!/opt/perl/bin/perl

# Change the shebang line above to point to your Perl installation

use strict;
use warnings;
use Carp;
use Getopt::Long;
use IO::Routine;

my $io = IO::Routine->new();

my ($LASTCHANGEDBY) = q$LastChangedBy: biocoder $ =~ m/.+?\:(.+)/;
my ($LASTCHANGEDDATE) = q$LastChangedDate: 2013-08-29 12:25:31 -0600 (Thu, 29 Aug 2013) $ =~ m/.+?\:(.+)/;
my ($VERSION) = q$LastChangedRevision: 19 $ =~ m/.+?(\d+)/;
my $AUTHORFULLNAME = 'Kranti Konganti';

my ($quiet, $refGeneTxt, $bedOut, $help);
my $is_valid_option = GetOptions ('help|?'     => \$help,
                                  'quiet'      => \$quiet,
                                  'bed=s'      => \$bedOut,
				  'refGene=s'  => \$refGeneTxt
				 );

$io->verify_options($is_valid_option);
$io->verify_files([$refGeneTxt], ['refGene.txt']);
my $refGeneTxt_fh = $io->open_file('<', $refGeneTxt);
my $bed_fh = $io->open_file('>', $bedOut);

$io->this_script_info($0, $VERSION, $AUTHORFULLNAME, $LASTCHANGEDBY, $LASTCHANGEDDATE, $quiet);


$io->execute_system_command(0,
			   'Validating refGene.txt File ...',
			   $quiet);

while (my $line = <$refGeneTxt_fh>) {
  $line = $io->strip_leading_and_trailing_spaces($line);
  my @refGeneCols = split/\t/, $line;

  $io->error("It seems like the $refGeneTxt file has invalid number of columns")
    if (scalar(@refGeneCols) != 16);

  print $bed_fh "$refGeneCols[2]\t$refGeneCols[4]\t$refGeneCols[5]\t$refGeneCols[12]\t$refGeneCols[11]\t$refGeneCols[3]\n";
}

close $bed_fh;
close $refGeneTxt_fh;

__END__

=head1 NAME

refGene_to_bed.pl

=head1 SYNOPSIS

This script will convert refGene.txt file to BED format.

Examples:

    perl refGene_to_bed.pl -h

    perl refGene_to_bed.pl -r refGene.txt -b bedOutput.txt

    perl refGene_to_bed.pl -r refGene.txt -b bedOutput.txt -q

=head1 DESCRIPTION

This script will use the following columns of the refGene table from UCSC to generate BED format.

  chrom
  strand
  txStart
  txEnd
  score
  name2

=head1 OPTIONS

refGene_to_bed.pl takes the following arguments:

=over 4

=item -h or --help (Optional)

  Displays this helpful message.

=item -q or --quiet (Optional)

  Providing this option suppresses the log messages to the shell.
  Default: disabled

=item -r or --refGene (Required)

  Path to the refGene.txt file you want to conver to BED.

=item -b or --bed (Required)

  Path to file to save the output.

=back

=head1 AUTHOR

Kranti Konganti, E<lt>konganti@tamu.eduE<gt>.

=head1 COPYRIGHT

This program is distributed under the Artistic License.

=head1 DATE

Aug-29-2013

=cut
