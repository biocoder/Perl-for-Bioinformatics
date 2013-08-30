#!/opt/perl/bin/perl

use strict;
use warnings;
use Getopt::Long;
use IO::Routine;
use File::Basename;

my $io = IO::Routine->new();

my ($help, $quiet, $sc1, $sc2, $cc1, $cc2,
    $sf, $cf, $chr_s, $chr_c,
    %seen_coord, %store_s_coords, %seen_s,
    $unique);

my $is_valid_option = GetOptions('source-column-1|sc-1|s1=s'  => \$sc1,
				 'source-column-2|sc-2|s2=s'  => \$sc2,
				 'compare-column-1|cc-1|c1=s' => \$cc1,
				 'compare-column-2|cc-2|c2=s' => \$cc2,
				 'source-file|sf|s=s'         => \$sf,
				 'compare-file|cf|c=s'        => \$cf,
				 'source-chr-column|scc=s'    => \$chr_s,
				 'compare-chr-column|ccc=s'   => \$chr_c,
				 'help'                       => \$help,
				 'quiet'                      => \$quiet,
				 'unique'                     => \$unique
				);

$io->verify_options([$is_valid_option, $sc1, $sc2, $cc1, $cc2,
                     $chr_s, $chr_c],
		    $help);

$io->verify_files([$sf, $cf],
                  ['Chromosome feature ( Source )',
                   'Chromosome feature ( Compare )']);

# In Perl index starts at 0.

$sc1--;
$sc2--;
$cc1--;
$cc2--;
$chr_s--;
$chr_c--;

my ($filename, $path, $suffix) = fileparse($cf, qr/\.[^.]*/);
if (defined($unique)) {
  $io->execute_system_command(0,
			      'Getting unique features...',
			      $quiet);
  $suffix = '.unique.txt';
}
else {
  $io->execute_system_command(0,
			      'Getting common features...',
			      $quiet);
  $suffix = '.common.txt';
}

my $new_f = $path . basename($sf) . '_' . basename($cf) . $suffix;

$io->execute_system_command(0,
                            "New file will be $new_f",
                            $quiet);

my $s_fh = $io->open_file('<', $sf);
my $c_fh = $io->open_file('<', $cf);
my $j_fh = $io->open_file('>', $new_f);

# Store source coordinates in memory.

while (my $line = <$s_fh>) {
  chomp $line;
  my @cols = split/\t/, $line;
  $io->error('Cannot find chromosome column') if ($cols[$chr_s] !~ m/^chr/i);
  $cols[$chr_s] = lc ($cols[$chr_s]);

  if (!exists $seen_s{$cols[$chr_s] . $cols[$sc1] . $cols[$sc2]}) {
      $seen_s{$cols[$chr_s] . $cols[$sc1] . $cols[$sc2]} = 1;
      push @{$store_s_coords{$cols[$chr_s]}{$cols[$sc1]}}, $cols[$sc2];
  }
}

# Now, extract either common or unique features.

while (my $line = <$c_fh>) {
    chomp $line;
    my @cols = split/\t/, $line;
    $io->error('Cannot find chromosome column') if ($cols[$chr_c] !~ m/^chr/i);
    $cols[$chr_c] = lc($cols[$chr_c]);
    #print "$cols[$chr_c]\t$cols[$cc1]\t$cols[$cc2]\n";

    if (exists $store_s_coords{$cols[$chr_c]}) {
        foreach my $left_coord (sort {$a <=> $b} keys %{$store_s_coords{$cols[$chr_c]}}) {
            foreach my $right_coord (sort {$a <=> $b} values @{$store_s_coords{$cols[$chr_c]}{$left_coord}}) {
                if ($cols[$cc1] <= $left_coord && $cols[$cc1] <= $right_coord && $cols[$cc2] >= $left_coord && $cols[$cc2] <= $right_coord) {
                    print $j_fh $line, "\n" if (!is_duplicate($line) && !defined $unique);
                }
                elsif ($cols[$cc1] >= $left_coord && $cols[$cc1] <= $right_coord && $cols[$cc2] >= $left_coord && $cols[$cc2] <= $right_coord) {
                    print $j_fh $line, "\n" if (!is_duplicate($line) && !defined $unique);
                }
                elsif ($cols[$cc1] >= $left_coord && $cols[$cc1] <= $right_coord && $cols[$cc2] >= $left_coord && $cols[$cc2] >= $right_coord ) {
                    print $j_fh $line, "\n" if (!is_duplicate($line) && !defined $unique);
                }
            }
        }
        foreach my $left_coord (sort {$a <=> $b} keys %{$store_s_coords{$cols[$chr_c]}}) {
            foreach my $right_coord(sort {$a <=> $b} values @{$store_s_coords{$cols[$chr_c]}{$left_coord}}) {
                if (!is_duplicate($line) && defined $unique) {
                    print $j_fh $line, "\n";
                }
            }
        }
    }
}

close $s_fh;
close $c_fh;
close $j_fh;

# Avoid duplicates

sub is_duplicate {
  my $line = shift;
  if (!exists $seen_coord{$line}) {
    $seen_coord{$line} = 1;
    return 0;
  }
  return 1;
}

__END__

=head1 NAME

get_unique_features.pl

=head1 SYNOPSIS

This script will extract unique or common features between 2 chromosome feature files of different formats like
GTF, GFF, BED etc...

Examples:

    perl get_unique_features.pl -h

    perl get_unique_features.pl get_unique_features.pl -scc 1 -ccc 1 -s1 2 -s2 3 -c1 2 -c2 3 -sf refGene.bed -cf novel_ncRNA.bed


=head1 DESCRIPTION

This script makes it easy to find out the common overlapping or non-overlapping features of the genome given two gene feature files in different file formats. For example, the BED format file stores the gene boundary information in first, second and third columns respectively, whereas the refGene.txt file stores the reference gene boundary information in second, fourth and fifth columns respectively. Now, if you want to find out whether the features in the BED file overlap or do not overlap with reference genes in refGene.txt file, you specify the refGene.txt file as "Source file" and BED file as "Compare file". Then the Source file's column information would be 2, 4, 5 (i.e -scc 2, -s1 4, -s2 5) respectively and the Compare file's column information will be 1, 2 and 3 (i.e -ccc 1, -c1 2, -c2 3) respectively. The output file will be created at the location of the "Compare" file. This script is strand agnostic, meaning that the features either common or unique are reported irrespective of the Source feature's strand information.


=head1 OPTIONS

get_unique_features.pl takes the following arguments:

=over 4

=item -h or --help (Optional)

  Displays this helpful message.

=item -q or --quiet (Optional)

  Providing this option suppresses the log messages to the shell.
  Default: disabled

=item -scc or --source-chr-column (Required)

  Column number of the Source file's chromosome information.

=item -s1 or --source-column-1 (Required)

 Column number of the Source file's chromosome left coordinate (start coordinate) information.

=item -s2 or --source-column-2 (Required)

Column number of the Source file's chromosome right coordinate (end coordinate) information.

=item -ccc or --compare-chr-column (Required)

  Column number of the Compare file's chromosome information.

=item -c1 or --compare-column-1 (Required)

  Column number of the Source file's chromosome left coordinate (start coordinate) information.

=item -c2 or --compare-column-2 (Required)

  Column number of the Source file's chromosome right coordinate (end coordinate) information.

=item -sf or --source-file (Required)

  Path to "Source" gene feature file.

=item -cf or --compare-file (Required)

  Path to "Compare" gene feature file.

=back

=head1 AUTHOR

Kranti Konganti, E<lt>konganti@tamu.eduE<gt>.

=head1 COPYRIGHT

This program is distributed under the Artistic License.

=head1 DATE

Aug-30-2013

=cut
