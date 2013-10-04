#!/opt/perl/bin/perl

use strict;
use warnings;
use Getopt::Long;
use IO::Routine;
use File::Basename;
use Set::IntervalTree;

my ($LASTCHANGEDBY) = q$LastChangedBy: konganti $ =~ m/.+?\:(.+)/;
my ($LASTCHANGEDDATE) = q$LastChangedDate: 2013-10-04 09:08:44 -0500 (Fri, 04 Oct 2013) $ =~ m/.+?\:(.+)/;
my ($VERSION) = q$LastChangedRevision: 62 $ =~ m/.+?(\d+)/;
my $AUTHORFULLNAME = 'Kranti Konganti';

my $io = IO::Routine->new();
my $s_time = $io->start_timer();

my ($help, $quiet, $sc1, $sc2, $cc1, $cc2,
    $sf, $cf, $chr_s, $chr_c,
    %seen_coord, %store_s_coords, %seen_s,
    $unique, $pipe2stdout, $j_fh, $overlap, 
    $tr_coords, $tr_start_col,
    $tr_end_col, $trans_bd_on_chr);

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
				 'unique'                     => \$unique,
				 'overlap=i'                  => \$overlap,
				 'no-exon-match=s'            => \$tr_coords,
				 'stdout'                     => \$pipe2stdout
				);

$io->verify_options([$is_valid_option, $sc1, $sc2, $cc1, $cc2,
                     $chr_s, $chr_c],
		    $help);

$io->verify_files([$sf, $cf],
                  ['Chromosome feature ( Source )',
                   'Chromosome feature ( Compare )']);

$io->this_script_info($0, $VERSION, $AUTHORFULLNAME, $LASTCHANGEDBY, $LASTCHANGEDDATE, $quiet);

# In Perl index starts at 0.

$sc1--;
$sc2--;
$cc1--;
$cc2--;
$chr_s--;
$chr_c--;

my ($cf_filename, $path, $suffix) = $io->file_basename($cf, 'all');

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

if (defined $pipe2stdout ) {
    $j_fh = *STDOUT;
}
else {
    my $new_f = $path . $io->file_basename($sf) . '_' . $cf_filename . $suffix;
    $j_fh = $io->open_file('>', $new_f);
    $io->execute_system_command(0,
				"New file will be $new_f",
				$quiet);
}

my $s_fh = $io->open_file('<', $sf);
my $c_fh = $io->open_file('<', $cf);

# If asked to find intronic and intergenic features, start a tree
if (defined($tr_coords) && 
    ($tr_coords ne '') &&
     ($tr_coords =~ m/^\d+\s*\,\s*\d+$/)) {
    $tr_coords =~ s/\s+//g;
}
elsif (defined($tr_coords) &&
       ($tr_coords ne '') &&
       ($tr_coords !~ m/\d+|\,/)) {
    $io->error('Transcript coordinates must be numeric and be separated by a comma')
}

# Store source coordinates in memory.

while (my $line = <$s_fh>) {
    
    chomp $line;
    my ($left_coords, $right_coords) = [];
    
    $line = $io->strip_leading_and_trailing_spaces($line);
    my @cols = split/\t/, $line;
    $io->error('Cannot find chromosome column in file [ ' . $sf . ' ]') if ($cols[$chr_s] !~ m/^chr/i);

    $cols[$chr_s] = lc ($cols[$chr_s]);

    if ($cols[$sc1] =~ m/\,/ && $cols[$sc2] =~ m/\,/) {
	chop $cols[$sc1] if ($cols[$sc1] =~ m/\,$/);
	chop $cols[$sc2] if ($cols[$sc2] =~ m/\,$/);
	@$left_coords = split/\,/, $cols[$sc1];
	@$right_coords = split/\,/, $cols[$sc2];
	$io->error('Number of left coordinates is not equal to number of right coordinates')
	    if (scalar(@$left_coords) != scalar(@$right_coords));
    }
    else {
	push @{$left_coords}, $cols[$sc1];
	push @{$right_coords}, $cols[$sc2];
    }
    
    for (0 .. $#$left_coords) {
	if (!exists $seen_s{$cols[$chr_s] . $left_coords->[$_] . $right_coords->[$_]}) {
	    $seen_s{$cols[$chr_s] . $left_coords->[$_] . $right_coords->[$_]} = 1;
	    push @{$store_s_coords{$cols[$chr_s]}{$left_coords->[$_]}}, $right_coords->[$_];
	}
    }

    if (defined ($tr_coords)) {
	($tr_start_col, $tr_end_col) = split/\,/, $tr_coords;
	$tr_start_col--;
	$tr_end_col--;
	push @{$trans_bd_on_chr->{$cols[$chr_s]}}, "$cols[$tr_start_col]|$cols[$tr_end_col]";
    }
}
# Now, extract either common or unique features.

while (my $line = <$c_fh>) {
    chomp $line;
    $line = $io->strip_leading_and_trailing_spaces($line);
    my @cols = split/\t/, $line;
    $io->error('Cannot find chromosome column in file [ ' . $cf . ' ]') if ($cols[$chr_c] !~ m/^chr/i);
    $cols[$chr_c] = lc($cols[$chr_c]);
    #print "$cols[$chr_c]\t$cols[$cc1]\t$cols[$cc2]\n";

    if (exists $store_s_coords{$cols[$chr_c]}) {
        foreach my $left_coord (sort {$a <=> $b} keys %{$store_s_coords{$cols[$chr_c]}}) {
            foreach my $right_coord (sort {$a <=> $b} values @{$store_s_coords{$cols[$chr_c]}{$left_coord}}) {
                if ($cols[$cc1] <= $left_coord && $cols[$cc1] <= $right_coord && $cols[$cc2] >= $left_coord && $cols[$cc2] <= $right_coord) {
		    last if (defined($overlap) && (($cols[$cc2] - $left_coord) <= $overlap));
                    print $j_fh $line, "\n" if (!is_duplicate($line) && !defined $unique);
                }
                elsif ($cols[$cc1] >= $left_coord && $cols[$cc1] <= $right_coord && $cols[$cc2] >= $left_coord && $cols[$cc2] <= $right_coord) {
		    last if (defined($overlap) && (($cols[$cc2] - $cols[$cc1]) <= $overlap));
                    print $j_fh $line, "\n" if (!is_duplicate($line) && !defined $unique);
                }
                elsif ($cols[$cc1] >= $left_coord && $cols[$cc1] <= $right_coord && $cols[$cc2] >= $left_coord && $cols[$cc2] >= $right_coord) {
		    last if (defined($overlap) && (($right_coord - $cols[$cc1]) <= $overlap));
                    print $j_fh $line, "\n" if (!is_duplicate($line) && !defined $unique);
                }
		elsif ($cols[$cc1] <= $left_coord && $cols[$cc1] <= $right_coord && $cols[$cc2] >= $left_coord && $cols[$cc2] >= $right_coord) {
		    last if (defined($overlap) && (($right_coord - $left_coord) <= $overlap));
		    print $j_fh $line, "\n" if (!is_duplicate($line) && !defined $unique);
		}
            }
        }


	if (!is_duplicate($line) && defined($unique) && defined($tr_coords)) {

	    my $tr_tree = Set::IntervalTree->new();
	    my $seen_tr_tree = {};

	    foreach my $tr_coord (values @{$trans_bd_on_chr->{$cols[$chr_c]}}) {

		my ($tr_start, $tr_end) = split /\|/, $tr_coord;

		if (!exists $seen_tr_tree->{"$cols[$chr_c]:$tr_start-$tr_end"}) {
		    $tr_tree->insert("$cols[$chr_c]:$tr_start-$tr_end",
				     $tr_start, $tr_end);
		    $seen_tr_tree->{"$cols[$chr_c]:$tr_start-$tr_end"} = 1;
		}
	    }
	    
	    my $tr_intersect = $tr_tree->fetch($cols[$cc1], $cols[$cc2]);
	    
	    if (scalar(@$tr_intersect) > 0) {
		print $j_fh $line, "\ti\n";
	    }
	    elsif (scalar(@$tr_intersect) == 0) {
		print $j_fh $line, "\tu\n";
	    }
	}

	if (!is_duplicate($line) && defined($unique)) {
	    print $j_fh $line, "\n";
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

$io->end_timer($s_time, $quiet);

__END__

=head1 NAME

get_unique_features.pl

=head1 SYNOPSIS

This script will extract unique or common features between 2 chromosome feature files of different formats like
GTF, GFF, BED etc...

Complete Description: 

    perldoc get_unique_features.pl

Examples:

    perl get_unique_features.pl -h

    perl get_unique_features.pl get_unique_features.pl -scc 1 -ccc 1 -s1 2 -s2 3 -c1 2 -c2 3 -sf refGene.bed -cf novel_ncRNA.bed


=head1 DESCRIPTION

This script makes it easy to find out the common overlapping or non-overlapping features of the genome given two gene feature files in different file formats. For example, the BED format file stores the gene boundary information in first, second and third columns respectively, whereas the refGene.txt file stores the reference gene boundary information in second, fourth and fifth columns respectively. Now, if you want to find out whether the features in the BED file overlap or do not overlap with reference genes in refGene.txt file, you specify the refGene.txt file as "Source file" and BED file as "Compare file". Then the Source file's column information would be 2, 4, 5 (i.e -scc 2, -s1 4, -s2 5) respectively and the Compare file's column information will be 1, 2 and 3 (i.e -ccc 1, -c1 2, -c2 3) respectively. You can directly use refGene.txt as source file and use exonStarts and exonEnds columns to extract common or unique features. The output file will be created at the location of the "Compare" file. The output can be printed to STDOUT with -stdout option. This script is strand agnostic, meaning that the features either common or unique are reported irrespective of the Source feature's strand information.

=head1 EXAMPLE

get_unique_features.pl -sf rn4_refGene.txt -cf H3K9me3_mapped-W200-G200-E100.bed -ccc 1 -c1 2 -c2 3 -scc 3 -s1 10 -s2 11 -u --no-exon-match '5,6'

The above command will extract H3K9me3 peaks that fall within entire reference intron or in intergenic regions when comparing the BED file [H3K9me3_mapped-W200-G200-E100.bed] with reference annotation. In this case, the reference annotation [refGene.txt] becomes the source file (-sf) and [H3K9me3_mapped-W200-G200-E100.bed] becomes the compare file (-cf), where -scc 1 -s1 10 -s2 11 --no-exon-match '5,6' are the column numbers of Chromosome (1), ExonStarts (10), ExonEnds (11), TranscriptStart (5) TranscriptEnd (6) features for refGene.txt file and -ccc 1 -c1 2 -c2 3 are the column numbers of Chromosome (1), PeakBoundaries (2, 3) for H3K9me3_mapped-W200-G200-E100.bed file.


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

=item -stdout (Optional)

  Print output to STDOUT instead of a file.

=item -ov or --overlap (Optional)
    
  Extract features that are unique or overlap by this many number of bases

=item --no-exon-match (Optional)
   
 Extract features that fall either within entire reference intron or in intergenic regions

=back

=head1 AUTHOR

Kranti Konganti, E<lt>konganti@tamu.eduE<gt>.

=head1 COPYRIGHT

This program is distributed under the Artistic License.

=head1 DATE

Aug-30-2013

=cut
