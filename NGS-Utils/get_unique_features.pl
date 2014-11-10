#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use IO::Routine;
use Set::IntervalTree;

my ($LASTCHANGEDBY) = q$LastChangedBy: konganti $ =~ m/.+?\:(.+)/;
my ($LASTCHANGEDDATE) = q$LastChangedDate: 2014-11-03 11:44:27 -0500 (Mon, 03 Nov 2014)  $ =~ m/.+?\:(.+)/;
my ($VERSION) = q$LastChangedRevision: 0515 $ =~ m/.+?(\d+)/;
my $AUTHORFULLNAME = 'Kranti Konganti';

my $io = IO::Routine->new();
my $s_time = $io->start_timer();

my ($help, $quiet, $sc1, $sc2, $cc1, $cc2,
    $sf, $cf, $chr_s, $chr_c,
    %seen_coord, %store_s_coords, %seen_s,
    %seen_cat_tr, $extract4pipeline, %logFoldCuff,
    $unique, $known, $pipe2stdout, $j_fh, $overlap, 
    $tr_coords, $tr_start_col,
    $tr_end_col, $trans_bd_on_chr, $keyword, $keyword_col,
    $sff, $cff, $logfold4cuff, $new_f, $output);

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
				 'known'                      => \$known,
				 'overlap=f'                  => \$overlap,
				 'no-exon-match=s'            => \$tr_coords,
				 'stdout'                     => \$pipe2stdout,
				 'compare-keyword|ck=s'       => \$keyword,
				 'keyword-col|kc=i'           => \$keyword_col,
				 'source-file-format|sff=s'   => \$sff,
				 'compare-file-format|cff=s'  => \$cff,
				 'extract'                    => \$extract4pipeline,
				 'fpkm-logfold'               => \$logfold4cuff,
				 'output=s'                   => \$output
				);

# Check if sff or cff is defined
($chr_s, $sc1, $sc2) = gtf_or_bed($sff);
($chr_c, $cc1, $cc2) = gtf_or_bed($cff);

$io->verify_options([$is_valid_option, $sc1, $sc2, $cc1, $cc2,
                     $chr_s, $chr_c],
		    $help);

$io->verify_files([$sf, $cf],
                  ['Chromosome feature ( Source )',
                   'Chromosome feature ( Compare )']);

$io->this_script_info($io->file_basename($0),
		      $VERSION,
		      $AUTHORFULLNAME,
		      $LASTCHANGEDBY,
		      $LASTCHANGEDDATE, '',
		      $quiet);

# In Perl index starts at 0.

$sc1--;
$sc2--;
$cc1--;
$cc2--;
$chr_s--;
$chr_c--;
$keyword_col-- if(defined($keyword_col) && ($keyword_col ne ''));

my ($cf_filename, $path, $suffix) = $io->file_basename($cf, 'all');

$io->warning('Both --known and --unique switches are on. Getting only known ncRNAs...')
    if (defined $unique && defined $known);
undef $unique if ( defined $known || (defined $known && defined $unique) );

if (defined $unique) {
  $io->execute_system_command(0,
			      'Getting unique features...',
			      $quiet);
  $suffix = '.unique.' . lc($cff);
}
else {
  $io->execute_system_command(0,
			      'Getting common features...',
			      $quiet);
  $suffix = '.common.' . lc($cff);
}
 
if (defined $pipe2stdout ) {
    $j_fh = *STDOUT;
}
elsif (!defined $pipe2stdout && defined $output) {
    $new_f = $output;
    $j_fh = $io->open_file('>', $new_f);
    $io->execute_system_command(0,
                                "New file will be $new_f",
                                $quiet);
}
elsif (!defined $pipe2stdout && !defined $extract4pipeline) {
    $new_f = $path . $io->file_basename($cf) . $suffix;
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
    $io->error('Transcript coordinates must be numeric and be separated by a comma');
}
elsif ($is_valid_option eq '') {
    $io->error("Invalid options. See $0 -h for help." );
}

if (( defined $logfold4cuff &&
      $sff !~ m/^gtf$/i &&
      $cff !~ m/^gtf$/) ||
    ( defined $unique && 
      defined $logfold4cuff )) {
    $io->error('Source and compare files must be in GTF format for log fold change of FPKM / RPKM values to be calculated between 2 files.' . 
	       "\nAlso, this option (-fpkm-logfold) can only be used when -unique option is not mentioned.");
}

# Store source coordinates in memory.
my $tr_feat_num = 0;
while (my $line = <$s_fh>) {
    
    chomp $line;
    my ($left_coords, $right_coords) = [];
    
    $line = $io->strip_leading_and_trailing_spaces($line);
    my @cols = split/\t/, $line;
    $io->error('Cannot find chromosome column in file [ ' . $sf . ' ]' .
	       "\nError occured on line:\n\n" . $line) if ($cols[$chr_s] !~ m/^chr/i);

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
	    
	    if ($logfold4cuff && !defined($unique)) {
		$logFoldCuff{$cols[$chr_s] . $left_coords->[$_] . $right_coords->[$_]} = $cols[$#cols];
	    }
	}
    }

    if (defined ($tr_coords)) {
	($tr_start_col, $tr_end_col) = split/\,/, $tr_coords;
	$tr_start_col--;
	$tr_end_col--;
	push @{$trans_bd_on_chr->{$cols[$chr_s]}}, "$cols[$tr_start_col]|$cols[$tr_end_col]";
    }
    $tr_feat_num++;
}

$io->c_time('Known non-coding RNAs found [ ' . $io->file_basename($sf, 'suffix') . ' ]: ' . $tr_feat_num)
    if (!$quiet && $tr_feat_num && !defined $extract4pipeline);

# Check if file has correct transcript-exon features.
my $is_correct_tr_exon_file = $io->execute_get_sys_cmd_output("grep -iP '\ttranscript\t' $cf | wc -l");
$io->warning('Correct transcript-exon features not present in the GTF file [ ' . $io->file_basename($cf, 'suffix') . 
	     ' ]. This may cause erratic behaviour and incorrect predictions.')
    if ($cff =~ m/gtf/i && $is_correct_tr_exon_file =~ m/STDERR/i);

# Now, extract either common or unique features. First print header if user wants log fold change values
print $j_fh "# GeneID\tTranscriptID\tLocus\tFPKM_S\tGeneID\tTranscriptID\tLocus\tFPKM_C\tlog2(FPKM_C/FPKM_S)\n"
    if (defined $logfold4cuff);

my $check_line = 0;
my $insert_line = 0;
while (my $line = <$c_fh>) {
    chomp $line;
    $line = $io->strip_leading_and_trailing_spaces($line);
    
    my @cols = split/\t/, $line;
    $io->error('Cannot find chromosome column in file [ ' . $cf . ' ]' .
	 "\nError occured on line:\n\n" . $line) if ($cols[$chr_c] !~ m/^chr/i);
    $cols[$chr_c] = lc($cols[$chr_c]);
    #print "$cols[$chr_c]\t$cols[$cc1]\t$cols[$cc2]\n";

    if (defined($keyword) && ($keyword ne '') &&
	defined($keyword_col) && ($keyword_col ne '')) {
	next if ($cols[$keyword_col] ne $keyword);
    }
       
    my ($from_cat_ncRNA_tr_id) = ($line =~ m/(transcript\_id\s+\".+?\")/);
    
    $overlap = 0 if (!defined($overlap) || $overlap eq '');
    #print $line, "\n";

    if (exists $store_s_coords{$cols[$chr_c]}) {
	my $uq_tr_tree = Set::IntervalTree->new() if (defined($unique));
	
        foreach my $left_coord (sort {$a <=> $b} keys %{$store_s_coords{$cols[$chr_c]}}) {
            foreach my $right_coord (sort {$a <=> $b} values @{$store_s_coords{$cols[$chr_c]}{$left_coord}}) {

		my $s_ex_len = $right_coord -  $left_coord;
		my $c_ex_len = $cols[$cc2] - $cols[$cc1];
		my $ex_ov_per = ($s_ex_len / $c_ex_len) * 100 if ($c_ex_len > 0);

		# Set new tree for "unique"
		if (defined $unique) {
		    $uq_tr_tree->insert($cols[$chr_c] . ":$left_coord-$right_coord | $insert_line", $left_coord, $right_coord);
		    $insert_line++;
		}


		if (defined $logfold4cuff && !defined($unique)) {
		    my ($s_fpkm_line, $c_fpkm_line) = get_cuff_values($logFoldCuff{$cols[$chr_c] . $left_coord . $right_coord}, $cols[$#cols]);
		    $line = "$s_fpkm_line->[0]\t$s_fpkm_line->[1]\t$cols[$chr_c]:$left_coord-$right_coord\t" .
		    "$s_fpkm_line->[2]\t$c_fpkm_line->[0]\t$c_fpkm_line->[1]\t$cols[$chr_c]:$cols[$cc1]-$cols[$cc2]\t$c_fpkm_line->[2]\t" .
		    sprintf('%.5f', log($c_fpkm_line->[2] / $s_fpkm_line->[2]) / log(2));
		}
		

                if ($cols[$cc1] <= $left_coord && $cols[$cc1] <= $right_coord && $cols[$cc2] >= $left_coord && $cols[$cc2] <= $right_coord) {
		    if (!defined($unique) && defined($overlap) && ($ex_ov_per >= $overlap) && !is_duplicate($line)) {
			if (defined $extract4pipeline) {
			    extract_tr($from_cat_ncRNA_tr_id);
			}
			else {
			    print $j_fh $line, "\n";
			}
		    }
		}
                elsif ($cols[$cc1] >= $left_coord && $cols[$cc1] <= $right_coord && $cols[$cc2] >= $left_coord && $cols[$cc2] <= $right_coord) {
                    if (!defined($unique) && defined($overlap) && ($ex_ov_per >= $overlap) && !is_duplicate($line)) {
                        if (defined $extract4pipeline) {
			    extract_tr($from_cat_ncRNA_tr_id);
			}
			else {
			    print $j_fh $line, "\n";
			}
                    }
                }
                elsif ($cols[$cc1] >= $left_coord && $cols[$cc1] <= $right_coord && $cols[$cc2] >= $left_coord && $cols[$cc2] >= $right_coord) {
                    if (!defined($unique) && defined($overlap) && ($ex_ov_per >= $overlap) && !is_duplicate($line)) {
                        if (defined $extract4pipeline) {
			    extract_tr($from_cat_ncRNA_tr_id);
			}
			else {
			    print $j_fh $line, "\n";
			}
                    }
                }
		elsif ($cols[$cc1] <= $left_coord && $cols[$cc1] <= $right_coord && $cols[$cc2] >= $left_coord && $cols[$cc2] >= $right_coord) {
                    if (!defined($unique) && defined($overlap) && ($ex_ov_per >= $overlap) && !is_duplicate($line)) {
                        if (defined $extract4pipeline) {
			    extract_tr($from_cat_ncRNA_tr_id);
			}
			else {
			    print $j_fh $line, "\n";
			}
                    }
		}
	    }
        }

	if (defined($unique)) {
	    next if ($line !~ m/\ttranscript\t/i);
	    my $is_unique_feat = $uq_tr_tree->fetch($cols[$cc1], $cols[$cc2]);
	    
	    #if (scalar(@$is_unique_feat) > 0) {
	    #	my $this_chr = (split /\:/, $is_unique_feat->[0])[0];
	    #	print @$is_unique_feat, "\t", "$cols[$chr_c]:$cols[$cc1]-$cols[$cc2]", "\n"
	    #	    if($this_chr eq $cols[$chr_c]);
	    #	exit 0;
	    #}
	    
	    next if (scalar(@$is_unique_feat));
	    $check_line++;
	    
	    if (defined $extract4pipeline) {
		extract_tr($from_cat_ncRNA_tr_id);
	    }
	    else {
		my $tr_lines = extract_tr($from_cat_ncRNA_tr_id, 'GET');
		print $j_fh $$tr_lines;
	    }
	}     
	
	if (defined($unique) && defined($tr_coords)) {

	    my $tr_tree = Set::IntervalTree->new();
	    my $seen_tr_tree = {};

	    foreach my $tr_coord (values @{$trans_bd_on_chr->{$cols[$chr_c]}}) {

		my ($tr_start, $tr_end) = split /\|/, $tr_coord;
		my $id = $cols[$chr_c] . ':' . $tr_start . '-' . $tr_end;

		if (!exists $seen_tr_tree->{$id}) {
		    $tr_tree->insert($id, $tr_start, $tr_end);
		    $seen_tr_tree->{$id} = 1;
		}
	    }
	    
	    my $tr_intersect = $tr_tree->fetch($cols[$cc1], $cols[$cc2]);
	    	    
	    if (scalar(@$tr_intersect) > 0) {
		if (defined $extract4pipeline) {
		    extract_tr($from_cat_ncRNA_tr_id);
		}
		else {
		    print $j_fh $line, "\ti\n";
		}
	    }
	    elsif (scalar(@$tr_intersect) == 0) {
		if (defined $extract4pipeline) {
		    extract_tr($from_cat_ncRNA_tr_id);
		}
		else {
		    print $j_fh $line, "\ti\n";
		}
	    }
	}
    }
}

$io->c_time('Unique transcript features found [ ' . $io->file_basename($new_f, 'suffix') . ' ]: ' . $check_line)
    if (!$quiet && $check_line && !defined $extract4pipeline);

$io->end_timer($s_time, $quiet);

close $s_fh;
close $c_fh;
close $j_fh if (defined $j_fh);


#################################### Functions ##########################################

# Get unique / common features from categorize_nRNAs module
sub extract_tr {
    my $tr_id = shift;
    my $get = shift;
    $tr_id =~ s/"/\\\"/g;
    $tr_id =~ s/\./\\\./g;
    #print "grep -P \"$tr_id\" $cf";
    
    if ($get && $get eq 'GET') {
	my $tr_lines = $io->execute_get_sys_cmd_output("grep -P \"$tr_id\" $cf");
	return \$tr_lines;
    }
    else {
	$io->execute_system_command("grep -P \"$tr_id\" $cf");
    }
    return;
}

sub get_cuff_values {
    my $line1 = shift;
    my $line2 = shift;
    
    my ($s_gene_id, $s_tr_id, $s_fpkm) = ($line1 =~ m/gene_id\s+\"(.+?)\".*?transcript_id\s+\"(.+?)\".+?[F|R]PKM\s+\"(.+?)\"/i);
    my ($c_gene_id, $c_tr_id, $c_fpkm) = ($line2 =~ m/gene_id\s+\"(.+?)\".*?transcript_id\s+\"(.+?)\".+?[F|R]PKM\s+\"(.+?)\"/i);
    
    return([$s_gene_id, $s_tr_id, $s_fpkm], [$c_gene_id, $c_tr_id, $c_fpkm]);
    #return "$s_gene_id\t$s_tr_id\t$cols[$chr_c]:$left_coord-$right_coord\t" .
	#"$s_fpkm\t$c_gene_id\t$c_tr_id\t$cols[$chr_c]:$cols[$cc1]-$cols[$cc2]\t$c_fpkm\t" .
	#sprintf('%.5f', log($c_fpkm / $s_fpkm) / log(2)) . "\n";
}

# Avoid duplicates

sub is_duplicate {
  my $line = shift;
  if (!exists $seen_coord{$line}) {
    $seen_coord{$line} = 1;
    return 0;
  }
  return 1;
}

# Return GTF or BED or GFF columns

sub gtf_or_bed {
    my $format = shift;
    if (defined($format) && $format ne '' && $format =~ m/^bed$/i) {
	return ("1", "2", "3");
    }
    elsif (defined($format) && $format ne '' && $format =~ m/^gtf|gff$/i) {
	return ("1", "4", "5");
    }
    elsif (defined($format) && $format ne '' && $format =~ m/^gtf|gff|bed$/i) {
	$io->error('Invalid file format [ ' . $format . ' ] specified!' .
		   'Currently, only GTF, GFF or BED file format is supported.');
    }
    return;
}

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

    perl get_unique_features.pl -scc 1 -ccc 1 -s1 2 -s2 3 -c1 2 -c2 3 -sf refGene.bed -cf novel_ncRNA.bed


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

=item -u or --unique (Optional)

  Get unique features between source and comparison file.
  Default: disabled

=item -sff or --source-file-format (Optional)

  You can directly specify one of the allowed file formats (gtf, gff or bed)
  for the source file and can skip -scc, -sc1 and -sc2.

=item -cff or --compare-file-format (Optional)

  You can directly specify one of the allowed file formats (gtf, gff or bed)
  for the compare file and can skip -ccc, -cc1 and -cc2.

=item -scc or --source-chr-column (Required)

  Column number of the Source file's chromosome information.
  This optional can be skipped if -sff is mentioned.
  
=item -s1 or --source-column-1 (Required)

  Column number of the Source file's chromosome left coordinate (start coordinate) information.
  This option can be skipped if -sff is mentioned.
 
=item -s2 or --source-column-2 (Required)

  Column number of the Source file's chromosome right coordinate (end coordinate) information.
  This option can be skipped if -sff is mentioned.

=item -ccc or --compare-chr-column (Required)

  Column number of the Compare file's chromosome information.
  This option can be skipped if -cff is mentioned.

=item -c1 or --compare-column-1 (Required)

  Column number of the Source file's chromosome left coordinate (start coordinate) information.
  This option can be skipped if -cff is mentioned.

=item -c2 or --compare-column-2 (Required)

  Column number of the Source file's chromosome right coordinate (end coordinate) information.
  This option can be skipped if -cff is mentioned.

=item -sf or --source-file (Required)

  Path to "Source" gene feature file.

=item -cf or --compare-file (Required)

  Path to "Compare" gene feature file.

=item -stdout (Optional)

  Print output to STDOUT instead of a file.

=item -ov or --overlap (Optional)
    
  Extract common features that overlap with source feature by at least this much percentage.

=item -out or --output (Optional)
    
  Print output to this file. If not mentioned, new output file will be created in current
  working directory.

=item -no or --no-exon-match (Optional)
   
 Extract features that fall either within entire reference intron or in intergenic regions.

=item -ck or --compare-keyword (Optional)

 The script stores the feature coordinates of the "Source file" in memory so that the feature
 coordinates of "Compare file" can be compared. If you want to compare features whose column 
 contains only this keyword, then use this option.
 
 Ex: if you want to just compare the transcript coordinates and not the exon coordinates 
     of a GTF file, you can use this option as -ck 'transcript'.

=item -kc or --keyword-col (Optional)

 The column number of the --compare-keyword.

 Ex: For GTF file, it is -ck 'transcript' -kc '3'

=item -extract

 Extract unique or common features using grep from the compare file.

=item -fpkm-logfold

 This is a special case option and works only on assembled transcripts' file in GTF format from
 Cufflinks since it contains proper gene_id, transcript_id and FPKM description values in the last
 column. Any GTF file that has this information can be used as input with this option. Using this 
 option will extract common features between 2 input GTF files and outputs a file similar to output
 from cuffdiff program excluding the test_statistic, p_value and q_value. This option was designed
 to be used with lncRNApipe pipeline. 

=back

=head1 AUTHOR

Kranti Konganti, E<lt>konganti@tamu.eduE<gt>.

=head1 COPYRIGHT

This program is distributed under the Artistic License.

=head1 DATE

Nov-03-2014

=cut
