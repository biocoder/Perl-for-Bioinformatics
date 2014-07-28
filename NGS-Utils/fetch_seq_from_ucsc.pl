#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long;
use LWP::Simple;
use XML::XPath;
use XML::XPath::XMLParser;
use IO::Routine;

my ($LASTCHANGEDBY) = q$LastChangedBy: konganti $ =~ m/.+?\:(.+)/;
my ($LASTCHANGEDDATE) = q$LastChangedDate: 2013-10-09 12:46:11 -0500 (Wed, 09 Oct 2013) $ =~ m/.+?\:(.+)/;
my ($VERSION) = q$LastChangedRevision: 64 $ =~ m/.+?(\d+)/;
my $AUTHORFULLNAME = 'Kranti Konganti';

# Declare initial global variables
my ($quiet, $tmap, $output, $help, $dbkey, $print_seq_fh,
    $chr_coords, $chr_info, $id_re, $skip_re, $file_format,
    $seq_desc, $seq_descs, $match_re);

my $is_valid_option = GetOptions ('help|?'     => \$help,
                                  'quiet'      => \$quiet,
                                  'output=s'   => \$output,
                                  'tmap=s'     => \$tmap,
                                  'dbkey=s'    => \$dbkey,
				  'chr-cols=s' => \$chr_coords,
				  'id-re=s'    => \$id_re,
				  'skip-re=s'  => \$skip_re,
				  'match-re=s' => \$match_re,
				  'ff=s'       => \$file_format,
				  'seq-desc=s' => \$seq_desc
                                  );

# Print info if not quiet
my $io = IO::Routine->new($help, $quiet);
my $s_time = $io->start_timer;

$io->this_script_info($io->file_basename($0),
                      $VERSION,
                      $AUTHORFULLNAME,
                      $LASTCHANGEDBY,
                      $LASTCHANGEDDATE);
		      


# Check for the validity of options
$io->verify_options([$is_valid_option, $dbkey,
		     $tmap, $output]);

$io->c_time('Analysis started...');

$io->c_time('Verifying file [ ' .
	    $io->file_basename($tmap) .
	    ' ] ...');

$io->verify_files([$tmap],
		  ['TMAP']);

$output = $io->validate_create_path($output, 'create', 'Output');

$io->c_time("FASTA files will be stored at $output ...");
$io->c_time('Fetching Sequences ...');

if (defined $seq_desc) {
    $seq_descs = [split/\,/, $seq_desc]; 
    for (0 .. $#$seq_descs) {
	$seq_descs->[$_]--;
    }
}

$chr_coords = join(',', gtf_or_bed($file_format));

if (defined $chr_coords && $chr_coords ne '') {
    $chr_coords = $io->strip_leading_and_trailing_spaces($chr_coords);
    $chr_coords =~ s/\s+//g;
    $chr_info = [split /\,/, $chr_coords];
    for (0 .. $#$chr_info) {
	$chr_info->[$_]--;
    }
}

my $tmap_fh = $io->open_file('<', $tmap);
my ($tmap_filename, $tmap_filepath, $suffix) = $io->file_basename($tmap, 'all');
my $fa = $output . $tmap_filename . '.fa';

if (!-e $fa && !-s $fa) {
    $print_seq_fh = $io->open_file('>', $fa)
}
else {
    $print_seq_fh = $io->open_file('>>', $fa);
}

$io->c_time('Checking for required GNU core utils...', $quiet);
$io->check_sys_level_cmds(['grep'],
                          ['2.6.3']);

# For some clean output
print STDOUT "\n";

my $seqs_fetched =  0;
while (my $line = <$tmap_fh>) {
    my ($u_seq_id, $unique_seq_id) = '';

    next if ($line =~ m/^ref\_gene\_id/ ||
	     $line =~ m/^$/);

    next if (defined $skip_re &&
	     $skip_re ne '' &&
	     $line =~ qr/$skip_re/);

    next if (defined $match_re &&
	     $match_re ne '' &&
	     $line !~ qr/$match_re/);

    ($u_seq_id) = ($line =~ qr/$id_re/) if (defined $id_re && $id_re ne '');
      
    chomp $line;
    my @cols = split/\t/, $line;
    
    if (ref($chr_info) eq 'ARRAY') {
	$cols[13] = lc($cols[$chr_info->[0]]);
	$cols[14] = $cols[$chr_info->[1]];
	$cols[15] = $cols[$chr_info->[2]];
	if ($u_seq_id) {
	    $u_seq_id = $io->strip_leading_and_trailing_spaces($u_seq_id);
	    $unique_seq_id = $u_seq_id . '|' . $cols[13] . ':' . $cols[14] . '-' . $cols[15]; 
	}
	else {
	    $unique_seq_id = $cols[13] . ':' . $cols[14] . '-' . $cols[15]; 
	}
    }
    else {
	$unique_seq_id = $cols[4] . '|' . "class_code:$cols[2]" . '_' . lc($cols[13]) . ':' . $cols[14] . '-' . $cols[15];
    }

    my $u_seq_desc = '';
    if (defined $seq_desc) {
	for (0 .. $#$seq_descs) {
	    $u_seq_desc .= "$cols[$seq_descs->[$_]]|";
	}
	chop $u_seq_desc if ($u_seq_desc =~ m/\|$/);
    }
    $unique_seq_id .= " $u_seq_desc";

    if ($cols[13] !~ m/^chr/i &&
	$cols[14] !~ m/^\d+$/ &&
	$cols[15] !~ m/^\d+$/) {
	$io->error("Cannot find chromosome information in columns...\n".
		   "Encountered line is:\n" .
		   $line);
    }

    $cols[13] =~ s/chr//;

    my $fetched_seq = $io->execute_get_sys_cmd_output("grep -A 1 \'$unique_seq_id\' $fa");
    chomp $fetched_seq;

    if ($fetched_seq !~ m/STDERR/i) {
	print STDOUT "Skipping $unique_seq_id ... Already fetched!\n" if (!$quiet);
        next;
    }
    else {
	print STDOUT "Querying $unique_seq_id against UCSC DAS ...\n" if (!$quiet);
    }

    my $xml = get("http://genome.ucsc.edu/cgi-bin/das/$dbkey/dna?segment=$cols[13]:$cols[14],$cols[15]");
    my $xpath = XML::XPath->new(xml=>$xml);
    my $xml_nodes = $xpath->find('/DASDNA/SEQUENCE/DNA/text()');
    print_seq(\$xml_nodes, $unique_seq_id);
    $seqs_fetched++;
}

$io->c_time("$seqs_fetched sequences fetched...");

$io->end_timer($s_time);
close $tmap_fh;
close $print_seq_fh;

###################### Functions ###############################################

# Print fasta sequences;
sub print_seq {
    my $xml_nodes = shift;
    my $seq_id = shift;

    foreach my $xml_node ($$xml_nodes->get_nodelist) {
        my $seq = $xml_node->getValue;
	chomp $seq;
        $seq = $io->strip_leading_and_trailing_spaces($seq);
        print $print_seq_fh ">$seq_id\n$seq\n";
    }
    return;
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
    return $chr_coords;
}

__END__

=head1 NAME

fetch_seq_from_ucsc.pl

=head1 SYNOPSIS

This script will fetch DNA seqeuences from UCSC database given a TMAP output
file that contains chromosome information.

Examples:

    perl fetch_seq_from_ucsc.pl -h

    perl fetch_seq_from_ucsc.pl --tmap proximal_vs_reference.cmp.transcripts.gtf.class_codes_u_and_i.tmap -out /home/konganti/fasta -dbkey rn4

=head1 DESCRIPTION

This script uses the TMAP file, which is generally an ouput format from Cufflinks
pipeline. This TMAP file must contain chromosome information in the last three
columns, mainly chromosome number, chromosome start and chromosome end. It will
use that information to query UCSC DAS to fetch DNA sequences.

Example line from a modified TMAP file with chromosome information:

ref_gene_id	ref_id	class_code	cuff_gene_id	cuff_id	FMI	FPKM	FPKM_conf_lo	FPKM_conf_hi	cov	len	major_iso_id	ref_match_len	chr	chr_start	chr_end

Lrp11	NM_001106217	i	CUFF.20	CUFF.20.1	100	1.052506	0.742139	1.362874	2.091452	1080	CUFF.20.1	3114	chr1	2257123	2258202

=head1 OPTIONS

fetch_seq_from_ucsc.pl takes the following arguments:

=over 4

=item -h or --help (Optional)

  Displays this helpful message.

=item -q or --quiet (Optional)

  Providing this option suppresses the log messages to the shell.
  Default: disabled

=item -t or --tmap (Required)

  Modified TMAP file from Cufflinks pipeline that contains chromosome information.

=item -d or --dbkey (Required)

  Database to query (ex: rn4, mm9, hg19 etc...).

=item -o or --output (Required)

  Path to output directory .i.e where should the fasta files be stored?

=item -c or --chr-cols (Optional)

  If you have chromosome information in custom columns, specify it with this option. 
  For example, if you have chromosome id, chromosome start, chromosome end in 
  columns 1, 2 and 3 respectively, specify as: -c '1,2,3'.

=item -id or --id-re (Optional)
    
  If you want to extract an ID based on pattern matched from the provided input file, 
  provide the regex with this option. For example, -id 'transcript_id.+?\"(.+?)\"' will
  extract the pattern within () and use it as FASTA sequence id.

=item -skip or --skip-re (Optional)

  If you want to skip lines containing this pattern, provide regex with this option.
  For example, -skip 'exon' will skip lines containing the word 'exon'.
  Regex is case sensitive.

=item -match or --match-re (Optional)

  If you want to match lines containing this pattern, provide regex with this option.
  For example, -match 'chr23' will parse lines containing the word 'chr23'.
  Regex is case sensitive.

=item -ff or --ff (Optional)

  You can directly specify one of the allowed file formats (gtf, gff or bed)
  to extract chromosome coordinate information.

=item --seq-desc

  You can describe which column values should be concatenated together to form 
  fasta sequence description id.

=back

=head1 AUTHOR

Kranti Konganti, E<lt>konganti@tamu.eduE<gt>.

=head1 COPYRIGHT

This program is distributed under the Artistic License.

=head1 DATE

Nov-27-2012

=cut
