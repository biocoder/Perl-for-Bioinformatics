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
my ($LASTCHANGEDDATE) = q$LastChangedDate: 2015-29-04 17:45:27 -0500 (Tue, 29 Apr 2015)  $ =~ m/.+?\:(.+)/;
my ($VERSION) = q$LastChangedRevision: 0708 $ =~ m/.+?\:\s*(.*)\s*.*/;
my $AUTHORFULLNAME = 'Kranti Konganti';

# Declare initial global variables
my ($quiet, $tmap, $output, $help, $dbkey, $print_seq_fh,
    $chr_coords, $chr_info, $id_re, $skip_re, $file_format,
    $seq_desc, $seq_descs, $match_re, $ncbi_gi, $source_db,
    $pause_ncbi, $local_ref_fa, @local_subseq_fa,
    $ref_fa_in, $ref_fa_out, $contigs, $transcripts,
    $strand, $contig_id);

my $is_valid_option = GetOptions ('help|?'          => \$help,
                                  'quiet'           => \$quiet,
                                  'output=s'        => \$output,
                                  'tmap=s'          => \$tmap,
                                  'dbkey=s'         => \$dbkey,
				  'chr-cols=s'      => \$chr_coords,
				  'id-re=s'         => \$id_re,
				  'skip-re=s'       => \$skip_re,
				  'match-re=s'      => \$match_re,
				  'ff=s'            => \$file_format,
				  'seq-desc=s'      => \$seq_desc,
				  'ncbi-gi-list=s'  => \$ncbi_gi,
				  'pause-ncbi=i'    => \$pause_ncbi,
				  'local-ref-fa=s'  => \$local_ref_fa
                                  );

# Print info if not quiet
my $io = IO::Routine->new($help, $quiet);
my $s_time = $io->start_timer;

$io->this_script_info($io->file_basename($0),
                      $VERSION,
                      $AUTHORFULLNAME,
                      $LASTCHANGEDBY,
                      $LASTCHANGEDDATE, '',
		      $quiet);
		      
# Check for the validity of options
$io->verify_options([$is_valid_option,
		     $tmap, $output]);

$io->error('DBKEY or Local FASTA file or NCBI GI list file not provided.')
    if (!defined $dbkey && !defined $ncbi_gi && !defined $local_ref_fa);

$io->error('Only one of the --dbkey, --local-ref-fa or --ncbi-gi-list options can be defined.')
    if (defined $dbkey && (defined $ncbi_gi || defined $local_ref_fa));

$io->c_time('Analysis started...');

$io->c_time('Verifying file [ ' .
	    $io->file_basename($tmap, 'suffix') .
	    ' ]...');

$io->verify_files([$tmap],
		  ['TMAP']);

$output = $io->validate_create_path($output, 'create', 'Output');

$io->c_time("FASTA files will be stored at $output");
$io->c_time('Fetching Sequences...');

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

if (defined $local_ref_fa) {
    
    require Bio::SeqIO;
    Bio::SeqIO->import();

    $io->error('FASTA reference does not exist.') if (!-e $local_ref_fa || !-s $local_ref_fa);

    $ref_fa_in = Bio::SeqIO->new(-file => $local_ref_fa,
				    -format => 'fasta');
    $contigs = {};
    $io->c_time("Storing reference FASTA in memory...\n");

    while (my $seq_in = $ref_fa_in->next_seq) {
	$contigs->{lc($seq_in->id)} = $seq_in;
	print STDOUT 'FASTA record for ' . lc($seq_in->id) . " stored in memory... \n" if (!$quiet);
    }
    $io->c_time('Now fetching FASTA for transcripts...');
}

# For some clean output
print STDOUT "\n" if (!defined $ncbi_gi);
$pause_ncbi = 1 if (!defined $pause_ncbi && defined $ncbi_gi);
$io->c_time("Program will pause for $pause_ncbi second(s) between queries to NCBI Entrez, " .
	    "so that our IP address is not blocked by NCBI...\n")
    if (defined $ncbi_gi);

while (my $line = <$tmap_fh>) {

    my $u_seq_id = my $unique_seq_id = my $u_seq_desc = '';
    my $chr_id = my $chr_start = my $chr_end = my $chr_strand = '';

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
	$chr_id = lc($cols[$chr_info->[0]]);
	$chr_start = $cols[$chr_info->[1]];
	$chr_end = $cols[$chr_info->[2]];
	$chr_strand = $cols[$chr_info->[3]] if ($cols[$chr_info->[3]]);
	
	if ($u_seq_id) {
	    $u_seq_id = $io->strip_leading_and_trailing_spaces($u_seq_id);
	    $unique_seq_id = $u_seq_id; 
	}
	else {
	    $io->error('Cannot create unique transcript ID with supplied regex ... Bailing out!');
	}
    }
    else {
	$io->error('Cannot create unique transcript ID ... Bailing out!');
    }

    if (defined $seq_desc) {
	for (0 .. $#$seq_descs) {
	    $u_seq_desc .= "$cols[$seq_descs->[$_]]|";
	}
	chop $u_seq_desc if ($u_seq_desc =~ m/\|$/);
    }
    
    $unique_seq_id .= " $u_seq_desc";
    $unique_seq_id = $io->strip_leading_and_trailing_spaces($unique_seq_id);

    if ($chr_start !~ m/^\d+$/ ||
	$chr_end !~ m/^\d+$/) {
	$io->error("Cannot find chromosome information in columns...\n".
		   "\nEncountered line is:\n\n" .
		   $line);
    }

    if (!exists $transcripts->{$unique_seq_id}->{$chr_start}) {
	$transcripts->{$unique_seq_id}->{$chr_start} = $chr_end;
	$contig_id->{$unique_seq_id} = $chr_id;
	$strand->{$unique_seq_id} = $chr_strand;
    }
    else {
	$io->error("Duplicate record found [ for $unique_seq_id ] with same feature start position, when trying to store: " . 
		   "\n\n" . $line . "\n");
    }
}

my $seqs_fetched =  0;
foreach my $unique_seq_id (keys %$transcripts) {
    my $seq = my $exon_coords = '';
    my @exon_str_id;
    
    my $unchanged_contig_id = $contig_id->{$unique_seq_id};
    $contig_id->{$unique_seq_id} =~ s/chr//;

    my $fetched_seq = $io->execute_get_sys_cmd_output("grep -A 1 -P \'\^\>$unique_seq_id\\s+\' $fa");
    chomp $fetched_seq;

    if (defined $ncbi_gi) {
	$source_db = 'NCBI Entrez';
    }
    elsif (defined $local_ref_fa) {
	$source_db = 'Local FASTA file';
    }
    else {
	$source_db = 'UCSC DAS'; 
    }
    
    if ($fetched_seq !~ m/STDERR/i) {
	print STDOUT "Skipping $unique_seq_id ... Already fetched!\n" if (!$quiet);
	next;
    }
    else {
	print STDOUT "Fetching $unique_seq_id from $source_db ...\n" if (!$quiet);
    }

    foreach my $exon_start (sort { $a <=> $b } keys %{$transcripts->{$unique_seq_id}}) {
	my $exon_end = $transcripts->{$unique_seq_id}->{$exon_start};
	
	if ($strand->{$unique_seq_id} eq '-') {
	    push @exon_str_id, "$exon_end-$exon_start";
	}
	else {
	    push @exon_str_id, "$exon_start-$exon_end";
	}
					
	if (defined $ncbi_gi) {
	    my $chr_gi = read_gi($ncbi_gi);
	    
	    $io->warning("GI ID does not exist for [ $unchanged_contig_id ] in the file [ " . $io->file_basename($ncbi_gi, 'suffix') . ' ] ...' .
			 "\nWill not be able to fetch sequence and therefore the transcript will not appear in final list."),
		next if (!exists $chr_gi->{$unchanged_contig_id});
	    
	    my $get_ncbi_seq = get('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?&db=nuccore&id=' .
				  $chr_gi->{$unchanged_contig_id} . "&seq_start=$exon_start&seq_stop=$exon_end&rettype=fasta&retmode=text");
	    
	    $io->warning('NCBI Entrez did not return any sequence from nuccore database.') if (!$get_ncbi_seq || $get_ncbi_seq eq '');
	    
	    $get_ncbi_seq =~ s/^>.+|\n*//g;
	    $seq .= $get_ncbi_seq;
	    sleep $pause_ncbi;
	}
	elsif (defined $local_ref_fa) {
	    if (!exists $contigs->{lc($unchanged_contig_id)}) {
		$io->warning("Contig / Chromosome ID [ $unchanged_contig_id ] does not exist in supplied reference FASTA.\n".
			     'Please make sure that the reference FASTA has those IDs.'); 
		next;
	    }
	    my $subseq_obj = $contigs->{lc($unchanged_contig_id)};
	    my $subseq = $subseq_obj->subseq($exon_start, $exon_end);
	    $seq .= $subseq;
	}
	else {
	    my $xml = get("http://genome.ucsc.edu/cgi-bin/das/$dbkey/dna?segment=" . $contig_id->{$unique_seq_id} . ":$exon_start,$exon_end");
	    $io->warning('UCSC DAS did not return any sequence for database ' . $dbkey .
		       ".\nQueried: http://genome.ucsc.edu/cgi-bin/das/$dbkey/dna?segment=" . $contig_id->{$unique_seq_id} . ":$exon_start,$exon_end")
		if (!$xml || $xml eq '' || $xml !~ m/.*?<\/DNA>/i);
	    my $xpath = XML::XPath->new(xml=>$xml);
	    my $xml_nodes = $xpath->find('/DASDNA/SEQUENCE/DNA/text()');
	    my $sub_seq = return_seq(\$xml_nodes, $unique_seq_id);
	    $sub_seq =~ s/\n*//g;
	    $seq .= $sub_seq;
	}
    }

    if ($strand->{$unique_seq_id} eq '-') {
	my @reverse_exon_coords = reverse @exon_str_id;
	$exon_coords = '[-], ' . join(', ', @reverse_exon_coords);
	$seq = revComp($seq);
    }
    else {
	$exon_coords = '[' . $strand->{$unique_seq_id} . '], ' . join(', ', @exon_str_id);
    }
    
    print $print_seq_fh ">$unique_seq_id $unchanged_contig_id: $exon_coords\n$seq\n";
    $seqs_fetched++;
}

$io->c_time("$seqs_fetched sequences fetched...");

$io->end_timer($s_time);
close $tmap_fh;
close $print_seq_fh;

###################### Functions ###############################################

# Store NCBI GI list
sub read_gi {
    my $gi_list = shift;
    my $gi_list_fh = $io->open_file('<', $gi_list);
    my %chr_gi;
    while (my $line = <$gi_list_fh>) {
	chomp $line;
	$io->error('Supplied GI list file does not seem to be in correct format' .
	    "\nIt should be a tab-separated file with chromosome ids in first column and corresponding NCBI GI ids in second column.\n\n" .
		   "Example:\n--------\n\nchr1\t240254421\nchr3\t332640072\n")
	    if ($line !~ m/^chr.+?\t\d+/i);
	my ($chr, $gi_id) = split/\t+/, $line;
	$chr_gi{lc($chr)} = $gi_id if (!exists $chr_gi{lc($chr)}); 
    }
    return \%chr_gi
}

# Print fasta sequences;
sub return_seq {
    my $xml_nodes = shift;
    my $seq_id = shift;
    my $seq = '';

    foreach my $xml_node ($$xml_nodes->get_nodelist) {
        $seq = $xml_node->getValue;
	chomp $seq;
        $seq = $io->strip_leading_and_trailing_spaces($seq);
    }
    return $seq;
}

# Return GTF or BED or GFF columns
sub gtf_or_bed {
    my $format = shift;
    if (defined($format) && $format ne '' && $format =~ m/^bed$/i) {
        return ("1", "2", "3");
    }
    elsif (defined($format) && $format ne '' && $format =~ m/^gtf|gff$/i) {
        return ("1", "4", "5", "7");
    }
    elsif (defined($format) && $format ne '' && $format =~ m/^gtf|gff|bed$/i) {
        $io->error('Invalid file format [ ' . $format . ' ] specified!' .
                   'Currently, only GTF, GFF or BED file format is supported.');
    }
    return;
}

# Reverse complement DNA
sub revComp {
    my $seq = shift;
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return scalar reverse $seq;
}

__END__

=head1 NAME

fetch_seq_from_ucsc.pl

=head1 SYNOPSIS

This script will fetch DNA seqeuences from UCSC or NCBI database or from
local FASTA genome sequence file, given a GTF file. The ID in the first column 
of the GTF file should match the contig / chromosome ID of the local FASTA
genome sequence file

Examples:

    perl fetch_seq_from_ucsc.pl -h

    perl fetch_seq_from_ucsc.pl --gtf transcripts.rat.putative.class.lncRNAs.unique.gtf -out /home/konganti/fasta -dbkey rn4

=head1 DESCRIPTION

This script uses the GTF file from the lncRNApipe pipeline to create FASTA sequence for 
each transcript with introns spliced out. The tool has the ability to fetch sequences
from UCSC or NCBI databases or from local FASTA genome reference. If "-local" option is
used, the ID in the first column of the GTF file should match the contig / chromosome ID 
of the local FASTA genome sequence file.

=head1 OPTIONS

fetch_seq_from_ucsc.pl takes the following arguments:

=over 4

=item -h or --help (Optional)

  Displays this helpful message.

=item -q or --quiet (Optional)

  Providing this option suppresses the log messages to the shell.
  Default: disabled

=item -tmap or --tmap (Required)

  Modified TMAP file from Cufflinks pipeline that contains chromosome information or
  a GTF file wherein, for each transcript of the GTF file, intron sequence is spliced out.

=item -d or --dbkey (Optional)

  Database to query (ex: rn4, mm9, hg19 etc...) if sequences should be fetched from UCSC.

=item -ncbi or --ncbi-gi-list (Required)

  In some instances, reference genome sequence for organism of your interest may not be available
  at UCSC. In such instances, you can use this option to retrieve DNA sequence from NCBI Entrez.
  You need to provide a tab-delimited file of GI IDs corresponding to the chromosome information
  of the GTF file.

  Example (Arabidopsis Thaliana):
  -------------------------------

  chr1     240254421
  chr2     330250293
  chr3     332640072
  chr4     240256243
  chr5     240256493

=item -p or --pause-ncbi (Optional)

  If fetching from NCBI, pause for this many seconds between queries to NCBI,
  so that our IP address is not blocked.

  Default: 1 second

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

=item -local or --local-ref-fa

  If you want to fetch sequences from a local reference FASTA instead of from UCSC or
  NCBI, provide the reference FASTA with this option.

  ** Requires Bio::SeqIO module be installed and available **

=back

=head1 AUTHOR

Kranti Konganti, E<lt>konganti@tamu.eduE<gt>.

=head1 COPYRIGHT

This program is distributed under the Artistic License 2.0.

=head1 DATE

Apr-29-2015

=cut
