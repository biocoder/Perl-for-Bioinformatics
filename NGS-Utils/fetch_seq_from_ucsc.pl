#!/opt/perl/bin/perl

use strict;
use warnings;
use Carp;

my ($CHANGEDBY) = q$LastChangedBy: konganti $ =~ m/.+?\:(.+)/;
my ($LASTCHANGEDDATE) = q$LastChangedDate: 2013-07-29 10:34:50 -0500 (Mon, 29 Jul 2013) $ =~ m/.+?\:(.+)/;
my ($VERSION) = q$LastChangedRevision: 48 $ =~ m/.+?(\d+)/;
my $AUTHORFULLNAME = 'Kranti Konganti';


# Being extra cautious for nice die message instead of Perl's message

check_and_load_modules(['Pod::Usage', 'Getopt::Long', 'File::Basename',
                        'LWP::Simple', 'XML::XPath', 'XML::XPath::XMLParser']);
this_script_info();

# Declare initial global variables

my ($quiet, $tmap, $output, $help, $dbkey, $print_seq_fh);
my $is_valid_option = GetOptions ('help|?' => \$help,
                                  'quiet' => \$quiet,
                                  'output=s' => \$output,
                                  'tmap=s' => \$tmap,
                                  'dbkey=s' => \$dbkey
                                  );

# Check for the validity of options
verify_options($is_valid_option, $help);
verify_input_files([$tmap],
		   ['TMAP ( ex: proximal_vs_reference.cmp.transcripts.gtf.class_codes_u_and_i.tmap )']);

error("DBKEY not mentioned.\nThis option is required to fetch sequences from UCSC Database.\nEx: rn4, hg19 etc..\n")
    if (!defined $dbkey);

$output = validate_create_path($output, 'create', 'Output');

execute_system_command(0,
                       "\nChromosome files will be stored at $output ...\n");

execute_system_command(0,
                       "\nFetching Sequences ...\n");

my $tmap_fh = open_file($tmap, '<');
my ($tmap_filename, $tmap_filepath, $suffix) = fileparse($tmap, qr/\.[^.]*/);
my $fa = $output . $tmap_filename . '.fa';

if (!-e $fa && !-s $fa) {
    $print_seq_fh = open_file($fa, '>')
}
else {
    $print_seq_fh = open_file($fa, '>>');
}

while (my $line = <$tmap_fh>) {
    chomp $line;
    next if ($line =~ m/^ref\_gene\_id/);
    my (@cols) = split/\t/, $line;
    error("TMAP file does not contain a total 16 columns\nCurrent column count $#cols + 1\n")
        if ($#cols != 15);
    my $unique_seq_id = $cols[4] . '_' . "class_code:$cols[2]" . '_' . $cols[13] . ':' . $cols[14] . '-' . $cols[15];
    $cols[13] =~ s/Chr/chr/i;

    my $fetched_seq = `grep -A 1 $unique_seq_id $fa`;
    chomp $fetched_seq;

    if ($fetched_seq || $fetched_seq ne '') {
        execute_system_command(0,
                               "Skipping $unique_seq_id ... Already fetched!\n");
        next;
    }
    else {
	execute_system_command(0,
                               "Querying $unique_seq_id against UCSC DAS ...\n");
    }

    my $xml = get("http://genome.ucsc.edu/cgi-bin/das/$dbkey/dna?segment=$cols[13]:$cols[14],$cols[15]");
    my $xpath = XML::XPath->new(xml=>$xml);
    my $xml_nodes = $xpath->find('/DASDNA/SEQUENCE/DNA/text()');
    print_seq(\$xml_nodes, $unique_seq_id);
}

close $tmap_fh;
close $print_seq_fh;

###################### Functions ###############################################

# Print fasta sequences;
sub print_seq {
    my $xml_nodes = shift;
    my $seq_id = shift;

    foreach my $xml_node ($$xml_nodes->get_nodelist) {
        my $seq = $xml_node->getValue;
        $seq = strip_white_spaces_and_newlines($seq);
        print $print_seq_fh ">$seq_id\n$seq\n";
    }
    return;
}

# Strip white spaces
sub strip_white_spaces_and_newlines {
    my $line = shift;
    $line =~ s/^\s+//;
    $line =~ s/\s+$//;
    $line = uc($line); # upper case all nucleotides
    $line =~ s/[\n\r]*//g;
    return $line;
}

# To check and load modules.
sub check_and_load_modules {
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

# Check if all options entered by user are valid
sub verify_options {
    my $valid_options = shift;
    my $help = shift;

    if (!$valid_options) {
        pod2usage(-exitval => 2,
		  -verbose => 2,
                  -msg => "\nSee $0 -h for options.\n");
    }

    if ($help) {
        pod2usage(-exitval => 1,
                  -sections => "OPTIONS",
                  -msg => "\n");
    }
    return;
}

# File checks
sub verify_input_files {
    my $files = shift;
    my $what_files = shift;
    my $file_no = 0;

    foreach my $file (@$files) {

        confess error("@$what_files[$file_no] file not specified: $!")
            if (!defined $file);

        confess error("@$what_files[$file_no] file ( $file ) does not exist: $!")
            if (!-e $file);

        confess error("@$what_files[$file_no] file ($file) is empty: $!")
            if (-s $file == 0);

        # Removing DOS Carriage returns
        `perl -i -p -e 's/\r/\n/g' $file`;
        `perl -i -p -e 's/^\n//' $file`;

        $file_no++;
    }
    return;
}

# Check if output path is mentioned, else create output path
sub validate_create_path {
    my $path = shift;
    my $create_dir = shift;
    my $msg = shift;

    confess error ("$msg path not defined or entered!")
        if (!defined $path || $path eq '');

    if ($create_dir eq 'create') {
        if (defined $path) {
            execute_system_command("mkdir -p $path",
                                   "\nAttempting to create $path ...\n\n",
                                   )
                if (defined ($path) && !-d $path);
        }
        else {
            $path = $ENV{'PWD'};
        }
    }

    confess error("Path ( $path ) does not exist: $!")
        if (!-d $path);

    $path .= '/'
        if ($path !~ m/\/$/);

    return $path;
}

# Subroutine to execute system command with verbosity
sub execute_system_command {
    my $command = shift;
    my $print_msg = shift;

    if (!$quiet) {
       print $print_msg if ($print_msg);
       system ("$command") if ($command);
    }
    elsif ($quiet) {
        system ("$command 1>/dev/null 2>/dev/null") if ($command);
    }
    return;
}

# Shell msg that differentiates log from error
sub error {
    my $msg = shift;
    print "\nERROR!\n------\n$msg\n\n";
    pod2usage(-exitval => 2,
	      -verbose => 2);
}

# Shell msg for warning
sub warning {
    my $msg = shift;
    print "\nWARNING!\n--------\n$msg\n\n";
    return;
}

# Subroutine to open files and return file handle
sub open_file {
    my $file = shift;
    my $mode = shift;
    open (my $file_handle, "$mode", $file) ||
        confess error("Cannot open ( $file ) in mode ( $mode ): $!");
    return $file_handle;
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

=back

=head1 AUTHOR

Kranti Konganti, E<lt>konganti@tamu.eduE<gt>.

=head1 COPYRIGHT

This program is distributed under the Artistic License.

=head1 DATE

Nov-27-2012

=cut
