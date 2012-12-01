#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

my ($CHANGEDBY) = q$LastChangedBy: konganti $ =~ m/.+?\:(.+)/;
my ($LASTCHANGEDDATE) = q$LastChangedDate: 2012-12-01 11:12:22 -0600 (Sat, 01 Dec 2012) $ =~ m/.+?\:(.+)/;
my ($VERSION) = q$LastChangedRevision: 44 $ =~ m/.+?(\d+)/;
my $AUTHORFULLNAME = 'Kranti Konganti';


# Being extra cautious for nice die message instead of Perl's message

check_and_load_modules(['Pod::Usage', 'Getopt::Long', 'File::Basename',
                        'Bio::SearchIO', 'Bio::Tools::GFF',
                        'Bio::SeqFeature::Generic']);
this_script_info();

# Declare initial global variables

my ($quiet, $blast_res_file, $output, $help, $num_hits, $identity_cutoff,
    $coverage_cutoff, $gff_feature_name);
my $is_valid_option = GetOptions ('help|?' => \$help,
                                  'quiet' => \$quiet,
                                  'output=s' => \$output,
                                  'blast-res-file=s' => \$blast_res_file,
                                  'num-hits=i' => \$num_hits,
                                  'identity-cutoff=i' => \$identity_cutoff,
                                  'coverage-cutoff=i' => \$coverage_cutoff
                                  );

# Check for the validity of options
verify_options($is_valid_option, $help);
verify_input_files([$blast_res_file], ['BLAST Result']);

# Set default values
$identity_cutoff = 0 if (!$identity_cutoff);
$coverage_cutoff = 0 if (!$coverage_cutoff);
$num_hits = 1 if (!$num_hits);
$output = $ENV{'PWD'} if (!$output);
$output = validate_create_path($output, 'create', 'Output');

my ($blast_res_file_basename, $path, $suffix) = fileparse($blast_res_file, qr/\.[^.]*/);
my $gff_out_file = $output . $blast_res_file_basename . '.gff';
my $gff_fh = open_file($gff_out_file, '>');

my $gff_out = Bio::Tools::GFF->new(#-file => $gff_out_file,
                                   -fh => $gff_fh,
                                   -gff_version => 3
                                   );

execute_system_command(0,
                       "\nNew GFF File is $gff_out_file ...\n");
execute_system_command(0,
                       "\nParsing BLAST Result from $blast_res_file ...\nand\nConverting to GFF ...\n\n");

my $blast_res_file_obj = new Bio::SearchIO->new(-file => $blast_res_file,
                                                -format => 'blast');

while (my $blast_res = $blast_res_file_obj->next_result) {
    next if ($blast_res->num_hits == 0);

    my $query_name = $blast_res->query_accession;
    my ($cufflinks_class_code, $chr, $chr_start, $chr_end) = ($query_name =~ m/class\_code\:(\w).+?(chr\w+)\:(\d+)-(\d+)/);
    $gff_feature_name = 'intron' if ($cufflinks_class_code =~ m/i/i);
    $gff_feature_name = 'unknown_intergenic' if ($cufflinks_class_code =~ m/u/i);

    my $hit_count = 0;
    while (my $hit = $blast_res->next_hit) {
        next if ($hit->num_hsps == 0);

        my $query_coverage = $hit->frac_aligned_query;
        my $percent_identity = $hit->frac_identical;
        my $evalue = $hit->significance;

        #print "$query_name\t$hit_name\tQL:$query_length\tHitCov:$hit_coverage\nQCOV:$query_coverage\tHI:$percent_identity\n\n";

        $query_coverage *= 100;
        $percent_identity *= 100;

        #print "$query_coverage\t$percent_identity\n";

        if ($query_coverage >= $coverage_cutoff && $percent_identity >= $identity_cutoff) {
            while (my $hsp = $hit->next_hsp) {
                my $corresponding_chr_start = $chr_start + $hsp->start - 1;
                my $corresponding_chr_end = $chr_end + $hsp->end - 1;
                my $gff_feature = new Bio::SeqFeature::Generic(-seq_id => $chr,
                                                               -source_tag => 'BLAST',
                                                               -strand => $hsp->strand('query'),
                                                               -primary => $gff_feature_name,
                                                               -display_name => $hit->name,
                                                               -display_id => $hit->name,
                                                               -seq_name => $hit->name,
                                                               -tag => {desc => $hit->name . $hit->description},
                                                               -score => $hsp->score,
                                                               -start => $corresponding_chr_start,
                                                               -end => $corresponding_chr_end,
                                                               -frame => '.'
                                                               );
                $gff_out->write_feature($gff_feature);

            }
        }

        $hit_count++;
        last if ($hit_count == $num_hits);
    }
}

$blast_res_file_obj->close;
close $gff_fh;

############################### Functions ##################################

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

        # This is not needed when parsing BLAST Result File
        # Removing DOS Carriage returns
        #`perl -i -p -e 's/\r/\n/g' $file`;
        #`perl -i -p -e 's/^\n//' $file`;

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
    print "\n", '@ ', '*' x 78, ' @', "\n";
    print "    Program Name       :  " , basename($0), "\n";
    print "    Version            :  $VERSION\n" if ($VERSION);
    print "    Author             :  $AUTHORFULLNAME\n" if ($AUTHORFULLNAME);
    print "    Last Changed By    : $CHANGEDBY\n" if ($CHANGEDBY);
    print "    Last Changed Date  : $LASTCHANGEDDATE\n";
    print '@ ', '*' x 78, ' @', "\n\n";
    return;
}

__END__

=head1 NAME

blast2gff.pl

=head1 SYNOPSIS

This script will parse the BLAST output and create a GFF file suitable for
upload into UCSC.

****************** !!! YOU HAVE BEEN WARNED !!! ************************

THIS SCRIPT, IN 99% OF THE CASES WILL ONLY WORK WITH BLAST
DATABASES AND QUERY SEQUENCES OBTAINED FROM split_fasta_seqs.pl and
fetch_seq_from_ucsc.pl. ONE CAN ALWAYS MODIFY THE CODE AS
THEY SEE APPROPRIATE TO FIT THEIR NEEDS.

************************************************************************

Examples:

    perl blast2gff.pl -h

    perl blast2gff.pl -b proximal_vs_reference.cmp.transcripts.gtf.class_codes_u_and_i.tmap_vs_ncRNA.blastn

=head1 DESCRIPTION

This script uses BioPerl to parse the text output of BLAST result file and create
GFF format file. Users can filter the features based on query percent coverage and
hit percent identity cut off values.

=head1 OPTIONS

blast2gff.pl takes the following arguments:

=over 4

=item -h or --help (Optional)

  Displays this helpful message.

=item -q or --quiet (Optional)

  Providing this option suppresses the log messages to the shell.
  Default: disabled

=item -b or --blast-res-file (Required)

  The BLAST result file in text format (i.e. default BLAST output).

=item -i or --identity-cutoff (Optional)

  Integer value to filter hits based on percent identity.
  Default: Disabled

=item -c or --coverage-cutoff (Optional)

  Integer value to filter hits based on the query coverage value. (i.e percent
  bases of the query sequence that matched in the database).
  Default: Disabled

=item -n or --num-hits (Optional)

  Integer value to report upto this number of hits for each query.
  Default: 1 hit

=item -o or --output (Optional)

  Path to output directory. If it is not mentioned, the GFF file will be created
  in current working directory.

=back

=head1 AUTHOR

Kranti Konganti, E<lt>konganti@tamu.eduE<gt>.

=head1 COPYRIGHT

This program is distributed under the Artistic License.

=head1 DATE

Nov-30-2012

=cut
