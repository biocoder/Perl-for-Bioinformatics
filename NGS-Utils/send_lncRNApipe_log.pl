#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use IO::Routine;
use MIME::Lite;
use Email::Valid;

my ($LASTCHANGEDBY) = q$LastChangedBy: konganti $ =~ m/.+?\:(.+)/;
my ($LASTCHANGEDDATE) = q$LastChangedDate: 2015-07-08 11:45:27 -0500 (Fri, 08 May 2015)  $ =~ m/.+?\:(.+)/;
my ($VERSION) = q$LastChangedRevision: 0709 $ =~ m/.+?\:\s*(.*)\s*.*/;
my $AUTHORFULLNAME = 'Kranti Konganti';

my ($help, $quiet, $log_file, $user_email, $msg);

my $is_valid_option = GetOptions('help:s'  => \$help,
                                 'quiet'   => \$quiet,
				 'log=s'   => \$log_file,
				 'email=s' => \$user_email,
				 'msg:s'   => \$msg);

my $io = IO::Routine->new($help, $quiet);
my $s_time = $io->start_timer;
$io->verify_options([$is_valid_option, $user_email]);
$io->verify_files([$log_file], ['lncRNApipe log']);

$io->c_time('Looking for lncRNApipe log file ...');
$io->c_time('Verifying Email');
unless(Email::Valid->address( -address => $user_email)) {
    $io->error('Invalid email address: ' . $user_email . ' !!');
}
$io->c_time('Composing email ...');

email("lncRNApipe Run Report from $ENV{'USER'}", 
      $log_file);

$io->c_time('Mail sent to lncRNApipe@outlook.com. Any replies will be sent to ' . $user_email . '.');
$io->end_timer($s_time);

######################### Functions ############################################
sub email {
    my $sbjt = shift;
    my $e_file = shift;

    $msg .= "\n\rThis is an automated message, Please do not reply.\n\n-----\nThanks,\nlncRNApipe\n\n";

    my $email = MIME::Lite->build (
        From     => $user_email,
	To       => 'lncRNApipe@outlook.com',
	Subject  => $sbjt,
        Type     => 'multipart/mixed',
        Debug    => 1
        );

    if ($msg) {
        $email->attach (
            Type => 'TEXT',
            Data => $msg
            );
    }

    if ($e_file) {
        $email->attach (
            Type => 'TEXT',
            Data => "\n\nNote: lncRNApipe log file attached.\n\n"
	    );

        $email->attach (
            Type => 'AUTO',
            Path => $e_file,
            Filename => $io->file_basename($e_file),
            Disposition => 'attachment'
	    ) || $io->error('Cannot attach lncRNApipe log file [ ' . $log_file . ' ] ...');
    }


    $email->send || die "Cannot send mail: $!";
    return;
}

__END__

=head1 NAME

send_lncRNApipe_log.pl

=head1 SYNOPSIS

This script will send lncRNApipe log to the developer.

=head2 DOCUMENTATION

    perldoc send_lncRNApipe_log.pl

=head3 EXAMPLES:

    perl send_lncRNApipe_log.pl -log lncRNApipe.log -email your_email@gmail.com

=head1 DESCRIPTION

This script will attach the specified log file and send email to me so that I can troubleshoot if you are facing any issues
with the pipeline. It is important that you specify your email address so that I can reply with any possible solutions.

=head1 OPTIONS

send_lncRNApipe_log.pl takes the following arguments:

=over 4

=item -h or --help (Optional)

    Displays this helpful message.

=item -q or --quiet (Optional)

    Turn off logging.

=item -e or --email (Required)

    Your email address.
    
    Ex: your_email@gmail.com

=item -l or --log (Required)

    Path to lncRNApipe log file.
    
    Ex: -l /data/lncRNApipe/run.log

=item -m or --msg (Optional)

    Any optional message you want to send with the report.
    Ex: -m 'I cannot run lncRNApipe because, I get the following error'

=back

=head1 AUTHOR

Kranti Konganti, E<lt>konganti@tamu.eduE<gt>.

=head1 COPYRIGHT

This program is distributed under the Artistic License.

=head1 DATE

May-08-2015

=cut
