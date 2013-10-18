package IO::Routine;

use 5.14.2;
use strict;
use warnings FATAL => 'all';
use Carp;
use Pod::Usage;
use Module::Load;
use Exporter;
use Time::HiRes qw(gettimeofday tv_interval);
use File::Basename;
use Cwd;

=head1 NAME

=over 3

IO::Routine - An attempt to provide a solution to avoid routine IO chores.

=back

=head1 VERSION

=over 3

Version 0.31

=back

=cut

our $VERSION = '0.31';

=head1 SYNOPSIS

=over 3

To avoid some of the repetetive stuff, this module:

=back

=over 4

=item * Checks file / directory paths.

=item * Verifies command line options used with Getopt module.

=item * Executes system level commands, while printing message to STDOUT.

=item * Checks for minimum version requirement for GNU Core Utils.

=item * Prints SVN version control information ( If requested ).

=item * Gets current memory usage of the program ( in GB(s) for Linux ).

=item * Starts timer to calculate elapsed time, prints time stamp.

=back

=head1 USAGE

=over 4

    use IO::Routine;
    use Getopt::Long;

    my $is_valid_option = $io->GetOptions('help|?' => \$help,
                                          'quiet' => \$quiet,
                                          'output=s' => \$output,
                                          'file1=s' => \$file1,
                                          'file2=s' => \$file2
                                         );

    # Get the new IO::Routine object. Passing quiet and help while instantiating the
    # the object sets the precedent for verbosity as requested by the user during execution.
    my $io = IO::Routine->new($help, $quiet);

    # If you are maitanining POD in separate file. It will search for POD file at the
    # current location of script file.
    my $io = IO::Routine->new($help, $quiet, 'podfilename');

    # Start the timer
    my $s_time = $io->start_timer();

    # Verify options from Getopt::Long and display POD's OPTIONS section
    # on die.
    $io->verify_options([$is_valid_option, $file1]);
    
    # Verify existense of non-empty files.
    $io->verify_files([$file2, $file2], ['file1 msg', 'file2 msg'])

    # Validate the output directory and make one if asked to.
    my $output = $io->validate_create_path($output, 'create', 'Output directory');
    my $output = $io->validate_create_path($output, 'do not create', 'Input directory');

    # With 0 as first argument, just display log message.
    $io->execute_system_command(0,
                                'Checking the version requirement for system level commands ...');
    
    # With a system command as first option, execute system command while displaying log message.
    $io->execute_system_command('ls',
                                'List of files ...');
    

    # Execute system command and get its output into a variable.
    my $unix_date = $io->execute_get_sys_cmd_output("/bin/date",
                                                        "\nGetting system date ... \n");

    # Check minimum version requirements of GNU Core Utils.
    $io->check_sys_level_cmds(['grep', 'ls', 'sort', 'basename'],
                              ['2.10', '8.13', '8.13', '8.13']);

    # On UNIX based systems, get memory usage. Useful while processing a large data structure.
    # You may refresh the command buffer (Ex: $|) to get up-to-date details.
    my ($v_mem, $r_mem) = $io->get_mem_usage();


    # End timer and print time elapsed in seconds, minutes, hours or days.
    $io->end_timer($s_time);

    # Print logging message with local time stamp.
    $io->c_time("\nAnalysis Finished on ");
    
    # Get different parts of the file name.
    my $filename = $io->file_basename($file);
    my $filename_w_suffix = $io->file_basename($file, 'suffix');
    my ($filename, $path, $suffix) = $io->file_basename($file, 'all');

    # Check if any system command exists
    $io->exist_sys_cmd(['cuffcompare -v', 'cufflinks', 'bwa']);

=back

=head1 SUBROUTINES/METHODS

=over 4

=item new()

=back

=over 5

Instantiates a new IO::Routine object and returns it.

=back

=over 4

=item file_basename()

=back

=over 5

Returns file's basename or suffix or both if requested

=back

=over 4

=item strip_leading_and_trailing_spaces()

=back

=over 5

Strips leading and traling white spaces.

=back

=over 4

=item verify_options()

=back

=over 5

Verifies options declared with Getopt::Long.

=back

=over 4

=item verify_files()

=back

=over 5

Verifies the existence of input files.

=back

=over 4

=item is_dir_empty()

=back

=over 5

Checks for existence of directory and whether or not it holds some files.

=back

=over 4

=item validate_create_path()

=back

=over 5

Validates and creates output path, if asked.

=back

=over 4

=item execute_system_command()

=back

=over 5

Method to enable / disable verbosity of various script actions, while executing system level commands.

=back

=over 4

=item execute_get_sys_cmd_output()

=back

=over 5

Execute system level command and capture everything from STDERR and STDOUT.

=back

=over 4

=item check_sys_level_cmds()

=back

=over 5

Checks the minimum version requirement of GNU Core Utils.

=back

=over 4

=item error()

=back

=over 5

Prints error with msg to STDERR and quits.

=back

=over 4

=item warning()

=back

=over 5

Prints warning message to STDOUT and continues.

=back

=over 4

=item open_file()

=back

=over 5

Opens file in requested mode and returns file handle.

=back

=over 4

=item this_script_info()

=back

=over 5

Prints Program's name, version control information (Supported for only SVN).

=back

=over 4

=item get_mem_usage()

=back

=over 5

Get memory usage on LINUX systems ( in GB(s) ) based on PID from /proc/PID and using vmmap utility for OS X systems.  Useful in data structure loops. You need to specify $|++; to flush buffer.

=back

=over 4

=item start_timer()

=back

=over 5

Start timer to calculate elapsed time and return reference.

=back

=over 4

=item end_timer()

=back

=over 5

End timer and print elapsed time in floating point seconds.

=back

=over 4

=item c_time()

=back

=over 5

Return ctime using localtime function.

=back

=over 4

=item exist_sys_cmd()

=back

=over 5

Abort if system level command does not exist.

=back

=head1 AUTHOR

=over 3

Kranti Konganti, E<lt>konganti@tamu.eduE<gt>.

=back

=head1 BUGS

=over 3

Please report any bugs or feature requests to C<bug-io-routine at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=IO-Routine>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.

=back

=head1 SUPPORT

=over 3

You can find documentation for this module with the perldoc command.

    perldoc IO::Routine


You can also look for information at:

=back

=over 5

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=IO-Routine>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/IO-Routine>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/IO-Routine>

=item * Search CPAN

L<http://search.cpan.org/dist/IO-Routine/>

=back

=head1 LICENSE AND COPYRIGHT

=over 3

Copyright 2013 Kranti Konganti.

This program is free software; you can redistribute it and/or modify it
under the terms of the the Artistic License (2.0). You may obtain a
copy of the full license at:

L<http://www.perlfoundation.org/artistic_license_2_0>

Any use, modification, and distribution of the Standard or Modified
Versions is governed by this Artistic License. By using, modifying or
distributing the Package, you accept this license. Do not use, modify,
or distribute the Package, if you do not accept this license.

If your Modified Version has been derived from a Modified Version made
by someone other than you, you are nevertheless required to ensure that
your Modified Version complies with the requirements of this license.

This license does not grant you the right to use any trademark, service
mark, tradename, or logo of the Copyright Holder.

This license includes the non-exclusive, worldwide, free-of-charge
patent license to make, have made, use, offer to sell, sell, import and
otherwise transfer the Package with respect to any patent claims
licensable by the Copyright Holder that are necessarily infringed by the
Package. If you institute patent litigation (including a cross-claim or
counterclaim) against any party alleging that the Package constitutes
direct or contributory patent infringement, then this Artistic License
to you shall terminate on the date that such litigation is filed.

Disclaimer of Warranty: THE PACKAGE IS PROVIDED BY THE COPYRIGHT HOLDER
AND CONTRIBUTORS "AS IS' AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES.
THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE, OR NON-INFRINGEMENT ARE DISCLAIMED TO THE EXTENT PERMITTED BY
YOUR LOCAL LAW. UNLESS REQUIRED BY LAW, NO COPYRIGHT HOLDER OR
CONTRIBUTOR WILL BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, OR
CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE OF THE PACKAGE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=back

=cut


# Will deal when I have figured out how to avoid code redundancy regarding
# check_and_load_modules

#BEGIN {
#  use vars qw[@ISA @EXPORT_OK];
#  @ISA = qw[Exporter];
#  @EXPORT_OK = qw[check_and_load_modules];
#}

# IO::Routine Constructor

my ($thisHelp, $thisQuiet, $podFilename);
$thisHelp = $thisQuiet = $podFilename = 0;

sub new {
    my $class = shift;
    $thisHelp = shift;
    $thisQuiet = shift;
    $podFilename = shift;

    my $self = {};
    bless $self, $class;

    $podFilename = file_basename($self, $podFilename) . '.pd' if ($podFilename);

    return $self;
}


############## Will deal with this block of code later ########################
#
# Check and load modules.
#
#sub check_and_load_modules {
#    my $self = shift;
#    my $module_list = shift;
#    my $is_module_loadable = 1;
#    my $req_modules = '';
#
#    foreach my $module (@$module_list) {
#      #eval { load $module };
#      eval "use $module; 1";
#      $is_module_loadable = 0,
#        $req_modules .= "$module, "
#	  if ( $@ );
#    }
#    $req_modules =~ s/\,\s+$//;
#
#    confess error($self,
#		  "Cannot load requested module(s):\n$req_modules\n\nMake sure that the Perl modules are in your path.\n\nHint: Use PERL5LIB environment #variable.\n")
#        if (!$is_module_loadable);
#
#    return $self;
#}


# Check if all options entered by user are valid with Getopt.

sub verify_options {
    my $self = shift;
    my $valid_options = shift;
    my $help = shift;

    $help = $thisHelp if (!$help);

    if (!$podFilename) {
	if ($help) {
	    pod2usage(-exitval => 1,
		      -sections => "OPTIONS",
		      -msg => "\n");
	}
	
	foreach my $valid_option (@$valid_options) {
	    if (!$valid_options || !defined($valid_option)) {
		pod2usage(-exitval => 2,
			  -verbose => 2,
			  -msg => "\nSee $0 -h for options.\n");
	    }
	}
    }
    elsif ($podFilename &&
	   -e $podFilename &&
	   -s $podFilename != 0) {
	if ($help) {
            pod2usage(-exitval => 1,
                      -sections => "OPTIONS",
                      -msg => "\n",
		      -input => $podFilename);
        }
	
	foreach my $valid_option (@$valid_options) {
            if (!$valid_options || !defined($valid_option)) {
                pod2usage(-exitval => 2,
                          -verbose => 2,
                          -msg => "\nSee $0 -h for options.\n",
			  -input => $podFilename);
            }
        }
    }
    elsif ($podFilename &&
           !-e $podFilename) {
	confess error($self,"POD file [ $podFilename ] does not exist: $!");
    }
    elsif ($podFilename &&
	   -s $podFilename == 0) {
	confess error($self,"POD file [ $podFilename ] is empty: $!");
    }

    return;
}


# Check the existence of files.

sub verify_files {
    my $self = shift;
    my $files = shift;
    my $what_files = shift;
    my $file_no = 0;

    foreach my $file (@$files) {

	confess error($self, "Empty argument!\nPath to file not provided?")
	    if (!$file || $file eq '');

        confess error($self, "@$what_files[$file_no] file not specified: $!")
            if (!defined $file);

        confess error($self, "@$what_files[$file_no] file ( $file ) does not exist: $!")
            if (!-e $file);

        confess error($self, "@$what_files[$file_no] file ($file) is empty: $!")
            if (-s $file == 0);

        $file_no++;
    }
    return;
}

# Check if the directory is empty

sub is_dir_empty {
    my $self = shift;
    my $dirname = shift;

    confess error("Directory $dirname does not exist")
        if (!-d $dirname);

    opendir (my $dh, $dirname) ||
        confess error($self, "$dirname not a directory: $!");

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
        return 1;
    }
    else {
        return 0;
    }
}

# Check if output path is mentioned, else create output path, if asked with 'create' message.

sub validate_create_path {
    my $self = shift;
    my $path = shift;
    my $create_dir = shift;
    my $msg = shift;

    confess error($self, "Invalid number of options for validate_create_path!\nLog Message not provided...")
	if (!$msg || $msg eq '');

    confess error ($self, "$msg path not defined or entered!")
        if (!defined $path || $path eq '');

    if ($create_dir eq 'create') {
        if (defined $path) {
            execute_system_command($self,
                                   "mkdir -p $path",
                                   "Attempting to make directory $path ...",
                                   )
                if (defined ($path) && !-d $path);
        }
        else {
            $path = $ENV{'PWD'};
        }
    }

    confess error($self, "Path ( $path ) does not exist: $!")
        if (!-d $path);

    $path .= '/'
        if ($path !~ m/\/$/);

    return $path;
}

# Subroutine to execute system command with verbosity.

sub execute_system_command {
    my $self = shift;
    my $command = shift;
    my $print_msg = shift;
    my $quiet = shift;

    $quiet = $thisQuiet if (!$quiet);

    if (!$quiet) {
       print "\n$print_msg\n\n" if ($print_msg);
       system ("$command") if ($command);
    }
    elsif ($quiet) {
        system ("$command 1>/dev/null 2>/dev/null") if ($command);
    }
    return;
}

# Subroutine to execute and return system command's output

sub execute_get_sys_cmd_output {
    my $self = shift;
    my $command = shift;
    my $print_msg = shift;
    my $quiet = shift;

    $quiet = $thisQuiet if (!$quiet);

    if (!$quiet) {
       print "\n$print_msg\n\n" if ($print_msg);
     }

    my $cmd_out = qx($command 2>&1);
    $cmd_out = 'Could not capture stream from either STDERR or STDOUT' if (!$cmd_out || $cmd_out eq '');

    return $cmd_out;
}

# Shell msg that differentiates log from error

sub error {
    my $self = shift;
    my $msg = shift;
    print STDERR "\nERROR!\n------\n$msg\n\n";
    pod2usage(-exitval => 2,
	      -verbose => 2);
    return;
}


# Shell msg for warning

sub warning {
    my $self = shift;
    my $msg = shift;
    print STDOUT "\nWARNING!\n--------\n$msg\n\n";
    return;
}

# Subroutine to open files and return file handle

sub open_file {
    my $self = shift;
    my $mode = shift;
    my $file = shift;

    confess error($self, "Empty argument!\nPath to file not provided?")
        if (!$file || $file eq '');

    return getcwd() if ($mode =~ m/^cwd$/i);
    return *STDOUT if ($mode =~ m/^stdout$/i);

    my $tainted_filename;

    if ($file =~ m/(.+?)\z/i) {
      $tainted_filename = $1;
    }

    open (my $file_handle, $mode, $tainted_filename) ||
        confess error($self, "Cannot open ( $file ) in mode ( $mode ): $!");
    return $file_handle;
 }


# Subroutine to print SCRIPT Version. Print this script's info

sub this_script_info {
    my ($self, $PRGNAME, $VERSION, $AUTHORFULLNAME, $CHANGEDBY, $LASTCHANGEDDATE, $quiet) = @_;

    $quiet = $thisQuiet if (!$quiet);
    return if $quiet;

    print "\n", '@ ', '*' x 78, ' @', "\n";
    print "    Program Name       :  " , $PRGNAME, "\n";
    print "    Version            :  $VERSION\n" if ($VERSION);
    print "    Author             :  $AUTHORFULLNAME\n" if ($AUTHORFULLNAME);
    print "    Last Changed By    : $CHANGEDBY\n" if ($CHANGEDBY);
    print "    Last Changed Date  : $LASTCHANGEDDATE\n";
    print '@ ', '*' x 78, ' @', "\n\n";
    
    return;
}

# Subroutine to check the least required version of system command.

sub check_sys_level_cmds {
    my $self = shift;
    my $cmds = shift;
    my $versions = shift;
    my ($curr_version_unf, $req_version);
    my $loop_limiter = 0;

    confess error($self, "Unmatched number of commands and versions list")
	if ($#$versions != $#$cmds);
    
    
    for (my $i=0; $i<scalar(@$cmds); $i++) {
        my $cmd_version_out = `@$cmds[$i] --version`;
        $req_version = @$versions[$i];
	my $req_version_parts = [split(/\./, @$versions[$i])];
                
	if ($cmd_version_out =~ m/^@$cmds[$i].*?(\d+[\.\d]*)/ ||
            $cmd_version_out =~ m/^.*?(\d+[\.\d]*)/) {
            my $curr_version = $curr_version_unf = $1;
	    my $curr_version_parts = [split(/\./, $curr_version)];
	    if ($#$curr_version_parts <= $#$req_version_parts) {
		$loop_limiter = $#$curr_version_parts 
	    }
	    else {
		$loop_limiter = $#$req_version_parts;
	    }
	    
	    for (my $j=0; $j<=$loop_limiter; $j++) {
		my $ver_diff = @$curr_version_parts[$j] - @$req_version_parts[$j];
		if ($ver_diff < 0) {
		    confess error($self,
				  "At least Version @$versions[$i] required for system level command: @$cmds[$i]\nCurrent Version: $curr_version_unf\n");
		}
	    }
	}
	else {
	    confess error($self,
			  "Supplied command [ @$cmds[$i] ] does not seem to be a GNU core utility.\n");
	}
    }
    return;
}

# Subroutine to strip leading and trailing white spaces from a line

sub strip_leading_and_trailing_spaces {
    my $self = shift;
    my $line = shift;
    $line =~ s/^\s+//;
    $line =~ s/\s+$//;
    return $line;
}

# Subroutine returns memory stats based on /proc/ on Linux machines

sub get_mem_usage {
    my $self = shift;
    chomp(my $host_machine = `uname`);
    if ($host_machine =~ m/darwin/i) {
        chomp(my $malloc =`vmmap -resident $$ | grep -i 'see malloc zone table below'`);
        my @mem_details = split/\s+/, $malloc;
        return ('Unable to get vmem', 'Unable to get rmem') if (!$mem_details[1] || !$mem_details[2]);
        return ($mem_details[1], $mem_details[2]);
    }
    elsif ($host_machine =~ m/linux/i) {
        return ('Unable to get vmem', 'Unable to get rmem') if (!-e "/proc/$$/status");
        my $curr_v_mem = `grep VmPeak /proc/$$/status`;
        my $curr_r_mem = `grep VmSize /proc/$$/status`;
        my ($curr_v_mem_usage) = $curr_v_mem =~ m/(\d+)/;
        my ($curr_r_mem_usage) = $curr_r_mem =~ m/(\d+)/;
        $curr_v_mem_usage = ( ($curr_v_mem_usage / 1024) / 1024 );
        $curr_r_mem_usage = ( ($curr_r_mem_usage / 1024) / 1024 );
        $curr_v_mem_usage = sprintf("%.2f", $curr_v_mem_usage) . 'GB';
        $curr_r_mem_usage = sprintf("%.2f", $curr_r_mem_usage) . 'GB';
        return ('Unable to get vmem', 'Unable to get rmem') if (!$curr_v_mem_usage || !$curr_r_mem_usage);
        return ($curr_v_mem_usage, $curr_r_mem_usage);
    }
    return;
}

# Subroutine to start timer

sub start_timer {
    my $self = shift;
    my $start_time = [gettimeofday()];
    return $start_time;
}

# Subroutine to end timer

sub end_timer {
    my $self = shift;
    my $start_time = shift;
    my $quiet = shift;

    $quiet = $thisQuiet if (!$quiet);

    confess error($self, "Unable to calculate time elapsed.\nDid not get start timer object...\n")
	if (!$start_time);
    
    if (!$quiet || !defined($quiet)) {
	if (sprintf("%.2f", tv_interval($start_time)) > 60) {
	    print "\nTime Elapsed: ", sprintf("%.2f", tv_interval($start_time) / 60), " Minute(s).\n\n";
	}
	elsif (( sprintf("%.2f", tv_interval($start_time)) / 60 ) > 60) {
	    print "\nTime Elapsed: ", sprintf("%.2f", tv_interval($start_time) / 3600), " Hour(s).\n\n";
	}
	elsif ((( sprintf("%.2f", tv_interval($start_time)) / 60 ) / 60 ) > 24) {
	    print "\nTime Elapsed: ", sprintf("%.2f", tv_interval($start_time) / 86400), " Day(s).\n\n";
	}
	elsif (sprintf("%.2f", tv_interval($start_time)) <= 60) {
	    print "\nTime Elapsed: ", sprintf("%.2f", tv_interval($start_time)), " Seconds.\n\n";
	}
    }
    return;
}

# Subroutine return current ctime

sub c_time {
    my $self = shift;
    my $msg = shift;
    my $quiet = shift;

    $quiet = $thisQuiet if (!$quiet);

    if (!$quiet || !defined($quiet) && ($msg ne '')) {
	print "\n", scalar(localtime(time)), "\t$msg\n";
    }
    elsif (!$quiet) {
	print "\n", scalar(localtime(time));
    }
    return;
}

# Subroutine to return filename parts

sub file_basename {
    my $self = shift;
    my $file = shift;
    my $mode = shift;
    
    my @file_attrs = fileparse($file, qr/\.[^.]*/);

    return $file_attrs[0] if (!defined($mode));
    return "$file_attrs[0]$file_attrs[2]" if ($mode eq 'suffix');
    return ($file_attrs[0], $file_attrs[1], $file_attrs[2]) if ($mode eq 'all');
}

# Subroutine to check the existence of a system level command

sub exist_sys_cmd {
    my $self = shift;
    my $cmds = shift;

    foreach my $cmd (@$cmds) {
	my $cmd_out = execute_get_sys_cmd_output(0, $cmd);
	error($self, "Cannot find the system level command [ $cmd  ]! ... Aborting.")
	    if ($cmd_out !~ m/$cmd$/);
    }
    return;
}

1; # End of IO::Routine
