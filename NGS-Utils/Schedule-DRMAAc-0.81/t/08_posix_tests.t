##############################################################################
#  Copyright (c) 2004, The Regents of the University of California.
#  Produced at the Lawrence Livermore National Laboratory.
#  Written by Tim Harsch <harsch1@llnl.gov>
#  UCRL-CODE-155918
#  All rights reserved.
#
#  This file is part of Schedule::DRMAAc. For details, see CPAN
#  Please also read LICENSE.txt which is found in this source distribution.
#
#  This program is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License (as published by the
#  Free Software Foundation) version 2, dated June 1991.
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY
#  OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the terms and conditions of the GNU General Public License for more
#  details.
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
##############################################################################
#  $Id: 08_posix_tests.t,v 1.2 2004/04/30 17:10:04 harsch Exp $
##############################################################################

use Schedule::DRMAAc qw/ :all /;

use strict;

use constant TESTS => 32;
use Test::More tests => TESTS;

use Cwd;

sub mydiag( $ $ ; $ );

my $remote_cmd = "/bin/sleep";

my( $error, $diagnosis );
my @vals = drmaa_strerror( 0 );
ok( @vals == 1, "drmaa_strerror returned " . scalar(@vals) . " of 1 args" );

@vals = ( $error, $diagnosis ) = drmaa_init( undef );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_init error?" )
	or mydiag( $error, $diagnosis );

my $jt;
@vals = ( $error, $jt, $diagnosis ) = drmaa_allocate_job_template();
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_allocate_job_template error?" )
	or mydiag( $error, $diagnosis );

@vals = ( $error, $diagnosis ) = drmaa_set_attribute( $jt, $DRMAA_REMOTE_COMMAND, $remote_cmd );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_set_attribute error?" )
	or mydiag( $error, $diagnosis );
( $error, $diagnosis ) = drmaa_set_attribute( $jt, $DRMAA_OUTPUT_PATH, ":$DRMAA_PLACEHOLDER_HD/DRMAA_posix_tests.out" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_set_attribute error?" )
	or mydiag( $error, $diagnosis );
( $error, $diagnosis ) = drmaa_set_attribute( $jt, $DRMAA_ERROR_PATH, ":$DRMAA_PLACEHOLDER_HD/DRMAA_posix_tests.err" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_set_attribute error?" )
	or mydiag( $error, $diagnosis );
( $error, $diagnosis ) = drmaa_set_vector_attribute( $jt, $DRMAA_V_ARGV, [ "1" ] );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_set_vector_attribute error?" )
	or mydiag( $error, $diagnosis );


my $jobid;
@vals = ( $error, $jobid, $diagnosis ) = drmaa_run_job( $jt );
ok( @vals == 3, "drmaa_run_job returned " . scalar(@vals) . " of 3 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_run_job error?" )
	or mydiag( $error, $diagnosis );

# this script tests that a job running a normal job such as /bin/sleep would return values
#  as expected in that case, from the POSIX functions (drmaa_w*).  However, the posix functions
#  are best used to determine the status of problem jobs.
my( $jobid_out, $stat, $rusage );
( $error, $jobid_out, $stat, $rusage, $diagnosis ) = @vals = drmaa_wait( $jobid, $DRMAA_TIMEOUT_WAIT_FOREVER );
ok( @vals == 5, "drmaa_wait returned " . scalar(@vals) . " of 5 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_wait error?" )
	or mydiag( $error, $diagnosis );
# job id should not change unless migrated, but why would it?
is( $jobid, $jobid_out, "drmaa_wait says jobid did not change?" );
ok( $stat == 1, "drmaa_wait should say there is more info available in POSIX funcs" );

my $aborted;
@vals = ( $error, $aborted, $diagnosis ) = drmaa_wifaborted( $stat );
ok( @vals == 3, "drmaa_wifaborted returned " . scalar(@vals) . " of 3 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_wifaborted error?" )
	or mydiag( $error, $diagnosis );
ok( $aborted == 0, "normal job should not abort." );

my $exited;
( $error, $exited, $diagnosis ) = drmaa_wifexited( $stat );
ok( @vals == 3, "drmaa_wifexited returned " . scalar(@vals) . " of 3 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_wifexited error?" )
	or mydiag( $error, $diagnosis );
ok( $exited == 1, "normal job should exit." );

my $exit_status;
( $error, $exit_status, $diagnosis ) = drmaa_wexitstatus( $stat );
ok( @vals == 3, "drmaa_wexitstatus returned " . scalar(@vals) . " of 3 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_wexitstatus error?" )
	or mydiag( $error, $diagnosis );
ok( $exit_status == 0, "normal job should have exit status of 0." );

my $signaled;
@vals = ( $error, $signaled, $diagnosis ) = drmaa_wifsignaled( $stat );
ok( @vals == 3, "drmaa_wifsignaled returned " . scalar(@vals) . " of 3 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_wifsignaled error?" )
	or mydiag( $error, $diagnosis );
ok( $exit_status == 0, "normal job should have not have terminatied due to signal." );
	
my $termsig;
@vals = ( $error, $termsig, $diagnosis ) = drmaa_wtermsig( $stat );
ok( @vals == 3, "drmaa_wtermsig returned " . scalar(@vals) . " of 3 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_wtermsig error?" )
	or mydiag( $error, $diagnosis );

my $core_dumped;
@vals = ( $error, $core_dumped, $diagnosis ) = drmaa_wcoredump( $stat );
ok( @vals == 3, "drmaa_wcoredumped returned " . scalar(@vals) . " of 3 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_wcoredumped error?" )
	or mydiag( $error, $diagnosis );
ok( $core_dumped == 0, "should be no core image of job." );

@vals = ( $error, $diagnosis ) = drmaa_delete_job_template( $jt );
ok( @vals == 2, "drmaa_delete_job_template returned " . scalar(@vals) . " of 2 args" );

( $error, $diagnosis ) = drmaa_exit();
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_exit error?" );


sub mydiag( $ $ ; $ ) {
	my $msg = "error: " . drmaa_strerror( $_[0] ) . "\n" . "diagnosis: " . $_[1];
	$msg .= $_[2] if defined $_[2];
	diag( $msg );
} # end sub

1; # Ancient Druid Custom
