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
#  $Id: rusage_example.pl,v 1.2 2004/04/27 20:50:35 harsch Exp $
##############################################################################

use strict;
use warnings;
use FindBin qw/ $RealBin /;
use Schedule::DRMAAc qw/ :all /;
use Data::Dumper qw/ Dumper /;

$|++;  # unbuffer STDOUT

use constant JOB_CHUNK => 1;
use constant NBULKS => 1;

sub create_job_template( $ $ $ );
sub print_job_results( $ );

my $job = "$RealBin/sleep_and_echo.csh";
my $rps_vals = {
$DRMAA_PS_UNDETERMINED         	=> 'DRMAA_PS_UNDETERMINED',
$DRMAA_PS_QUEUED_ACTIVE        	=> 'DRMAA_PS_QUEUED_ACTIVE',
$DRMAA_PS_SYSTEM_ON_HOLD       	=> 'DRMAA_PS_SYSTEM_ON_HOLD',
$DRMAA_PS_USER_ON_HOLD         	=> 'DRMAA_PS_USER_ON_HOLD',
$DRMAA_PS_USER_SYSTEM_ON_HOLD  	=> 'DRMAA_PS_USER_SYSTEM_ON_HOLD',
$DRMAA_PS_RUNNING              	=> 'DRMAA_PS_RUNNING',
$DRMAA_PS_SYSTEM_SUSPENDED     	=> 'DRMAA_PS_SYSTEM_SUSPENDED',
$DRMAA_PS_USER_SUSPENDED       	=> 'DRMAA_PS_USER_SUSPENDED',
$DRMAA_PS_USER_SYSTEM_SUSPENDED	=> 'DRMAA_PS_USER_SYSTEM_SUSPENDED',
$DRMAA_PS_DONE                 	=> 'DRMAA_PS_DONE',
$DRMAA_PS_FAILED               	=> 'DRMAA_PS_FAILED',
};


my( $error, $diagnosis ) = drmaa_init( undef );
die drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

warn "*** Send Jobs to the DRM ***";
eval {
	# create the job template
	my $jt = create_job_template( $job, 1, 1 );

	my @all_jobids;
	
	#
	# Submit some bulk jobs
	#
	warn "Run drmaa_run_bulk_jobs ";
	for( my $n=0; $n < &NBULKS; $n++ ) {
		my $jobids;
		while( ( $error, $jobids, $diagnosis ) = drmaa_run_bulk_jobs( $jt, 1, &JOB_CHUNK, 1 ) ) {
			if( $error == $DRMAA_ERRNO_DRM_COMMUNICATION_FAILURE ) {
				warn drmaa_strerror( $error ) . "\n$diagnosis\nCouldn't contact qmaster... retrying: ";
				sleep 1;
				next;
			} # end if
			if( $error != $DRMAA_ERRNO_SUCCESS ) {
				die drmaa_strerror( $error ) . "\n$diagnosis\nError: drmaa_run_bulk_jobs failed.  Aborting... ";
			} # end if
			last;
		} # end while

		print "submitted bulk job with jobids:\n";
		for (my $j=0; $j<JOB_CHUNK; $j++) {
			my( $error, $jobid ) = drmaa_get_next_job_id( $jobids );
			die "\$error='$error' \$jobid='$jobid'" if $error != $DRMAA_ERRNO_SUCCESS;
			push @all_jobids, $jobid;
			print "\t \"$jobid\"\n";
		}
		drmaa_release_job_ids($jobids);
	} # end for

	#
	# Submit some sequential jobs
	#
	$jt = create_job_template( $job, 1, 0 );

	warn "Run drmaa_run_job ";
	for (my $j=0; $j<JOB_CHUNK; $j++) {
		my $jobid;
		while( ( $error, $jobid, $diagnosis ) = drmaa_run_job( $jt ) ) {
			if( $error == $DRMAA_ERRNO_DRM_COMMUNICATION_FAILURE ) {
				warn "Couldn't contact qmaster... retrying:  $error, $diagnosis ";
				sleep 1;
				next;
			} # end if
			if( $error != $DRMAA_ERRNO_SUCCESS ) {
				die "drmaa_run_bulk_jobs failed:  $error, $diagnosis ";
			} # end if
			last;
		} # end while

		print "\t \"$jobid\"\n";
		push @all_jobids, $jobid;
	} # end while

	( $error, $diagnosis ) = drmaa_delete_job_template( $jt );
	die drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

	#
	#   synchronize with all jobs
	#
	# simple synchronize
	( $error, $diagnosis ) = drmaa_synchronize( \@all_jobids, $DRMAA_TIMEOUT_NO_WAIT, 0 );
	die drmaa_strerror( $error ) . "\n" . $diagnosis if $error and $error != $DRMAA_ERRNO_EXIT_TIMEOUT;
	warn("synchronized with all jobs\n");

	#sleep 1;
	#( $error, $diagnosis ) = drmaa_control( $all_jobids[$#all_jobids], $DRMAA_CONTROL_TERMINATE );
	#die drmaa_strerror( $error ) . "\n" . $diagnosis if $error;
	
	print_job_results( \@all_jobids );


}; # end eval

die $@ if $@;


sub create_job_template( $ $ $ ) {
	my $job_path = shift;
	my $seconds = shift;
	my $as_bulk_job = shift;

	my( $error, $jt, $diagnosis ) = drmaa_allocate_job_template();
	die drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

	warn "run in users home directory";
	( $error, $diagnosis ) = drmaa_set_attribute( $jt, $DRMAA_WD, $DRMAA_PLACEHOLDER_HD );
	die drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

	warn "set the job to be run";
	( $error, $diagnosis ) = drmaa_set_attribute( $jt, $DRMAA_REMOTE_COMMAND, $job );
	die drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

	warn "Set job parameters (vector attributes)";
	drmaa_set_vector_attribute( $jt, $DRMAA_V_ARGV, [ "1", "OUTPUT", "ABORT" ] );  # 1st argv param = 1 (to run for 1 second),
															# 2nd argv param = OUTPUT (which getys echo'd to STDOUT)
	warn "join output/error file";
	( $error, $diagnosis ) = drmaa_set_attribute( $jt, $DRMAA_JOIN_FILES, "y" );
	die drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

	warn "path for output";
	if( ! $as_bulk_job ) {
		( $error, $diagnosis ) = drmaa_set_attribute( $jt, $DRMAA_OUTPUT_PATH,
										":" . $DRMAA_PLACEHOLDER_HD . "/DRMAA_JOB" );
	} else {
		( $error, $diagnosis ) = drmaa_set_attribute( $jt, $DRMAA_OUTPUT_PATH,
											":" . $DRMAA_PLACEHOLDER_HD . "/DRMAA_JOB." . $DRMAA_PLACEHOLDER_INCR );
	} # end if
	die drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

	return $jt;
} # end sub

sub print_job_results( $ ) {
	my $all_jobids = shift;

	foreach my $jobid ( @$all_jobids ) {
          my( $error, $job_id_out, $stat, $rusage, $diagnosis ) = drmaa_wait( $jobid, $DRMAA_TIMEOUT_WAIT_FOREVER );
		die drmaa_strerror( $error ) . "\n" . $diagnosis if $error and $error != $DRMAA_ERRNO_NO_RUSAGE;
		
		if( $error != $DRMAA_ERRNO_NO_RUSAGE ) {
			warn 'print \$rusage ' . Dumper( \$rusage );
			while( my( $retval, $value ) = drmaa_get_next_attr_value( $rusage ) ) {
				last if $retval != $DRMAA_ERRNO_SUCCESS;
	
				print "\"value\" => \"$value\"\n";
			} # end while
			drmaa_release_attr_values( $rusage );
		} # end if

		my $aborted;
		( $error, $aborted, $diagnosis ) = drmaa_wifaborted( $stat );
		die drmaa_strerror( $error ) . "\n" . $diagnosis if $error;
	
		if( $aborted ) {
			print "job \"$jobid\" never ran\n";
		} else {
			my $exited;
			( $error, $exited, $diagnosis ) = drmaa_wifexited( $stat );
			die drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

			if( $exited ) {
				my $exit_status;
				( $error, $exit_status, $diagnosis ) = drmaa_wexitstatus( $stat );
				die drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

				print "job \"$jobid\" finished regularly with exit status $exit_status\n";
			} else {
				my $signaled;
				( $error, $signaled, $diagnosis ) = drmaa_wifsignaled( $stat );
				die drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

				if( $signaled ) {
					my $termsig;
					( $error, $termsig, $diagnosis ) = drmaa_wtermsig( $stat );
					die drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

					print "job \"$jobid\" finished due to signal \"$termsig\".\n";
				} else {
					print "job \"$jobid\" finished with unclear conditions\n";
				} # end if( $signaled )
			} # end if( $exited )
		} # end if( $aborted )
	} # end foreach

} # end sub

1; # Ancient Druid Custom
