#!/usr/local/bin/perl32

use strict;
use warnings;
use Data::Dumper qw/ Dumper /;
use DRMAA;

$|++;  # unbuffer STDOUT

use constant JOB_CHUNK => 2;
use constant NBULKS => 1;

sub create_job_template( $ $ $ );
sub identify_system();

my $job = "/home/harsch/tmp/DRMAA_TEST/SWIG/sl.csh";
my $rps_vals = {
$DRMAA::DRMAA_PS_UNDETERMINED         	=> 'DRMAA_PS_UNDETERMINED',
$DRMAA::DRMAA_PS_QUEUED_ACTIVE        	=> 'DRMAA_PS_QUEUED_ACTIVE',
$DRMAA::DRMAA_PS_SYSTEM_ON_HOLD       	=> 'DRMAA_PS_SYSTEM_ON_HOLD',
$DRMAA::DRMAA_PS_USER_ON_HOLD         	=> 'DRMAA_PS_USER_ON_HOLD',
$DRMAA::DRMAA_PS_USER_SYSTEM_ON_HOLD  	=> 'DRMAA_PS_USER_SYSTEM_ON_HOLD',
$DRMAA::DRMAA_PS_RUNNING              	=> 'DRMAA_PS_RUNNING',
$DRMAA::DRMAA_PS_SYSTEM_SUSPENDED     	=> 'DRMAA_PS_SYSTEM_SUSPENDED',
$DRMAA::DRMAA_PS_USER_SUSPENDED       	=> 'DRMAA_PS_USER_SUSPENDED',
$DRMAA::DRMAA_PS_USER_SYSTEM_SUSPENDED	=> 'DRMAA_PS_USER_SYSTEM_SUSPENDED',
$DRMAA::DRMAA_PS_DONE                 	=> 'DRMAA_PS_DONE',
$DRMAA::DRMAA_PS_FAILED               	=> 'DRMAA_PS_FAILED',
};


# let's identify ourselves, before initializing the DRM...
warn "*** DRM System identification ***";
identify_system();
	
my( $error, $diagnosis ) = &DRMAA::drmaa_init( undef );
die &DRMAA::drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

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
		while( ( $error, $jobids, $diagnosis ) = DRMAA::drmaa_run_bulk_jobs( $jt, 1, &JOB_CHUNK, 1 ) ) {
			if( $error == $DRMAA::DRMAA_ERRNO_DRM_COMMUNICATION_FAILURE ) {
				warn &DRMAA::drmaa_strerror( $error ) . "\n$diagnosis\nCouldn't contact qmaster... retrying: ";
				sleep 1;
				next;
			} # end if
			if( $error != $DRMAA::DRMAA_ERRNO_SUCCESS ) {
				die &DRMAA::drmaa_strerror( $error ) . "\n$diagnosis\nError: drmaa_run_bulk_jobs failed.  Aborting... ";
			} # end if
			last;
		} # end while

		print "submitted bulk job with jobids:\n";
		for (my $j=0; $j<JOB_CHUNK; $j++) {
			my( $error, $jobid ) = DRMAA::drmaa_get_next_job_id( $jobids );
			die "\$error='$error' \$jobid='$jobid'" if $error != $DRMAA::DRMAA_ERRNO_SUCCESS;
			push @all_jobids, $jobid;
			print "\t \"$jobid\"\n";
		}
		DRMAA::drmaa_release_job_ids($jobids);
	} # end for

	( $error, $diagnosis ) = DRMAA::drmaa_delete_job_template( $jt );
	die &DRMAA::drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

	#
	# Submit some sequential jobs
	#
	$jt = create_job_template( $job, 1, 0 );

	for (my $j=0; $j<JOB_CHUNK; $j++) {
		my $jobid;
		while( ( $error, $jobid, $diagnosis ) = DRMAA::drmaa_run_job( $jt ) ) {
			if( $error == $DRMAA::DRMAA_ERRNO_DRM_COMMUNICATION_FAILURE ) {
				warn "Couldn't contact qmaster... retrying:  $error, $diagnosis ";
				sleep 1;
				next;
			} # end if
			if( $error != $DRMAA::DRMAA_ERRNO_SUCCESS ) {
				die "drmaa_run_bulk_jobs failed:  $error, $diagnosis ";
			} # end if
			last;
		} # end while

		print "\t \"$jobid\"\n";
		push @all_jobids, $jobid;
	} # end while

	( $error, $diagnosis ) = DRMAA::drmaa_delete_job_template( $jt );
	die &DRMAA::drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

	#
	#   synchronize with all jobs
	#
	if( 0 ) {
		( $error, $diagnosis ) = DRMAA::drmaa_synchronize( \@all_jobids, $DRMAA::DRMAA_TIMEOUT_WAIT_FOREVER, 0 );
		die &DRMAA::drmaa_strerror( $error ) . "\n" . $diagnosis if $error;
		warn("synchronized with all jobs\n");
	} else {
		do {
			my @job_constant = ( $DRMAA::DRMAA_JOB_IDS_SESSION_ALL );
			# TODO: setting dispose to 1 causes seg fault
			( $error, $diagnosis ) = DRMAA::drmaa_synchronize( \@job_constant, 3, 1 );
			die &DRMAA::drmaa_strerror( $error ) . "\n" . $diagnosis
				if $error and $error != $DRMAA::DRMAA_ERRNO_EXIT_TIMEOUT;

			foreach my $jobid ( @all_jobids ) {
				my ( $error, $j, $remoteps, $diagnosis ) = &DRMAA::drmaa_job_ps( $jobid );
				die &DRMAA::drmaa_strerror( $error ) . "\n" . $diagnosis if $error;
				warn "\$jobid = $jobid, \$remoteps = $rps_vals->{$remoteps} ";
			}
		} while ( $error == $DRMAA::DRMAA_ERRNO_EXIT_TIMEOUT );
		warn("synchronized with all jobs\n");
	} # end if

}; # end eval

die $@ if $@;


sub create_job_template( $ $ $ ) {
	my $job_path = shift;
	my $seconds = shift;
	my $as_bulk_job = shift;

	my( $error, $jt, $diagnosis ) = &DRMAA::drmaa_allocate_job_template();
	die &DRMAA::drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

	warn "run in users home directory";
	( $diagnosis ) = &DRMAA::drmaa_set_attribute( $jt, $DRMAA::DRMAA_WD, $DRMAA::DRMAA_PLACEHOLDER_HD );
	die &DRMAA::drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

	warn "the job to be run";
	( $diagnosis ) = &DRMAA::drmaa_set_attribute( $jt, $DRMAA::DRMAA_REMOTE_COMMAND, $job );
	die &DRMAA::drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

	my $value;
	( $error, $value, $diagnosis ) = &DRMAA::drmaa_get_attribute( $jt, $DRMAA::DRMAA_REMOTE_COMMAND );
	die &DRMAA::drmaa_strerror( $error ) . "\n" . $diagnosis if $error;
	warn "... was submitted as <$value> ";

	warn "Set job parameters ";
	&DRMAA::drmaa_set_vector_attribute( $jt, $DRMAA::DRMAA_V_ARGV, [ "1" ] );  # 1st argv param = 1 (to run for 1 second)

	warn "join output/error file";
	( $diagnosis ) = &DRMAA::drmaa_set_attribute( $jt, $DRMAA::DRMAA_JOIN_FILES, "y" );
	die &DRMAA::drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

	warn "path for output";
	if( ! $as_bulk_job ) {
		( $diagnosis ) = DRMAA::drmaa_set_attribute( $jt, $DRMAA::DRMAA_OUTPUT_PATH,
										":" . $DRMAA::DRMAA_PLACEHOLDER_HD . "/DRMAA_JOB" );
	} else {
		( $diagnosis ) = DRMAA::drmaa_set_attribute( $jt, $DRMAA::DRMAA_OUTPUT_PATH,
											":" . $DRMAA::DRMAA_PLACEHOLDER_HD . "/DRMAA_JOB." . $DRMAA::DRMAA_PLACEHOLDER_INCR );
	}
	die &DRMAA::drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

	return $jt;
} # end sub


sub identify_system() {
	my( $error, $major, $minor, $diagnosis ) = &DRMAA::drmaa_version();
	die &DRMAA::drmaa_strerror( $error ) . "\n" . $diagnosis if $error;
	warn "Using DRMAA version \"$major.$minor\" ";
	
	my $DRMAA_impl;
	( $error, $DRMAA_impl, $diagnosis ) = &DRMAA::drmaa_get_DRMAA_implementation();
	die &DRMAA::drmaa_strerror( $error ) . "\n" . $diagnosis if $error;
	warn "\"$DRMAA_impl\" will be the DRMAA implementation used to distribute your jobs ";

	my $DRM_system;
	( $error, $DRM_system, $diagnosis ) = &DRMAA::drmaa_get_DRM_system();
	die &DRMAA::drmaa_strerror( $error ) . "\n" . $diagnosis if $error;
	warn "\"$DRM_system\" will be the system used to distribute your jobs ";

	my $contact;
	( $error, $contact, $diagnosis ) = &DRMAA::drmaa_get_contact();
	die &DRMAA::drmaa_strerror( $error ) . "\n" . $diagnosis if $error;
	warn "\"$contact\" will be notified by the DRM if needed ";

} # end sub

