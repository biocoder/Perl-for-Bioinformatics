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
#  $Id: 06_simple_job.t,v 1.1 2004/04/27 20:50:36 harsch Exp $
##############################################################################

use Schedule::DRMAAc qw/ :all /;

use strict;

use constant TESTS => 12;
use Test::More tests => TESTS;

sub mydiag( $ $ ; $ );

my $vec_args = [ 1 ];
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
ok( @vals == 2, "drmaa_set_attribute returned " . scalar(@vals) . " of 2 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_set_attribute error?" )
	or mydiag( $error, $diagnosis );

# test setting one vector attribute and retrieving its values
@vals = ( $error, $diagnosis ) = drmaa_set_vector_attribute( $jt, $DRMAA_V_ARGV, $vec_args );
ok( @vals == 2, "drmaa_set_vector_attribute returned " . scalar(@vals) . " of 2 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_set_vector_attribute error?" )
	or mydiag( $error, $diagnosis );

my $jobid;
@vals = ( $error, $jobid, $diagnosis ) = drmaa_run_job( $jt );
ok( @vals == 3, "drmaa_run_job returned " . scalar(@vals) . " of 3 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_run_job error?" )
	or mydiag( $error, $diagnosis );

my @job_constant = ( $DRMAA_JOB_IDS_SESSION_ALL );
@vals = ( $error, $diagnosis ) = drmaa_synchronize( \@job_constant , $DRMAA_TIMEOUT_WAIT_FOREVER, 0 );
ok( @vals == 2, "drmaa_synchronize returned " . scalar(@vals) . " of 2 args" );

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
