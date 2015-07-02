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
#  $Id: simple.pl,v 1.2 2004/04/27 20:50:35 harsch Exp $
##############################################################################

use Schedule::DRMAAc qw/ :all /;

( $error, $diagnosis ) = drmaa_init( undef );
die drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

( $error, $jt, $diagnosis ) = drmaa_allocate_job_template();
die drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

( $error, $diagnosis ) = drmaa_set_attribute( $jt, $DRMAA_REMOTE_COMMAND, '/bin/sleep' );
die drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

( $error, $diagnosis ) = drmaa_set_vector_attribute( $jt, $DRMAA_V_ARGV, [ "1" ] );  # 1st argv param = 1 (to run for 1 second)
die drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

( $error, $jobid, $diagnosis ) = drmaa_run_job( $jt );
die drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

my @job_constant = ( $DRMAA_JOB_IDS_SESSION_ALL );
( $error, $diagnosis ) = drmaa_synchronize( \@job_constant , $DRMAA_TIMEOUT_WAIT_FOREVER, 0 );
die drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

( $error, $diagnosis ) = drmaa_delete_job_template( $jt );
die drmaa_strerror( $error ) . "\n" . $diagnosis if $error;

