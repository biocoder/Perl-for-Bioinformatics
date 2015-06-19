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
#  $Id: 03_system_queries.t,v 1.2 2004/04/27 20:50:36 harsch Exp $
##############################################################################

use Schedule::DRMAAc qw/ :all /;

use strict;

use constant TESTS => 15;
use Test::More tests => TESTS;

sub mydiag( $ $ ; $ );

my( $error, $diagnosis );
my @vals = drmaa_strerror( 0 );
ok( @vals == 1, "drmaa_strerror returned " . scalar(@vals) . " of 1 args" );

@vals = ( $error, $diagnosis ) = drmaa_init( undef );
ok( @vals == 2, "drmaa_init returned " . scalar(@vals) . " of 2 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_init error?" )
	or mydiag( $error, $diagnosis );

my $jt;
@vals = ( $error, $jt, $diagnosis ) = drmaa_allocate_job_template();
ok( @vals == 3, "drmaa_init returned " . scalar(@vals) . " of 3 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_allocate_job_template error?" )
	or mydiag( $error, $diagnosis );

my( $major, $minor );
@vals = ( $error, $major, $minor, $diagnosis ) = drmaa_version();
ok( @vals == 4, "drmaa_version returned " . scalar(@vals) . " of 4 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_version error?" )
	or mydiag( $error, $diagnosis );

my $DRMAA_impl;
@vals = ( $error, $DRMAA_impl, $diagnosis ) = drmaa_get_DRMAA_implementation();
ok( @vals == 3, "drmaa_version returned " . scalar(@vals) . " of 3 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_version error?" )
	or mydiag( $error, $diagnosis );

my $DRM_system;
@vals = ( $error, $DRM_system, $diagnosis ) = drmaa_get_DRM_system();
ok( @vals == 3, "drmaa_get_DRM_system returned " . scalar(@vals) . " of 3 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_get_DRM_system error?" )
	or mydiag( $error, $diagnosis );

my $contact;
@vals = ( $error, $contact, $diagnosis ) = drmaa_get_contact();
ok( @vals == 3, "drmaa_get_contact returned " . scalar(@vals) . " of 3 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_get_contact error?" )
	or mydiag( $error, $diagnosis );

@vals = ( $error, $diagnosis ) = drmaa_exit();
ok( @vals == 2, "drmaa_exit returned " . scalar(@vals) . " of 2 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_exit error?" )
	or mydiag( $error, $diagnosis );


sub mydiag( $ $ ; $ ) {
	my $msg = "error: " . drmaa_strerror( $_[0] ) . "\n" . "diagnosis: " . $_[1];
	$msg .= $_[2] if defined $_[2];
	diag( $msg );
} # end sub

1; # Ancient Druid Custom
