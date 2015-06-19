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
#  $Id: 04_standard_attrs.t,v 1.2 2004/04/27 20:50:36 harsch Exp $
##############################################################################

use Schedule::DRMAAc qw/ :all /;

use strict;

use constant TESTS => 17;
use Test::More tests => TESTS;

sub mydiag( $ $ ; $ );

my $remote_cmd = "/bin/sleep";


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

@vals = ( $error, $diagnosis ) = drmaa_set_attribute( $jt, $DRMAA_REMOTE_COMMAND, $remote_cmd );
ok( @vals == 2, "drmaa_set_attribute returned " . scalar(@vals) . " of 2 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_set_attribute error?" )
	or mydiag( $error, $diagnosis );

my $value;
@vals = ( $error, $value, $diagnosis ) = drmaa_get_attribute( $jt, $DRMAA_REMOTE_COMMAND );
ok( @vals == 3, "drmaa_get_attribute returned " . scalar(@vals) . " of 3 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_get_attribute error?" )
	or mydiag( $error, $diagnosis );
is( $value, $remote_cmd, "DRM return the same value as we set?" );

my $names;
@vals = ( $error, $names, $diagnosis) = drmaa_get_attribute_names();
ok( @vals == 3, "drmaa_get_attribute_names returned " . scalar(@vals) . " of 3 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_get_attribute_names error?" )
	or mydiag( $error, $diagnosis );

my( $retval, $name );
@vals = ( $retval, $name ) = drmaa_get_next_attr_name( $names );
ok( @vals == 2, "drmaa_get_next_attr_name returned " . scalar(@vals) . " of 2 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_get_next_attr_name error?" )
	or mydiag( $error, $diagnosis );

my @names = ( $name );
while( my( $retval, $name ) = drmaa_get_next_attr_name( $names ) ) {
	last if $retval != $DRMAA_ERRNO_SUCCESS;

	push @names, $name; # record all DRM attrs
} # end while

my %legit_std_attrs = ( # 1 mandatory, 2 optional
			    $DRMAA_REMOTE_COMMAND => 1,
			    $DRMAA_JS_STATE => 1,
                   $DRMAA_WD => 1,
			    $DRMAA_JOB_CATEGORY => 1,
                   $DRMAA_NATIVE_SPECIFICATION => 1,
                   $DRMAA_BLOCK_EMAIL => 1,
                   $DRMAA_START_TIME => 1,
                   $DRMAA_JOB_NAME => 1,
                   $DRMAA_INPUT_PATH => 1,
                   $DRMAA_OUTPUT_PATH => 1,
                   $DRMAA_ERROR_PATH => 1,
                   $DRMAA_JOIN_FILES => 1,
                   $DRMAA_TRANSFER_FILES => 2,
                   $DRMAA_DEADLINE_TIME => 2,
                   $DRMAA_WCT_HLIMIT => 2,
                   $DRMAA_WCT_SLIMIT => 2,
                   $DRMAA_DURATION_HLIMIT => 2,
                   $DRMAA_DURATION_SLIMIT => 2,
                   );
my $all_std_defined = 0;
my $all_opt_defined = 0;
foreach my $attr ( @names ) {
	if( defined $legit_std_attrs{$name} ) {
		$all_std_defined++ if $legit_std_attrs{$name} == 1;
		$all_opt_defined++ if $legit_std_attrs{$name} == 2;
	} # end if\
} # end foreach
is( $all_std_defined, scalar( grep{ $legit_std_attrs{$_} == 1} keys %legit_std_attrs), "DRM returned enough mandatory standard attributes?" );
ok( $all_opt_defined >= 0, "DRM returned $all_opt_defined of " . scalar( grep{ $legit_std_attrs{$_} == 2} keys %legit_std_attrs) . " possible optional attributes" );

@vals = drmaa_release_attr_names( $names );
ok( @vals == 0, "drmaa_release_attr_names returned " . scalar(@vals) . " of 0 args" );

sub mydiag( $ $ ; $ ) {
	my $msg = "error: " . drmaa_strerror( $_[0] ) . "\n" . "diagnosis: " . $_[1];
	$msg .= $_[2] if defined $_[2];
	diag( $msg );
} # end sub

1; # Ancient Druid Custom
