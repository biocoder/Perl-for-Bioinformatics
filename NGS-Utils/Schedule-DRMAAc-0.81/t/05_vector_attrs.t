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
#  $Id: 05_vector_attrs.t,v 1.2 2004/04/27 20:50:36 harsch Exp $
##############################################################################

use Schedule::DRMAAc qw/ :all /;

use strict;

use constant TESTS => 15;
use Test::More tests => TESTS;

sub mydiag( $ $ ; $ );

my $vec_args = [ 1 ];

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

# check correct vector attributes available from DRM
my $vec_names;
@vals = ( $error, $vec_names, $diagnosis) = drmaa_get_vector_attribute_names();
ok( @vals == 3, "drmaa_get_attribute_names returned " . scalar(@vals) . " of 3 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_get_attribute_names error?" )
	or mydiag( $error, $diagnosis );

my @vec_names;
while( my( $error, $vec_name ) = drmaa_get_next_attr_name( $vec_names ) ) {
	last if $error != $DRMAA_ERRNO_SUCCESS;

	push @vec_names, $vec_name; # record all DRM attrs
} # end while

my %legit_vec_attrs = ( # 1 mandatory, 2 optional
			    $DRMAA_V_ARGV => 1,
			    $DRMAA_V_ENV => 1,
                   $DRMAA_V_EMAIL => 1,
                   );
my $all_man_defined = 0;
my $all_opt_defined = 0;
foreach my $vec_name ( @vec_names ) {
	if( defined $legit_vec_attrs{$vec_name} ) {
		$all_man_defined++ if $legit_vec_attrs{$vec_name} == 1;
		$all_opt_defined++ if $legit_vec_attrs{$vec_name} == 2;
	} # end if
} # end foreach
is( $all_man_defined, scalar( grep{ $legit_vec_attrs{$_} == 1} keys %legit_vec_attrs), "DRM returned enough mandatory vector attributes?" );
ok( $all_opt_defined >= 0, "DRM returned $all_opt_defined of " . scalar( grep{ $legit_vec_attrs{$_} == 2} keys %legit_vec_attrs) . " possible optional attributes" );

@vals = drmaa_release_attr_names( $vec_names );

# test setting one vector attribute and retrieving its values

@vals = ( $error, $diagnosis ) = drmaa_set_vector_attribute( $jt, $DRMAA_V_ARGV, $vec_args );
ok( @vals == 2, "drmaa_set_vector_attribute returned " . scalar(@vals) . " of 2 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_set_vector_attribute error?" )
	or mydiag( $error, $diagnosis );

my $vec_vals;
@vals = ( $error, $vec_vals, $diagnosis ) = drmaa_get_vector_attribute( $jt, $DRMAA_V_ARGV );
ok( @vals == 3, "drmaa_get_next_attr_name returned " . scalar(@vals) . " of 3 args" );
ok( $error == $DRMAA_ERRNO_SUCCESS, "drmaa_get_next_attr_name error?" )
	or mydiag( $error, $diagnosis );

my $vec_val;
@vals = ( $error, $vec_val ) = drmaa_get_next_attr_value( $vec_vals );
ok( @vals == 2, "drmaa_get_next_attr_value returned " . scalar(@vals) . " of 2 args" );

my @vec_vals = ( $vec_val );
while( my( $error, $vec_val ) = drmaa_get_next_attr_value( $vec_vals ) ) {
	last if $error != $DRMAA_ERRNO_SUCCESS;

	# should never reach this line if this module is properly working.  We'll check in a bit...
	push @vec_vals, $vec_val; # record all values of the vector
} # end while

is( @vec_vals, 1, "Proper number of values in the vector attribute we set?" );
is( $vec_vals[0], $vec_args->[0], "Correct value in the vector we set?" );

@vals = drmaa_release_attr_values( $vec_vals );
ok( @vals == 0, "drmaa_release_attr_values returned " . scalar(@vals) . " of 0 args" );


sub mydiag( $ $ ; $ ) {
	my $msg = "error: " . drmaa_strerror( $_[0] ) . "\n" . "diagnosis: " . $_[1];
	$msg .= $_[2] if defined $_[2];
	diag( $msg );
} # end sub

1; # Ancient Druid Custom
