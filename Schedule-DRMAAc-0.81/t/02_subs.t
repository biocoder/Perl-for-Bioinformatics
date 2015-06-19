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
#  $Id: 02_subs.t,v 1.2 2004/04/27 20:50:36 harsch Exp $
##############################################################################

use Schedule::DRMAAc qw/ :all /;

use Test::More;

my @funcs =  qw/
	drmaa_get_next_attr_name drmaa_get_next_attr_value drmaa_get_next_job_id drmaa_release_attr_names
	drmaa_release_attr_values drmaa_release_job_ids drmaa_init drmaa_exit drmaa_allocate_job_template
	drmaa_delete_job_template drmaa_set_attribute drmaa_get_attribute drmaa_set_vector_attribute
	drmaa_get_vector_attribute drmaa_get_attribute_names drmaa_get_vector_attribute_names drmaa_run_job
	drmaa_run_bulk_jobs drmaa_control drmaa_synchronize drmaa_wait drmaa_wifexited drmaa_wexitstatus
	drmaa_wifsignaled drmaa_wtermsig drmaa_wcoredump drmaa_wifaborted drmaa_job_ps drmaa_strerror
	drmaa_get_contact drmaa_version drmaa_get_DRM_system drmaa_get_DRMAA_implementation
	/;

plan tests => scalar( @funcs );

foreach my $x ( @funcs ) {
	ok( defined *{$x}, "$x defined? " );
}

1; # Ancient Druid Custom
