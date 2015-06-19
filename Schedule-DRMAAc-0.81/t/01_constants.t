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
#  $Id: 01_constants.t,v 1.3 2004/04/27 20:50:36 harsch Exp $
##############################################################################

use Schedule::DRMAAc qw/ :all /;

use Test::More;

my @constants = qw/
	DRMAA_TIMEOUT_WAIT_FOREVER DRMAA_TIMEOUT_NO_WAIT DRMAA_JOB_IDS_SESSION_ANY DRMAA_JOB_IDS_SESSION_ALL
	DRMAA_SUBMISSION_STATE_ACTIVE DRMAA_SUBMISSION_STATE_HOLD DRMAA_PLACEHOLDER_INCR DRMAA_PLACEHOLDER_HD
	DRMAA_PLACEHOLDER_WD DRMAA_REMOTE_COMMAND DRMAA_JS_STATE DRMAA_WD DRMAA_JOB_CATEGORY DRMAA_NATIVE_SPECIFICATION
	DRMAA_BLOCK_EMAIL DRMAA_START_TIME DRMAA_JOB_NAME DRMAA_INPUT_PATH DRMAA_OUTPUT_PATH DRMAA_ERROR_PATH
	DRMAA_JOIN_FILES DRMAA_TRANSFER_FILES DRMAA_V_ARGV DRMAA_V_ENV DRMAA_V_EMAIL DRMAA_ERRNO_SUCCESS
	DRMAA_ERRNO_INTERNAL_ERROR DRMAA_ERRNO_DRM_COMMUNICATION_FAILURE DRMAA_ERRNO_AUTH_FAILURE
	DRMAA_ERRNO_INVALID_ARGUMENT DRMAA_ERRNO_NO_ACTIVE_SESSION DRMAA_ERRNO_NO_MEMORY DRMAA_ERRNO_INVALID_CONTACT_STRING
	DRMAA_ERRNO_DEFAULT_CONTACT_STRING_ERROR DRMAA_ERRNO_DRMS_INIT_FAILED DRMAA_ERRNO_ALREADY_ACTIVE_SESSION
	DRMAA_ERRNO_DRMS_EXIT_ERROR DRMAA_ERRNO_INVALID_ATTRIBUTE_FORMAT DRMAA_ERRNO_INVALID_ATTRIBUTE_VALUE
	DRMAA_ERRNO_CONFLICTING_ATTRIBUTE_VALUES DRMAA_ERRNO_TRY_LATER DRMAA_ERRNO_DENIED_BY_DRM
	DRMAA_ERRNO_INVALID_JOB DRMAA_ERRNO_RESUME_INCONSISTENT_STATE DRMAA_ERRNO_SUSPEND_INCONSISTENT_STATE
	DRMAA_ERRNO_HOLD_INCONSISTENT_STATE DRMAA_ERRNO_RELEASE_INCONSISTENT_STATE DRMAA_ERRNO_EXIT_TIMEOUT
	DRMAA_ERRNO_NO_RUSAGE DRMAA_NO_ERRNO DRMAA_PS_UNDETERMINED DRMAA_PS_QUEUED_ACTIVE DRMAA_PS_SYSTEM_ON_HOLD
	DRMAA_PS_USER_ON_HOLD DRMAA_PS_USER_SYSTEM_ON_HOLD DRMAA_PS_RUNNING DRMAA_PS_SYSTEM_SUSPENDED
	DRMAA_PS_USER_SUSPENDED DRMAA_PS_USER_SYSTEM_SUSPENDED DRMAA_PS_DONE DRMAA_PS_FAILED DRMAA_CONTROL_SUSPEND
	DRMAA_CONTROL_RESUME DRMAA_CONTROL_HOLD DRMAA_CONTROL_RELEASE DRMAA_CONTROL_TERMINATE
	/;	
plan tests => scalar( @constants );

foreach my $x ( @constants ) {
	ok( defined *{$x}, "$x defined? " );
}

1; # Ancient Druid Custom
