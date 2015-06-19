/*****************************************************************************
 *
 *  Copyright (c) 2004, The Regents of the University of California.
 *  Produced at the Lawrence Livermore National Laboratory.
 *  Written by Tim Harsch <harsch1@llnl.gov>
 *  UCRL-CODE-155918
 *  All rights reserved.
 *
 *  This file is part of Schedule::DRMAAc. For details, see CPAN
 *  Please also read LICENSE.txt which is found in this source distribution.
 *
 *  This program is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License (as published by the
 *  Free Software Foundation) version 2, dated June 1991.
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY
 *  OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the terms and conditions of the GNU General Public License for more
 *  details.
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software Foundation,
 *  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 *****************************************************************************
 *  $Id: Schedule_DRMAAc.i,v 1.2 2004/04/27 20:50:30 harsch Exp $
 *****************************************************************************/

%module "Schedule::DRMAAc"

%include "typemaps.i"

%{
#include "drmaa.h"

%}

%typemap(in,numinputs=0) (char * error_diagnosis, size_t error_diag_len ) {
	// %typemap(in,numinputs=0) ($1_type $1_name, $2_type $2_name )
	$1 = malloc( sizeof(char) * DRMAA_ERROR_STRING_BUFFER );
	$2 = DRMAA_ERROR_STRING_BUFFER;
}
%apply (char * error_diagnosis, size_t error_diag_len ) { ( char *value, size_t value_len ) }
%typemap(in,numinputs=0) ( char *job_id, size_t job_id_len ) {
	// %typemap(in,numinputs=0) ($1_type $1_name, $2_type $2_name )
	$1 = malloc( sizeof(char) * DRMAA_JOBNAME_BUFFER );
	$2 = DRMAA_JOBNAME_BUFFER;
}
%apply (char *job_id, size_t job_id_len ) { ( char *job_id_out, size_t job_id_out_len ) }
%typemap(in,numinputs=0) ( char *drm_system, size_t drm_system_len ) {
	// %typemap(in,numinputs=0) ($1_type $1_name, $2_type $2_name )
	$1 = malloc( sizeof(char) * DRMAA_DRM_SYSTEM_BUFFER );
	$2 = DRMAA_DRM_SYSTEM_BUFFER;
}
%apply (char *drm_system, size_t drm_system_len ) { ( char *drmaa_impl, size_t drmaa_impl_len ) }
%typemap(in,numinputs=0) ( char *contact, size_t contact_len ) {
	// %typemap(in,numinputs=0) ($1_type $1_name, $2_type $2_name )
	$1 = malloc( sizeof(char) * DRMAA_CONTACT_BUFFER );
	$2 = DRMAA_CONTACT_BUFFER;
}
%typemap(in,numinputs=0) ( char *signal, size_t signal_len ) {
	// %typemap(in,numinputs=0) ($1_type $1_name, $2_type $2_name )
	$1 = malloc( sizeof(char) * DRMAA_SIGNAL_BUFFER );
	$2 = DRMAA_SIGNAL_BUFFER;
}


%typemap(argout) char *error_diagnosis {
	// %typemap(argout) $1_type $1_name
	$result = sv_newmortal();
	sv_setpv($result, $1);
	argvi++;                     /* Increment return count -- important! */
}
%apply char *error_diagnosis { char *value }
%apply char *error_diagnosis { char *job_id }
%apply char *error_diagnosis { char *job_id_out }
%apply char *error_diagnosis { char *contact }
%apply char *error_diagnosis { char *drm_system }
%apply char *error_diagnosis { char *drmaa_impl };

%typemap(argout) const char * "pass";

%typemap(in,numinputs=0) drmaa_job_template_t ** ( $*1_type keepme )  {
	// %typemap(in) $1_type $1_name
	keepme = NULL;

	$1 = &keepme;
}

%typemap(argout) drmaa_job_template_t ** {
	// %typemap(argout) $1_type $1_name
	ST(argvi) = sv_newmortal();
	SWIG_MakePtr(ST(argvi++), (void *)*$1, $*1_descriptor, 0);
}

%apply drmaa_job_template_t ** { drmaa_job_ids_t ** }
%apply drmaa_job_template_t ** { drmaa_attr_names_t ** }
%apply drmaa_job_template_t ** { drmaa_attr_values_t ** }


%typemap(in) const char *[] {
	// %typemap(argout) $1_type $1_name
	AV *tempav;
	I32 len;
	int i;
	SV  **tv;
	char *myConst;

	if (!SvROK($input)) {
	    croak("Error: $symname: Argument $argnum is not a reference ");
	}
        if (SvTYPE(SvRV($input)) != SVt_PVAV) {
	    croak("Error: $symname: Argument $argnum is not an array ");
	}
        tempav = (AV*)SvRV($input);
	len = av_len(tempav);
	$1 = (char **) malloc((len+2)*sizeof(char *));
	for (i = 0; i <= len; i++) {
	    tv = av_fetch(tempav, i, 0);	
	    $1[i] = (char *) SvPV(*tv,PL_na);
        }
	$1[i] = NULL;
};

%apply unsigned int *OUTPUT { unsigned int * }
%apply int *OUTPUT { int * }

%pragma(perl5) include="Schedule_DRMAAc_ext.pm"
%pragma(perl5) include="Schedule_DRMAAc.pod"

%include "drmaa.h"

