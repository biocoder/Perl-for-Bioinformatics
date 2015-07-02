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
#  $Id: Makefile.PL,v 1.4 2004/04/27 20:50:30 harsch Exp $
##############################################################################

use ExtUtils::MakeMaker;
use Config;

die <<EOF unless my $SGE_ROOT = $ENV{SGE_ROOT};
SGE_ROOT environment variable not defined.  Currently linking
support to other DRM's has not been added to Schedule::DRMAAc.
If you are attempting to use this module with another DRM.  Try
manually editing the Makefile.PL file.  ( And inform harsch1\@llnl.gov
what you needed to do :-)
EOF
die "You forgot to create a link to drmaa.h, see README" unless -e 'drmaa.h';

my $arch = `$SGE_ROOT/util/arch`;
chomp $arch;
my $sgelib = "$SGE_ROOT/lib/$arch";
if( $Config{longsize} == 4 ) {
    if( $arch eq 'sol-sparc64' ) {
        if( -e "$SGE_ROOT/lib/sol-sparc" ) {
		$sge_lib = "$SGE_ROOT/lib/sol-sparc";
        } else {
		warn "You are trying to compile 64bit SGE libs onto a 32 bit perl.  Not likely to work.  Try compiling 32bit SGE from sources";
        } # end if
    } # end if
} # end if

if( $Config{longsize} == 8 ) {
    if( $arch eq 'sol-sparc' ) {
        if( -e "$SGE_ROOT/lib/sol-sparc64" ) {
		$sge_lib = "$SGE_ROOT/lib/sol-sparc64";
        } else {
		warn "You are trying to compile 32bit SGE libs onto a 64 bit perl.  Not likely to work.  Try compiling 64bit SGE from sources";
        } # end if
    } # end if
} # end if



# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'		=> 'Schedule::DRMAAc',
    'VERSION_FROM'	=> 'DRMAAc.pm', # finds $VERSION
    'PREREQ_PM'		=> {}, # e.g., Module::Name => 1.1
    ($] >= 5.005 ?    ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM => 'DRMAAc.pm', # retrieve abstract from module
       AUTHOR     => 'Tim Harsch <harsch1@llnl.gov>') : ()),
    'LIBS'		=> ["-L$sgelib -ldrmaa -lsocket -lnsl -lm -lpthread"],
    'DEFINE'		=> '', # e.g., '-DHAVE_SOMETHING'
	# Insert -I. if you add *.h files later:
	'INC'		=> "-I$SGE_ROOT/include",
	
	# Un-comment this if you add C files to link with later:
	'OBJECT'	=> '$(O_FILES)', # link all the C files too
);
