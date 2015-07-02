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
 *  $Id: abort.c,v 1.2 2004/04/27 20:50:34 harsch Exp $
 *****************************************************************************/

#include <stdio.h>

int *t;
main( ) {
    t = NULL;
    // should cause this program to segfault
    *t = 123;
}
