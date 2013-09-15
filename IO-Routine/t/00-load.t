#!perl -T
use 5.14.2;
use strict;
use warnings FATAL => 'all';
use Test::More;

plan tests => 1;

BEGIN {
    use_ok( 'IO::Routine' ) || print "Bail out!\n";
}

diag( "Testing IO::Routine $IO::Routine::VERSION, Perl $], $^X" );
