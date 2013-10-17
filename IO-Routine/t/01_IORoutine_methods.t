#!perl -T

use 5.14.2;
use strict;
use warnings FATAL => 'all';
use Test::More tests => 30;

my ($pathTest, $errorTest, $openFileTest, $tmpLoc) = 0;
$|=0;

diag("Tests begin ...");
BEGIN {
    use_ok( 'IO::Routine' ) || print "Bail out! Cannot load IO::Routine for testing!\n";
}

diag( "Testing IO::Routine $IO::Routine::VERSION, Perl $], $^X" );

my $io = IO::Routine->new();
diag( "Testing constructor new()" );
ok( defined $io, 'IO::Routine->new() returns an object.');

diag( "Testing for blessed object" );
ok( $io->isa( 'IO::Routine' ), 'IO::Routine object verified.' );

SKIP: {
    skip(' ... Skipping methods: check_and_load_modules()... Incomplete code.', 1);
    ok( defined( $io->check_and_load_modules( ['File::Basename', 'Pod::Usage']) ), 'check_and_load_modules() loads core / optional modules.');
}

diag( "Testing method verify_options()" );
ok( defined( $io ) && can_ok($io, 'verify_options'), 'verify_options() - Requires Pod::Usage when using IO::Routine.');

SKIP: {
  skip(' ... Skipping methods: verify_input_files() ... Unable to test without providing UNIX paths.', 1)
    if (!$pathTest);
  ok( defined( $io->verify_input_files([''], [''])) );
}

SKIP: {
  skip(' ... Skipping methods: validate_create_path() ... Unable to test without providing UNIX paths.', 1)
    if (!$pathTest);
  ok( defined( $io->validate_create_path('', '', '')) );
}

SKIP: {
  skip(' ... Skipping methods: error() ... since error() prints fatal errors to STDERR.', 1)
    if (!$errorTest);
  ok( $io->error( "Testing fatal errors ... " ) , 'error() - Is loadable and since fatal errors are printed to STDERR.');
}

SKIP: {
  skip(' ... Skipping methods: check_sys_level_cmds() ... since fatal errors are printed to STDERR.', 1)
    if (!$errorTest);
  ok( $io->check_sys_level_cmds( ['grep'], ['2.10'] ) , 'check_sys_level_cmds() - Is loadable and since prints fatal errors to STDERR.');
}

SKIP: {
  skip(' ... Skipping methods: get_mem_usage() ... since it uses system level command grep.', 1)
    if (!$errorTest);
  ok( $io->get_mem_usage() ne '', 'get_mem_usage() - Prints current memory usage to STDOUT. Need to auto flush the command line buffer.');
}

$ENV{'PATH'} = '/bin/';

diag( "Testing method execute_system_command()" );
ok( defined( $io ) && can_ok($io, 'execute_system_command'), 'execute_system_command() - Is loadable and prints debug messages to STDOUT.');

diag( "Testing method execute_get_sys_cmd_output()" );
ok( defined( $io ) && can_ok($io, 'execute_get_sys_cmd_output'), 'execute_get_sys_cmd_output() - Is loadable and prints debug messages to STDOUT.');

diag( "Testing method open_file()" );
if ( $ENV{'TMP'} ) {
  diag( "Using user preferences and setting temporary location to $ENV{'TMP'}" );
  $tmpLoc = $ENV{'TMP'}
}
elsif ( $ENV{'TEMP'} ) {
  diag( "Using user preferences and setting temporary location to $ENV{'TEMP'}" );
  $tmpLoc = $ENV{'TEMP'}
}
elsif ( $ENV{'TMPDIR'} ) {
  diag( "Using user preferences and setting temporary location to $ENV{'TMPDIR'}" );
  $tmpLoc = $ENV{'TMPDIR'}
}
else {
  $tmpLoc = '/tmp';
}

$tmpLoc .= '/' if ( $tmpLoc !~ m/\/$/ );
my $tmpLocFile = $tmpLoc . 'IORoutineOpenFileTest.txt';

ok( $io->open_file('>', $tmpLocFile) ne '', 'open_file() - Returns non-empty filehandle for IO operations.');

diag( "Testing method is_dir_empty()" );
ok( defined( $io->is_dir_empty("$tmpLoc") ), 'is_dir_empty() Is loadable, Checks for empty directory.');

diag( "Testing method warning()" );
ok( defined( $io ) && can_ok($io, 'warning'), 'warning() - Is loadable, Prints warning messages to STDOUT.');

diag( "Testing method this_script_info()" );
ok( defined( $io ) && can_ok($io, 'this_script_info'), 'this_script_info() - Is loadable, Prints the script banner to STDOUT.');

diag( "Testing method start_timer()" );
ok( defined( $io ) && can_ok($io, 'start_timer'), 'start_timer() - Is loadable, Returns reference to started timer.');

diag( "Testing method end_timer()" );
ok( defined( $io ) && can_ok($io, 'end_timer'), 'end_timer() - Is loadable, Returns elapsed time in floating point seconds.');

diag( "Testing method c_time()" );
ok( defined( $io ) && can_ok($io, 'c_time'), 'c_time - Is loadable, Returns ctime from localtime().');

diag( "Testing method file_basename()" );
ok( $io->file_basename($tmpLocFile, 'suffix') =~ m/txt$/i , 'file_basename() - Returns suffix of a file.');

diag(  "Testing method exist_sys_cmd()" );
ok( defined( $io ) && can_ok($io, 'exist_sys_cmd'), 'exist_sys_cmd - Is loadable, Returns if true.');

diag("\n\nTesting Complete!\n\n");
