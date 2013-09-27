use 5.14.2;
use strict;
use warnings FATAL => 'all';
use ExtUtils::MakeMaker;

WriteMakefile(
	      NAME             => 'IO::Routine',
	      AUTHOR           => q{Kranti Konganti <konganti@tamu.edu>},
	      VERSION_FROM     => 'lib/IO/Routine.pm',
	      ABSTRACT_FROM    => 'lib/IO/Routine.pm',
	      LICENSE          => 'Artistic_2_0',
	      PL_FILES         => {},
	      MIN_PERL_VERSION => 5.14.2,
	      CONFIGURE_REQUIRES => {
				     'ExtUtils::MakeMaker' => 0,
				    },
	      BUILD_REQUIRES => {
				 'Test::More' => 0,
				},
	      PREREQ_PM => {
		            'strict'         => 0,
		            'warnings'       => 0,
		            'Carp'           => 0,
			    'Module::Load'   => 0,
			    'Pod::Usage'     => 0,
			    'Exporter'       => 0,
			    'Time::HiRes'    => 0,
			    'Cwd'            => 0
			   },
	      dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
	      clean => { FILES => 'IO-Routine-*' },
	     );
