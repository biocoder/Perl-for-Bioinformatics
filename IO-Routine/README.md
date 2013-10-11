IO::Routine [![Build Status](https://travis-ci.org/biocoder/Perl-for-Bioinformatics.png?branch=master)](https://travis-ci.org/biocoder/Perl-for-Bioinformatics)
===========

An attempt to provide a solution to avoid routine IO chores.

Version
-------
0.28

Installation
------------

To install this module, run the following commands:

	perl Makefile.PL
	make
	make test
	make install

To install IO::Routine at custom location, run the following commands:

	perl Makefile.PL PREFIX=/path/to/custom/perllib LIB=/path/to/custom/perllib/lib
	make
	make test
	make install

and then in your script add that path to your @INC variable as:

`BEGIN{ push (@INC, '/path/to/custom/perllib/lib') };`

Or, for UNIX systems, edit your `.bash_profile` by adding custom location to PERL5LIB

`export PERL5LIB=$PERL5LIB:/path/to/custom/perllib/lib`

Support and Documentation
-------------------------

After installing, you can find documentation for this module with the
perldoc command.

    perldoc IO::Routine

You can also look for information at:

1. RT, CPAN's request tracker (report bugs here)  
  * http://rt.cpan.org/NoAuth/Bugs.html?Dist=IO-Routine  

2. AnnoCPAN, Annotated CPAN documentation  
  * http://annocpan.org/dist/IO-Routine

3. CPAN Ratings  
  * http://cpanratings.perl.org/d/IO-Routine

4. Search CPAN  
  * http://search.cpan.org/dist/IO-Routine/

