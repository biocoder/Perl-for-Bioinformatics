Bioinformatics tools written in Perl [![Bitdeli Badge](https://d2weczhvl823v0.cloudfront.net/biocoder/perl-for-bioinformatics/trend.png)](https://bitdeli.com/free "Bitdeli Badge")
====================================
Most of the scripts here were written while I was working on different projects, which I think will be useful to others and can be extended / modified per there needs.

IO::Routine [![Build Status](https://travis-ci.org/biocoder/Perl-for-Bioinformatics.png?branch=master)](https://travis-ci.org/biocoder/Perl-for-Bioinformatics)
----------------------------
* The scripts use custom [IO::Routine](https://github.com/biocoder/Perl-for-Bioinformatics/tree/master/IO-Routine) Perl Module.

* Please see the installation instructions by browsing the [IO::Routine](https://github.com/biocoder/Perl-for-Bioinformatics/tree/master/IO-Routine) directory.

* If you are installing **ncRNAScan** Pipeline, `IO::Routine` module is automatically installed.

ncRNAScan [![Build Status](https://travis-ci.org/biocoder/Perl-for-Bioinformatics.png?branch=master)](https://travis-ci.org/biocoder/Perl-for-Bioinformatics)
--------------------------
* A pipeline to extract putative novel ncRNAs ab initio, given a list of transcripts in GTF format assembled from deep sequencing data (ex: RNA-Seq) and annotation data.

* Head on to [NGS-Utils](https://github.com/biocoder/Perl-for-Bioinformatics/tree/master/NGS-Utils) directory for script list.

* Install ncRNAScan (Mac and Linux):

          cd /to/your/preferred/install/path
          curl -O https://raw.github.com/biocoder/Perl-for-Bioinformatics/master/NGS-Utils/ncRNAScan
          perl ncRNAScan -setup

* Documentation:
          
          perl ncRNAScan -h
or

          perldoc ncRNAScan
or to get help documentation for individual modules, do:

      	  perl ncRNAScan -h cuff
      	  perl ncRNAScan -h cat
      	  perl ncRNAScan -h get
      	  perl ncRNAScan -h fetch
      	  perl ncRNAScan -h cpc

=========
Cheers,

BioCoder
