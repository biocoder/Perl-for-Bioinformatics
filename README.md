Bioinformatics tools written in Perl
====================================
Most of the scripts here were written while I was working on different projects, which I think will be useful to others and can be extended / modified per there needs.

IO::Routine [![Build Status](https://api.travis-ci.org/biocoder/Perl-for-Bioinformatics.png?branch=master)](https://travis-ci.org/biocoder/Perl-for-Bioinformatics)
----------------------------
* The scripts use custom [IO::Routine](https://github.com/biocoder/Perl-for-Bioinformatics/tree/master/IO-Routine) Perl Module.

* Please see the installation instructions by browsing the [IO::Routine](https://github.com/biocoder/Perl-for-Bioinformatics/tree/master/IO-Routine) directory.

* If you are installing **ncRNAScan** Pipeline, `IO::Routine` module is automatically installed.

☲☴ ncRNAScan [![Build Status](https://api.travis-ci.org/biocoder/Perl-for-Bioinformatics.png?branch=master)](https://travis-ci.org/biocoder/Perl-for-Bioinformatics)
--------------------------
* A pipeline to extract putative novel ncRNAs ab initio, given a list of transcripts in GTF format assembled from deep sequencing data (ex: RNA-Seq) and annotation data.

* Head on to [NGS-Utils](https://github.com/biocoder/Perl-for-Bioinformatics/tree/master/NGS-Utils) directory for script list.

* Install ncRNAScan and all its dependencies (Mac and Linux):

          cd /to/your/preferred/install/path
          curl -O https://raw.githubusercontent.com/biocoder/Perl-for-Bioinformatics/master/NGS-Utils/ncRNAScan
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
      	  perl ncRNAScan -h rna
      	  perl ncRNAScan -h inf
      	  
* Known issues:

     * If pipeline setup fails due to `XML::Parser` module, you need to install XML parser C libraries.
     * On Ubuntu / Debian based Linux distributions, as `root` user, do:
     
            apt-get install libexpat1 libexpat1-dev
                    
     * On RedHat / Fedora / CentOS based Linux distributions, as `root` user do:
     
            yum install expat expat-devel
            
* Caveats:
    * The pipeline script uses a lot of inherent Linux core utils and has been only tested in BASH shell. 
 
           


=========
Cheers,

BioCoder
