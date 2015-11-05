Bioinformatics tools written in Perl 
====================================

[![Join the chat at https://gitter.im/biocoder/Perl-for-Bioinformatics](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/biocoder/Perl-for-Bioinformatics?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

Most of the scripts here were written while I was working on different projects, which I think will be useful to others and can be extended / modified per there needs.

IO::Routine [![Build Status](https://api.travis-ci.org/biocoder/Perl-for-Bioinformatics.png?branch=master)](https://travis-ci.org/biocoder/Perl-for-Bioinformatics)
----------------------------
* The scripts use custom [IO::Routine](https://github.com/biocoder/Perl-for-Bioinformatics/tree/master/IO-Routine) Perl Module.

* Please see the installation instructions by browsing the [IO::Routine](https://github.com/biocoder/Perl-for-Bioinformatics/tree/master/IO-Routine) directory.

* If you are installing **lncRNApipe** Pipeline, `IO::Routine` module is automatically installed.

* Requires `Bio::SeqIO` module be installed and available.

☲☴ lncRNApipe [![Build Status](https://api.travis-ci.org/biocoder/Perl-for-Bioinformatics.png?branch=master)](https://travis-ci.org/biocoder/Perl-for-Bioinformatics)
--------------------------
* A pipeline to extract putative novel lncRNAs ab initio, given a list of transcripts in GTF format assembled from deep sequencing data (ex: RNA-Seq) and annotation data.

* Head on to [NGS-Utils](https://github.com/biocoder/Perl-for-Bioinformatics/tree/master/NGS-Utils) directory for script list.

* Install lncRNApipe and all its dependencies (Mac and Linux):

          cd /to/your/preferred/install/path
          curl -kO https://raw.githubusercontent.com/biocoder/Perl-for-Bioinformatics/master/NGS-Utils/lncRNApipe
          perl lncRNApipe --setup

* Documentation:
          
          perl lncRNApipe --h
or

          perldoc lncRNApipe
or to get help documentation for individual modules, do:

      	  perl lncRNApipe --h cuff
      	  perl lncRNApipe --h cat
      	  perl lncRNApipe --h get
      	  perl lncRNApipe --h fetch
      	  perl lncRNApipe --h cpc
      	  perl lncRNApipe --h rna
      	  perl lncRNApipe --h inf
      	  
* Known issues:

     * If pipeline setup fails due to `XML::Parser` module, you need to install XML parser C libraries.
     * On Ubuntu / Debian based Linux distributions, as `root` user, do:
     
            apt-get install libexpat1 libexpat1-dev
                    
     * On RedHat / Fedora / CentOS based Linux distributions, as `root` user do:
     
            yum install expat expat-devel
            
* Caveats:
    * The pipeline script uses a lot of inherent Linux core utils and has been only tested in BASH shell.
    * Please use absolute full PATH names. Instead of using `lncRNApipe --run ./lncRNApipe_output ...`, use
          `lncRNApipe --run /data/lncRNApipe_output ...`

 
           


=========
Cheers,

BioCoder
