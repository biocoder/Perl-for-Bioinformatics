#!/usr/bin/env perl

# BioCoder
# May-8-2012
# This script introduces control structures, arrays and hashes.

use strict;
use warnings;

my @first_array = ('DNA', 'ATGCGTGC', 5, 'RNA', 'AUGC');
print $first_array[0], "\n\n";

# Scalar gives actual size of an array

my $size_of_array = scalar(@first_array);
my $another_way_of_getting_size_of_array = @first_array; # implicit way
print "Scalar of array: $size_of_array\n\n";
print "Perl's index size of array: $#first_array\n\n";
print "Another way of getting size: $another_way_of_getting_size_of_array\n\n";

# Control Loop: for

for (my $i=0; $i<=$#first_array; $i++) {
    print "Perl's array index: $i\n\n";
    print "$first_array[$i] \n\n"; 
}

my %sequence = ('DNA' => 'ATCGATGCT',
                'RNA' => 'AUGC',
                'Number of seqs' => 2
                );

print $sequence{'DNA'}, "\n";
# Control Loop: foreach

foreach my $key (sort (keys %sequence)) {
    print "Key of hash is $key\tValue of hash is $sequence{$key}\n\n";
}