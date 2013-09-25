#!/usr/bin/perl

# BioCoder
# 11/17/2011
# This script contains examples of different perl operations

use strict;
use warnings;

my $x = 27;
my $y = 10;
my $string_1 = 'My name is';
my $string_2 = 'BioCoder';

# Add
my $add = $x + $y;

# Substract
my $substract = $x - $y;

# Multiply
my $multiply = $x * $y;

# Divide
my $divide = $x / $y;

# Power
my $exp = $x ** $y;

# Reminder
my $rem = $x % $y;

# Print above results

print "\nSum: $x + $y = $add\n";
print "Substraction: $x - $y = $substract\n";
print "Multiplication: $x x $y = $multiply\n";
print "Division: $x / $y = $divide\n";
print "Reminder: Reminder of $x divided by $y = $rem\n";
print "Exp: $x to the power of $y = $exp\n\n";

# Increment and Decrement Operators

# Post increment the value of x
my $total_value = $x++;

print "\nTotal Value: $total_value\n";
print "\nIncremented Value of x: $x\n";

# Post decrement the value of x

my $post_decremented_value = $x--;

print "Post Decremented Value: $post_decremented_value\n";
print "\nx: $x\n";

# Pre Increment and Pre Decrement

my $pre_increment_value = ++$x;
print "Pre Increment Value: $pre_increment_value\n";
print "x: $x\n\n";

my $pre_decrement_value = --$x;
print "Pre Decrement value = $pre_decrement_value\n";
print "x: $x\n\n";

# Concatenation
print "$string_1\n";
print "$string_2\n";

my $concatenated_string = $string_1 . ' ' . $string_2;
print "\nConcatenated String is: $concatenated_string\n";




