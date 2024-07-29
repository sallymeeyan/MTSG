#!/usr/bin/env perl
use strict;
use Getopt::Long;
use Array::Transpose;
use Data::Dumper;

my @array;
my $opt_delim="\t";

GetOptions(
    "d|delimiter" =>\$opt_delim,
    );

while(<>){
    chomp;
    push @array,[(split/$opt_delim/,$_,-1)];
}

my @transposed=transpose(\@array);
foreach (@transposed){
    print join $opt_delim,@{$_};
    print "\n";
}
