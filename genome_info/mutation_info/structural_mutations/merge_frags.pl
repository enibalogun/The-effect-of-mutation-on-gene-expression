#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use Data::Dumper;

open (IN1, "s_CC4532_muts.bed") or die;

my %lift;

while (my $l1 = <IN1>) {
	chomp $l1;
	my @c1 = split(/\t/, $l1);
	push @{ $lift{$c1[3]} }, $l1;
}

close IN1;

open (OUT, ">s_CC4532_muts.10bp.bed") or die;

foreach my $mut (keys %lift) {
	open (OUT1, ">temp.bed") or die;
    foreach my $frags (@{$lift{$mut}}) {
        print OUT1 "$frags\n";
    }
    close OUT1;
    system("bedtools merge -i temp.bed -d 10 -s -c 6 -o distinct > merged.bed");
    open (IN2, "merged.bed") or die;
    while (my $l2 = <IN2>) {
    	chomp $l2;
    	my @c2 = split(/\t/, $l2);
    	print OUT "$c2[0]\t$c2[1]\t$c2[2]\t$mut\t0\t$c2[3]\n";
    }
    close IN2;
    system("rm temp.bed");
    system("rm merged.bed")
}

close OUT;

exit;
