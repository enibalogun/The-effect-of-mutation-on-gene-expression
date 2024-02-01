#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use Data::Dumper;

open (IN, "anc_muts.tsv") or die;
open (OUT, ">anc_muts.bed") or die;

while (my $line = <IN>) {
	chomp $line;
	my @cols = split(/\t/, $line);
	my $chr;
	my $id;
	if ($cols[0] < 10) {
		$chr = "chromosome_0" . $cols[0];
	}
	else {
		$chr = "chromosome_" . $cols[0];
	}
	if ($cols[5] eq "Insertion Mobile") {
		my $s1 = $cols[1] - 5;
		my $e1 = $cols[2] + 5;
		$id = $cols[4] . "_insertion_" . $cols[6] . "_" . $chr . "_" . $cols[1];
		print OUT "$chr\t$s1\t$e1\t$id\t0\t+\n";
	}
	elsif ($cols[5] eq "Excision Mobile") {
		$id = $cols[4] . "_excision_" . $cols[6] . "_" . $chr . "_" . $cols[1] . "_" . $cols[2] . "_left";
		my $s = $cols[1] - 5;
		print OUT "$chr\t$s\t$cols[1]\t$id\t0\t+\n";
		$id = $cols[4] . "_excision_" . $cols[6] . "_" . $chr . "_" . $cols[1] . "_" . $cols[2] . "_right";
		my $e = $cols[2] + 5;
		print OUT "$chr\t$cols[2]\t$e\t$id\t0\t+\n";
	}  
	elsif ($cols[5] eq "Duplication") {
		$id = $cols[4] . "_duplication_" . $chr . "_" . $cols[1] . "_" . $cols[2];
		print OUT "$chr\t$cols[1]\t$cols[2]\t$id\t0\t+\n";
	}
	elsif ($cols[5] eq "Deletion") {
		$id = $cols[4] . "_deletion_" . $chr . "_" . $cols[1] . "_" . $cols[2];
		print OUT "$chr\t$cols[1]\t$cols[2]\t$id\t0\t+\n";
	}
	elsif ($cols[5] eq "Inversion") {
		$id = $cols[4] . "_inversion_" . $chr . "_" . $cols[1] . "_" . $cols[2] . "_left";
		my $s4 = $cols[1] - 5;
		print OUT "$chr\t$s4\t$cols[1]\t$id\t0\t+\n";
                $id = $cols[4] . "_inversion_" . $chr . "_" . $cols[1] . "_" . $cols[2] . "_right";
                my $e4 = $cols[2] + 5;
                print OUT "$chr\t$cols[2]\t$e4\t$id\t0\t+\n";

	}
	elsif ($cols[5] eq "Translocation") {
		$id = $cols[4] . "_translocation_" . $chr . "_" . $cols[1];
		my @i = split(/\:/, $cols[2]);
		my $s2 = $cols[1] - 5;
		my $e2 = $cols[1] + 5;
		print OUT "$chr\t$s2\t$e2\t$id\t0\t+\n";
		$id = $cols[4] . "_translocation_" . $i[0] . "_" . $i[1];
		my $s3 = $i[1] - 5;
		my $e3 = $i[1] + 5;
		print OUT "$i[0]\t$s3\t$e3\t$id\t0\t+\n";
	}
	else {
		print "$cols[5]\n";
	}
}

close IN;
close OUT;

exit;
