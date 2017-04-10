#!/usr/bin/perl

use warnings;
use strict;

#Command line processing.
use Getopt::Long;

my $arbre;
my $param;
my $length;
Getopt::Long::Configure ('bundling');
GetOptions ('a|arbre=s' => \$arbre,
	    'p|param=s' => \$param,
	    'l|length=i' => \$length)
or die("Error in command line arguments\n");


open(F2,"$param") or die ("Erreur de creation de F2\n");
open(F3,">param.output") or die ("Erreur de creation de F3\n");
while (my $li = <F2>){
	chomp($li);
	if ($li =~ m/^SequenceNum=\d+/){
		$li =~ s/.+/SequenceNum=$length/g;
		print F3 "$li\n";	
	}
	elsif ($li =~ m/^TheTree=\(.+\;/){
		open(F1, "$arbre") or die("erreur ouverture fichier inputFile\n");
		while (my $li1 = <F1>){
			chomp($li1);
			$li =~ s/.+/TheTree=$li1/g;
			print F3 "$li\n";
		}
		close F1;
	}
	else{
		print F3 "$li\n";
	}	
}
close F2;
close F3;
