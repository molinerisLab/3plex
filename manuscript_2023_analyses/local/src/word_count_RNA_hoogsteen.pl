#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
$,="\t";
$\="\n";

$SIG{__WARN__} = sub {die @_};

my $usage = "$0 [-h|help] < ssRNA.fa 
count the frequency of each word of L length [default 6], indipendently for each fasta block,
if -c is not given count only exixting words,
if -c has a value like \"ACGT\" then count all possible words of len LEN composed of such characters,
	even if not present in the fasta block (0 frequency)
\n";

my $help=0;
my $verbose=0;
GetOptions (
	'h|help' => \$help,
) or die($usage);

if($help){
	print $usage;
	exit(0);
}

my $ssRNA="";
my $header=undef;
my $block="";
while(<>){
	chomp;
	next if length==0;
	if(/^>/){
		my $post_header=$_;
		if(defined($header)){
			(my $s_rule_MP, my $s_rule_YP, my $s_rule_RA, my $s_rule_MA) = hoogsteen_RNA_full_single_strand($block);
			print $header."MP\n".$s_rule_MP;
			print $header."YP\n".$s_rule_YP;
			print $header."RA\n".$s_rule_RA;
			print $header."MA\n".$s_rule_MA;
			$block = "";
		}
		$header=$post_header;
		next;
	}
	$block.=$_;
}
if(length($block)>0){
	(my $s_rule_MP, my $s_rule_YP, my $s_rule_RA, my $s_rule_MA) = hoogsteen_RNA_full_single_strand($block);
	print $header."_MP\n".$s_rule_MP;
	print $header."_YP\n".$s_rule_YP;
	print $header."_RA\n".$s_rule_RA;
	print $header."_MA\n".$s_rule_MA;
}

sub hoogsteen_RNA_full_single_strand{
	my $s=shift;
	my $s_rule_MP=$s;
	my $s_rule_YP=$s;
	my $s_rule_RA=reverse($s);
	my $s_rule_MA=$s_rule_RA;

	# http://docs.google.com/presentation/d/1qDIVgZI-TZd9SLyOtSRUROcL0qhrEhp816pzFUmbnBQ/
	#
	# [T,C] (pyrimidine motif)  		Y
	# [G,A] (purine motif)      		R 
	# [G,T] (purine–pyrimidine motif) 	M 
	# Triplexator paper and (Morgan and Wells 1968; Cooney et al. 1988; Beal and Dervan 1991). 
	# 
	# In the pyrimidine motif, 
	# T:AT and C+:GC 
	# triads are formed in Hoogsteen configuration, 
	# (the cytosine of the third strand needs to be protonated 
	# in order to form the second hydrogen bond, acidic pH)
	#
	# In the purine motif, reverse Hoogsteen bonds are formed
	# G:GC and A:AT triads
	# 
	# [G,T] motif allows a mixed purine–pyrimidine TFO and forms 
	# G:GC and T:AT triads 
	# in either Hoogsteen or reverse Hoogsteen configuration. 
	#
	# (T refers to uracil in case a strand is made of RNA.)
	#
	#
	

	$s_rule_MP=~tr/ACGT/NNGA/;	# M P 	TA-U & CG-G   
	$s_rule_YP=~tr/ACGT/NGNA/;	# Y P	TA-U & CG-C+
	$s_rule_RA=~tr/ACGT/ANGN/;	# R A	TA-A & CG-G
	$s_rule_MA=~tr/ACGT/NNGA/;	# M A	TA-U & CG-G
	return ($s_rule_MP, $s_rule_YP, $s_rule_RA, $s_rule_MA)
}
