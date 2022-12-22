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
			(my $s_rule_f, my $s_rule_r) = hoogsteen_RNA_full_single_strand($block);
			print $header."f\n".$s_rule_f;
			print $header."r\n".$s_rule_r;
			$block = "";
		}
		$header=$post_header;
		next;
	}
	$block.=$_;
}
if(length($block)>0){
	(my $s_rule_f, my $s_rule_r) = hoogsteen_RNA_full_single_strand($block);
	print $header."_f\n".$s_rule_f;
	print $header."_r\n".$s_rule_r;
}

sub hoogsteen_RNA_full_single_strand{
	my $s=shift;
	my $s_rule_f=$s;
	my $s_rule_r=reverse($s);
	$s_rule_f=~tr/ACGT/NGGA/;   	#TA-U & CG-G & CG-C+          http://docs.google.com/presentation/d/1qDIVgZI-TZd9SLyOtSRUROcL0qhrEhp816pzFUmbnBQ/
	$s_rule_r=~tr/ACGT/ANGA/;	#TA-A & CG-G TA-U
	return ($s_rule_f, $s_rule_r)
}
