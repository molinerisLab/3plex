#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
$,="\t";
$\="\n";

$SIG{__WARN__} = sub {die @_};

my $usage = "$0 [-h|help] [-l|length L] ssRNA.fa < cCRE.fa
count the frequency of each word of L length [default 6], indipendently for each fasta block,
if -c is not given count only exixting words,
if -c has a value like \"ACGT\" then count all possible words of len LEN composed of such characters,
	even if not present in the fasta block (0 frequency)
\n";

my $help=0;
my $len=6;
my $verbose=0;
GetOptions (
	'h|help' => \$help,
	'v|verbose' => \$verbose,
	'l|length=i' => \$len
) or die($usage);

if($help){
	print $usage;
	exit(0);
}

my $ssRNA="";
my $ssRNA_file_name = shift(@ARGV);
die $usage if not $ssRNA_file_name;

open(my $ssRNA_file, "<", $ssRNA_file_name) or die "Can't open file ($ssRNA_file_name): $!";
while(<$ssRNA_file>){
	next if m/^>/;
	chomp;
	$ssRNA.=$_;
}

my %ssRNA_words = words($ssRNA,$len);
if($verbose){
	my $i=0;
	for(keys(%ssRNA_words)){
		$i++;
		print $i,$_
	}
}



my $header=undef;
my $block="";
while(<>){
	chomp;
	next if length==0;
	if(/^>/){
		my $post_header=$_;
		if(defined($header)){
			print $header if hoogsteen_match($block);
			$block = "";
		}
		$header=$post_header;
		next;
	}
	$block.=$_;
}
if(length($block)>0){
	print $header if hoogsteen_match_both_strand($block);
}





sub words{
	my $s=shift;
	my $l=shift;
	my $L = length($s);
	my %a=();
	my $word = substr($s,0,$len);
	while(length($word)==$len){
		$a{$word}=1;
		substr($s,0,1,"");
		$word = substr($s,0,$len);
	}
	return %a
}

sub hoogsteen_match_both_strand{
	my $s_FORWARD=shift;
	my $L = length($s_FORWARD);
	my $s_REVCOMP=reverse($s_FORWARD);
	   $s_REVCOMP=~tr/ACGT/TGCA/;
	
	my $retval=hoogsteen_match_single_strand($s_FORWARD);
	if($retval){
		return $retval
	}else{
		return hoogsteen_match_single_strand($s_REVCOMP)
	}
}

sub hoogsteen_match_single_strand{
	my $s=shift;
	my $word = substr($s,0,$len);
	while(length($word)==$len){
		#print ">$word";
		if($word!~m/[CT]/){
			my $word_rule_f1=$word;
			my $word_rule_f2=$word;
			my $word_rule_r1=reverse($word);
			my $word_rule_r2=$word_rule_r1;
			$word_rule_f1=~tr/A/T/;   	#TA-U & CG-G          http://docs.google.com/presentation/d/1qDIVgZI-TZd9SLyOtSRUROcL0qhrEhp816pzFUmbnBQ/
			$word_rule_f2=~tr/AG/TC/;	#TA-U & CG-C+
			#$word_rule_r1=$word; 		#TA-A & CG-G
			$word_rule_r2=~tr/A/T/;		#TA-U & CG-G
			#print $word_rule_f1,$word_rule_f2,$word_rule_r1,$word_rule_r2;
			for($word_rule_f1,$word_rule_f2,$word_rule_r1,$word_rule_r2){
				if($ssRNA_words{$_}){
					print $_ if $verbose;
					return 1
				}
			}
		}
		substr($s,0,1,"");
		$word = substr($s,0,$len);
	}
	return 0
}
