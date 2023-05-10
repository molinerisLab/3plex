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

my %ssRNA_hoogsteen_words = ();
hoogsteen_words_both_strand($ssRNA);

if($verbose){
	my $i=0;
	while(my ($k,$v) = each(%ssRNA_hoogsteen_words)){
		$i++;
		print $i,$k,$v,
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
			print $header;
			match($block);
			$block = "";
		}
		$header=$post_header;
		next;
	}
	$block.=$_;
}
if(length($block)>0){
	print $header;
	match($block);
}




sub hoogsteen_words_both_strand{
	my $s_FORWARD=shift;
	my $L = length($s_FORWARD);
	my $s_REVCOMP=reverse($s_FORWARD);
	   $s_REVCOMP=~tr/ACGT/TGCA/;
	
	hoogsteen_words_single_strand($s_FORWARD,"F");
	hoogsteen_words_single_strand($s_REVCOMP,"R");
}

sub hoogsteen_words_single_strand{
	my $s=shift;
	my $FR=shift;
	my $l=shift;
	my $L = length($s);
	my $i=0;
	my $word = substr($s,0,$len);
	while(length($word)==$len){
		if($word!~m/[CT]/){
			my $word_rule_f1=$word;
			my $word_rule_f2=$word;
			my $word_rule_r1=reverse($word);
			my $word_rule_r2=$word_rule_r1;
			$word_rule_f1=~tr/A/T/;   	#TA-U & CG-G          http://docs.google.com/presentation/d/1qDIVgZI-TZd9SLyOtSRUROcL0qhrEhp816pzFUmbnBQ/
			$word_rule_f2=~tr/AG/TC/;	#TA-U & CG-C+
			#$word_rule_r1=$word; 		#TA-A & CG-G
			$word_rule_r2=~tr/A/T/;		#TA-U & CG-G
			my $rule=0;
			for($word_rule_f1, $word_rule_f2, $word_rule_r1, $word_rule_r2){
				$rule++;
				$ssRNA_hoogsteen_words{$_}="$word;$FR;$rule;$i"
			}
		}
		substr($s,0,1,"");
		$word = substr($s,0,$len);
		$i+=$len;
	}
}


sub match{
	my $s=shift;
	my $i=0;
	my $word = substr($s,0,$len);
	while(length($word)==$len){
		my $v = $ssRNA_hoogsteen_words{$word};
		if($v){
			print $v,$i
		}
		substr($s,0,1,"");
		$word = substr($s,0,$len);
		$i+=$len
	}
	return 0
}