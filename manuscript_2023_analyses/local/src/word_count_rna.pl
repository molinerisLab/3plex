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
	print join("\n",hoogsteen_full_both_strand($ssRNA));
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
			print $header if match($block);
			$block = "";
		}
		$header=$post_header;
		next;
	}
	$block.=$_;
}
if(length($block)>0){
	print $header if match($block);
}




sub hoogsteen_words_both_strand{
	my $s_FORWARD=shift;
	my $L = length($s_FORWARD);
	my $s_REVCOMP=reverse($s_FORWARD);
	   $s_REVCOMP=~tr/ACGT/TGCA/;
	
	hoogsteen_words_single_strand($s_FORWARD);
	hoogsteen_words_single_strand($s_REVCOMP);
}

sub hoogsteen_words_single_strand{
	my $s=shift;
	my $l=shift;
	my $L = length($s);
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
			my $i=0;
			for($word_rule_f1, $word_rule_f2, $word_rule_r1, $word_rule_r2){
				$i++;
				#$ssRNA_hoogsteen_words{$_}=$word."_$i"
				$ssRNA_hoogsteen_words{$_}=1
			}
		}
		substr($s,0,1,"");
		$word = substr($s,0,$len);
	}
}

sub hoogsteen_full_both_strand{
	my $s_FORWARD=shift;
	my $s_REVCOMP=reverse($s_FORWARD);
	   $s_REVCOMP=~tr/ACGT/TGCA/;
	return (hoogsteen_full_single_strand($s_FORWARD),hoogsteen_full_single_strand($s_REVCOMP))
	
}

sub hoogsteen_RNA_full_single_strand{
	my $s=shift;
	my $s_rule_f1=$s;
	my $s_rule_f2=$s;
	my $s_rule_r1=reverse($s);
	my $s_rule_r2=$s_rule_r1;
	$s_rule_f1=~tr/ACGT/NNGA/;   	#TA-U & CG-G          http://docs.google.com/presentation/d/1qDIVgZI-TZd9SLyOtSRUROcL0qhrEhp816pzFUmbnBQ/
	$s_rule_f2=~tr/ACGT/NGNA/;	#TA-U & CG-C+
	$s_rule_r1=~tr/ACGT/ANGN/;	#TA-A & CG-G
	$s_rule_r2=~tr/ACGT/NNGA/;	#TA-U & CG-G
	return ($s_rule_f1, $s_rule_f1, $s_rule_r1, $s_rule_r2)
}


sub match{
	my $s=shift;
	my $word = substr($s,0,$len);
	while(length($word)==$len){
		return 1 if $ssRNA_hoogsteen_words{$word};
		substr($s,0,1,"");
		$word = substr($s,0,$len);
	}
	return 0
}