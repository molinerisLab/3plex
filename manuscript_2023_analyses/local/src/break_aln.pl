#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
$,="\t";
$\="\n";

$SIG{__WARN__} = sub {die @_};

my $usage = "$0 [-h|help] [-g glue]]\n";

my $help=0;
my $max_consecutive_errors=1;
my $min_len=8;
my $error_rate=0.2;
GetOptions (
	'h|help' => \$help,
	'c|max_consecutive_errors=i' => \$max_consecutive_errors,
	'l|min_len=i' => \$min_len,
	'e|error_rate=f' => \$error_rate
) or die($usage);

if($help){
	print $usage;
	exit(0);
}

my $break_at=$max_consecutive_errors+1;

while(<>){
	chomp;
	(	my $query_id,
		my $target_id,
		my $query_start,
		my $query_end,	
		my $target_start,
		my $target_end,
		my $query_strand,
		my $target_strand,	
		my $query_fragment,	
		my $target_fragment, 
		my $similarity
	) = split /\t/,$_,-1;

	print ">$_";
	my @tpxs = split /(\s{$break_at,})/,$similarity; 
	my $pos=0;
	for(@tpxs){
		my $l=length($_);
		my $match = $_ =~ tr/://;
		my $e=($l-$match)/$l;
		if($e<=$error_rate and $l>=$min_len){
			print 	$query_id,
				$query_start+$pos,
				$query_start+$pos+$l,
				$target_id,
				$target_start+$pos,
				$target_start+$pos+$l,
				$match,$l,$e,$_,$pos;
		}
		$pos+=length($_)
	}
}
