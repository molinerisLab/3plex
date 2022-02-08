#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
$,="\t";
$\="\n";

$SIG{__WARN__} = sub {die @_};

my $usage = "$0 -c 1 -e 0.2 -l 10 \n";

my $help=0;
my $max_consecutive_errors=1;
my $min_len=8;
my $error_rate=0.2;
my $verbose=0;
my $no_trim=0;
GetOptions (
	'h|help' => \$help,
	'c|max_consecutive_errors=i' => \$max_consecutive_errors,
	'l|min_len=i' => \$min_len,
	'e|error_rate=f' => \$error_rate,
	'v|verbose' => \$verbose,
	'n|no_trim' => \$no_trim

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
	#if($query_strand==-1){
	#	$pos=$query_end;
	#}

	for(@tpxs){
		my $l=length($_);
		my $match = $_ =~ tr/://;
		my $e=($l-$match)/$l;
		if($l>=$min_len){
			my $tfo_start = -1;
			my $tfo_end   = -1;
			my $tts_start = -1;
			my $tts_end   = -1;
			if($query_strand>0){
				$tfo_start = $query_start+$pos;
				$tfo_end   = $query_start+$pos+$l;
				$tts_start = $target_start+$pos;
				$tts_end   = $target_start+$pos+$l;
			}else{
				$tfo_start = $query_end-($pos+$l);
				$tfo_end   = $query_end-$pos;
				$tts_start = $target_start+$pos;
				$tts_end   = $target_start+$pos+$l;
			}
			if($e<=$error_rate){
				print 	$query_id,
					$tfo_start,
					$tfo_end,
					$target_id,
					$tts_start,
					$tts_end,
					$match,$l,$e,$_,$pos;
			}elsif(not $no_trim){
				trim_error_rate(
					$query_id,
                                        $tfo_start,
                                        $tfo_end,
                                        $target_id,
                                        $tts_start,
                                        $tts_end,
                                        $match,$l,$e,$_,$pos
				);
				
			}
		}
		$pos+=length($_)
	}
}


sub trim_error_rate{
	my $query_id = shift;
 	my $tfo_start = shift;
	my $tfo_end = shift;
	my $target_id = shift;
	my $tts_start = shift;
	my $tts_end = shift;
	my $match = shift;
	my $l = shift;
	my $e = shift;
	my $_ = shift;
	my $pos = shift;

	#trim_error_rate_F(0,0,$_);
	#trim_error_rate_R(0,0,$_);
	my $trimmed=trim_error_rate_norecursion($_);
	print "best",@{$trimmed}
}


#	my $best_side=undef;
#	if(defined($F_e) and defined($R_e)){
#		$best_side="F";
#		if($R_e>$F_e){
#			$best_side="R";
#		}
#	}elsif(defined($F_e)){
#		$best_side="F";
#	}elsif(defined($R_e)){
#		$best_side="R";
#	}
#	
#	
#	if($best_side eq "F" and $query_strand>0){
#		$tfo_star+=$F_shift;
#		$
#	}


sub trim_error_rate_norecursion{
	my $_ = shift;
	my $l = length($_);
	
	my @trimming=();

	for(my $start=0; $start<=$l-$min_len; $start++){
		next if substr($_,$start,1) eq " ";
		for(my $end=$l; $end>=$min_len; $end--){
			my $s=$_;
			$s=substr $_, $start, $end-$start;
			next if substr($s,-1,1) eq " ";
			my $match = () = $s =~ /:/g;
			my $s_l=length($s);
			$start+=length($1) if(s/^(\s+)//g);
			$end  -=length($1) if(s/(\s+)$//g);
			last if $end - $start < $min_len;

			my $e=($s_l-$match)/$s_l;
			my @retval=($start,$end,"($s)",$s_l,$e);
			print @retval if $verbose;

			if($e<=$error_rate){
				push @trimming, \@retval;
			}
		}
	}

	my $t=scalar(@trimming);
	if($t==0){
		return undef;
	}
	if($t>1){
		@trimming = sort { $b->[3] <=> $a->[3] } @trimming
	}

	return $trimming[0]
}

sub trim_error_rate_F{
	my $start_shift = shift;
	my $end_shift = shift;
	my $_ = shift;

	my @trimming=();

	s/:/ /;
	s/^(\s+)//g;
	die("someting wrong in recursion") if not $1;
	$start_shift+=length($1);
	my $match = () = $_ =~ /:/g;
	my $l=length;

	return if $l<=$min_len;

	my $e=($l-$match)/$l;
	if($e<=$error_rate){
		my @retval=($start_shift,$end_shift,$_,$l,$e);
		print @retval if $verbose;
		push @trimming, \@retval;
		return
	}else{
		trim_error_rate_F($start_shift,$end_shift,$_);
		trim_error_rate_R($start_shift,$end_shift,$_);
	}
}

sub trim_error_rate_R{
	my $start_shift = shift;
	my $end_shift = shift;
	my $_ = shift;
	$_=reverse($_);
	s/:/ /;
	$_=reverse($_);
	s/(\s+)$//g;
	$end_shift+=length($1);
	die("sometin wrong in recursion") if not $1;
	my $match = () = $_ =~ /:/g;
	my $l=length;

	return if $l<=$min_len;

	my $e=($l-$match)/$l;
	if($e<=$error_rate){
		my @retval=($start_shift,$end_shift,$_,$l,$e);
		print @retval if $verbose;
		push @trimming, \@retval;
		return
	}else{
		trim_error_rate_F($start_shift,$end_shift,$_);
		trim_error_rate_R($start_shift,$end_shift,$_);
	}
}
