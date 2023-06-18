#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
$,="\t";


# TODO add tool general Description above Usage
# TODO add option '-h show this message and exit'
# TODO Fix better 'Option syntax: ..' # call it Option combos syntax?
# put it above Options section
# Option syntax: [-i|ignore-missing-columns] [-s|split_separator] [-w|split_each_word] [-a|append] [-d [-g glue] [-j]] [-m glue] [-v|allow_missing_translations [-e VALUE ]] [-k|kill] [-n|invert-dictionary] [-f|key_field N] [-p|pass REGEXP]
# TODO substitute column indication 1 2 3 with..?
#my $usage = "Usage: $0 [-i|ignore-missing-columns] [-s|split_separator] [-w|split_each_word] [-a|append] [-d [-g glue] [-j]] [-m glue] [-v|allow_missing_translations [-e VALUE ]] [-k|kill] [-n|invert-dictionary] [-f|key_field N] [-p|pass REGEXP] DICTIONARY 1 2 3 < tab_file\
my $usage = "\nUsage: translate [options] DICTIONARY 1 2 3 < TAB_FILE \

Options:
	-f N		specify what is the column (by its number) of DICTIONARY that contains the keys to be translated in TAB_FILE [default is column 1]
	-a		append translation: add new columns with translation on the right of the translated column (instead of substituting it)
	-r		put the added columns at the end of all TAB_FILE columns
	-i		ignore empty fields in column to be translated in TAB_FILE and keep the input row as is [NB: that row will have a different number of fields if used with -a option. Consider instead to suppress the row using -k]
	-k		kill untranslated rows (because of KEY missing in DICTIONARY)
	-d		allow duplicated keys in DICTIONARY
	  -g GLUE	indicate the separator for multiple values of the same key in output (when there are duplicated translations for the SAME KEY, list them in the same row, separated by GLUE)
	  -j		join like out (when there are duplicated translations for the SAME KEY, generate multiple rows clustered together, one for each translation)
	-m GLUE		if the DICTIONARY file has more than 2 columns, the translation is multi-column and the separator is GLUE (different from -j because double tanslations are on the same row here) [default glue is ';']
	-b FILE		print the killed rows of STDIN/TAB_FILE in FILE
	-v		allow missing translations in the dictionary (print also rows of TAB_FILE whose value to be translated is not in the dictionary)
	  -e VALUE	use VALUE as translation when a key is not present in dictionary [must be used together with -v]
	-p REGEXP	ignore input rows matching REGEXP (in Perl syntax), e.g. use -p '^>' to skip the translation of FASTA headers
	-c		case insensitive (ignore case in the comparison with the dictionary)
	-w		translate each word (word = each string of letters, numbers and/or underscores ('_'). Words can be separated by spaces (' '), hyphens ('-'), dots('.'), tabs etc [NOT compatible with column indication nor with options -a, -k]
	-s		split separator
	-z		allow dictionary to be empty
	-n 		invert columns of DICTIONARY (same as -f 2) [compatibile with -m]


";

# -s instead of tabs, split terms on.. #right?
#TODO: options: doubts:
#-z also non existent DICTIONARY ok?
#

# old help:
#Options:
#	-d	allow duplicated keys in DICTIONARY
#	-g	indicate the separator for multiple values of the same key in output
#	-m	if the dictionary file has more than 2 columns, the transaltion is multi column and the separator is glue
#	-k	kill untranslated rows
#	-b FILE print the killed rows of STDIN in FILE
#	-e	use VALUE as translation when a key si not present in dictionary
#	-j	join lyke out (when there duplicated translation for the same key generate multiple rows)
#	-p	ignore input rows matching REGEXP, use -p '^>' to skip the translation of fasta headers
#	-c	ignore case in the comparison with the dictionary
#	-w	translate each word
#	-r	put the added columns at the end of rows
#	-z	allow empty dictionary

my $key_field = undef;
my $ignore_missing_columns=0;
my $append=0;
my $duplicated_keys=0;
my $glue=undef;
my $split_separator="\t";
my $split_each_word=0;
my $kill=0;
my $print_killed=undef;
my $invert_dictionary=0;
my $allow_missing_translations=0;
my $join_lyke_out=0;
my $empty_key = undef;
my $multi_column_separator = undef;
my $skip_input_regexp = undef;
my $ignore_case=0;
my $append_at_the_R_of_rows=0;
my $allow_empty_dictionary=0;
GetOptions (
	'f|key_field=i' => \$key_field,
	'i|ignore-missing-columns' => \$ignore_missing_columns,
	'n|invert-dictionary' => \$invert_dictionary,
	'a|append' => \$append,
	'd|duplicated_keys' => \$duplicated_keys,
	's|split_separator=s' => \$split_separator,
	'e|empty=s' => \$empty_key,
	'w|split_each_word' => \$split_each_word,
	'g|glue=s' => \$glue,
	'j|join' => \$join_lyke_out,
	'k|kill' => \$kill,
	'b|print_killed=s' => \$print_killed,
	'v|allow_missing_translations' => \$allow_missing_translations,
	'm|multi_column_separator=s' => \$multi_column_separator,
	'p|pass=s' => \$skip_input_regexp,
	'c|case' => \$ignore_case,
	'r|append_at_the_R_of_rows' => \$append_at_the_R_of_rows,
	'z|allow_empty_dictionary' => \$allow_empty_dictionary
) or die($usage);

$SIG{__WARN__} = sub {die @_};

die("ERROR: -g option is meaningless without -d option") if defined $glue  and !$duplicated_keys;

$glue = ';' if not defined $glue;

#$allow_missing_translations = 1 if $kill;

my $filename = shift @ARGV;
if ($filename =~ /.gz$/) {
    open(FH, "zcat $filename|") or die("ERROR: Can't open file ($filename)");
} elsif ($filename =~ /.xz$/) {
    open(FH, "xzcat $filename|") or die("ERROR: Can't open file ($filename)");
} else {
    open(FH,$filename) or die("ERROR: Can't open file ($filename)");
}

open KILLED,">$print_killed" or die("ERROR: Can't open file ($print_killed)") if $print_killed;

my @columns=@ARGV;

die("ERROR: no column indication\n$usage") if scalar(@columns) == 0 and !$split_each_word;

for(@columns){
	die("ERROR: invalid columns ($_)") if !m/^\d+$/;
	$_--;
}

die("ERROR: -w option incompatible with column indication") if($split_each_word and scalar(@columns) > 0);
die("ERROR: -j option requires -d and -a options and conflicts with -w") if $join_lyke_out and ($split_each_word or (!$duplicated_keys or !$append));
die("ERROR: -w option conflicts with -k") if ($split_each_word and $kill);
die("ERROR: -b|print_killed option require -k") if ($print_killed and not $kill);
die("ERROR: -e meaningless without -v") if defined($empty_key) and not $allow_missing_translations;
die("ERROR: -n meaningless with -f") if defined($key_field) and $invert_dictionary;
die("ERROR: -f require a parameter >=1") if defined($key_field) and $key_field < 1;
die("ERROR: -r require -a") if $append_at_the_R_of_rows and not $append;
die("ERROR: -r not compatible with -j") if $append_at_the_R_of_rows and $join_lyke_out;

if (!defined($multi_column_separator)) {
    $multi_column_separator = "\t";
}

$key_field-- if defined($key_field);
$key_field = 0 if !defined($key_field);




my $columns_added_by_translation = undef;
my %hash=();
my $empty_map = 1;
while(<FH>){
	$empty_map = 0;
	chomp;

	my $k = undef;
	my $v = undef;
	if($key_field == 0){
		die("ERROR: At least 2 columns required in dictionary file (".$filename.")") if !m/\t/;
		m/([^\t]+)[\t](.*)/;
		if(!$invert_dictionary){
			$k = $1;
			$v = $2;
		}else{
			$v = $1;
			$k = $2;
		}
		if(defined($multi_column_separator)){
			$v =~ s/\t/$multi_column_separator/g;
		}
		die("ERROR: Malformed input in dictionary ($_)") if not defined $v;
		if(not defined $columns_added_by_translation){
			my @F=split(/\t/, $v, -1);
			$columns_added_by_translation = scalar(@F);
		}
	}else{
		my @F = split(/\t/, $_, -1);
		$k = splice @F, $key_field, 1;
		$v = join( $multi_column_separator, @F);
		if(defined($columns_added_by_translation)){
			warn "WARNING: The dictionary file has rows with different number of fields." if scalar @F != $columns_added_by_translation;
		}
		$columns_added_by_translation = scalar(@F) if not defined $columns_added_by_translation;
	}

	$k = uc($k) if $ignore_case;
	if(not defined($k)){
		die("ERROR: Malformed input in dictionary ($_)");
	}

	if(defined $hash{$k}){
		if(!$duplicated_keys){
			die("ERROR: Duplicated key in dictionary ($k)");
		}else{
			if($join_lyke_out){
				if(ref($hash{$k}) eq 'ARRAY'){
					push(@{$hash{$k}},$v);
				}else{
					my @tmp=($hash{$k},$v);
					$hash{$k}=\@tmp;
				}
			}else{
				$hash{$k}.=$glue.$v;
			}
		}
	}else{
		$hash{$k}=$v;
	}
}

if($empty_map){
	if($allow_empty_dictionary){
		$columns_added_by_translation=1
	}else{
		if($kill){
			#return(0); # will raise a broken pipe
			while(<STDIN>){
				#consume input
			}
			exit(0);
		}
		die("WARNING: The dictionary is empty");
	}
}

if($columns_added_by_translation > 1){
	if(not defined $empty_key){
		$empty_key = "\t" x ($columns_added_by_translation-1) if $allow_missing_translations and not defined $empty_key;
	}else{
		$empty_key = "\t$empty_key" x $columns_added_by_translation if $allow_missing_translations;
		$empty_key =~ s/^\t//;
	}
}else{
	$empty_key = "" if not defined $empty_key;
}



my $warning=0;

while(<STDIN>){

	if(defined $skip_input_regexp and m/$skip_input_regexp/){
		print;
		next;
	}

	if(!$split_each_word){
		chomp;
		my @F = split /$split_separator/;
		my @G=@F;
		my $print=1;
		for(@columns){
			$a=$F[$_];
			if(defined($a)){
				my ($val, $translated)=@{&translate($a)};
				$print = 0 if (!$translated and $kill);
				if( not $join_lyke_out){
					if(not $append_at_the_R_of_rows){
						$F[$_] = $val;
					}else{
						push(@F,$val);
					}
				}else{ 
					#$append_at_the_R_of_rows not allowed in -j mode
					my @tmp= ($F[$_],$val); 
					$F[$_] = \@tmp; #nella colonna da tradurre metto una ref ad un array col valore iniziale
							#e quello tradotto (che sara` quello iniziale tab traduzione in caso di -a
							# e non join_lyke_out e solo la traduzione negli altri casi)
				}
			}else{
				die("ERROR: column $_+1 not defined in standard input (-i to ignore)") if !$ignore_missing_columns;
				if(!$warning){
					print STDERR "WARNING: $0, column not defined\n";
					$warning=1;
				}
			}
		}
		if($print){
			if(!$join_lyke_out){
				print @F;
				print "\n";
			}else{
				die("ERROR: only one column is allowed when using option '-j' (used to gather multiple rows in output, each with a translation for the same key)") if scalar(@columns) > 1;
				# original message: "only one column allowed when join like output enabled"
				my $c = $columns[0];
				my $a = $G[$c];
				my $val=$F[$c]->[1];
				if(ref($val) ne 'ARRAY'){
					$F[$c] = $a .$split_separator.$val;
					print @F;
					print "\n";
				}else{
					for(@{$val}){
						$F[$c] = $a .$split_separator.$_;
						print @F;
						print "\n";
					}
				}
			}
		}else{
			if($print_killed){
				#print KILLED @F;
				#print KILLED "\n";
				if(!$join_lyke_out){
					print KILLED @F;
					print KILLED "\n";
				}else{ 
					die("ERROR: only one column is allowed when using option '-j' (used to gather multiple rows in output, each with a translation for the same key)") if scalar(@columns) > 1;
					my $c = $columns[0];
					my $a = $G[$c];
					my $val=$F[$c]->[1];
					if(ref($val) ne 'ARRAY'){
						$F[$c] = $val;
						print KILLED @F;
						print KILLED "\n";
					}else{
						for(@{$val}){
							$F[$c] = $split_separator.$_;
							print KILLED @F;
							print KILLED "\n";
						}
					}
				}

			}
		}
	}else{
		s/^([\W]+)//;
		print $1 if $1;
		while(s/([^\W]+)([\W]+)//){
			my ($val, $translated)=@{&translate($1)};
			print $val;
			print $2;
		}
	}
}

sub translate
{
	my $a=shift;
	my $b=undef;
	if($ignore_case){
		$b=$hash{uc($a)};
	}else{
		$b=$hash{$a};
	}
	if(!defined($b)){
		if($allow_missing_translations){
			if(!$split_each_word){
				if(defined($empty_key)){
					$b=$empty_key
				}
			}else{
				$b=$a
			}
		}elsif(not $kill){
			warn "WARNING: missing translation for key ($a)\n" if !$allow_missing_translations;
		}
	}
	if(defined($b)){
		if($append and not $join_lyke_out and not $append_at_the_R_of_rows){
			$a .= "\t$b";
		}else{
			$a = $b;
		}
	}

	my @tmp=($a, defined($b));
	return \@tmp;
}
