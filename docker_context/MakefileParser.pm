#package MakefileParser;
#our @EXPORT = qw(readmakefile);
use warnings;
use strict;

my $makefile_name = "";
sub set_makefile_name{
	$makefile_name = shift;
}

sub readmakefile
{
	my $out='';
	my $fh=shift or die("readmakefile require a filehandler param");
	while(<$fh>){
		if( m/^-?include\s+([^\s]+)/){
			my $raw_path = $1;
			my $path = $raw_path;
			if($path =~ /\$\((PRJ|TASK)_ROOT\)/){
				my $prj_root = $makefile_name;
				$prj_root =~ s/^$ENV{'BIOINFO_ROOT'}//;
				$prj_root =~ s|([^/]+/[^/]+)/.*|$1|;
				$path =~ s/\$\(PRJ_ROOT\)/$ENV{'BIOINFO_ROOT'}\/$prj_root/;
				$path =~ s/\$\(TASK_ROOT\)/$ENV{'BIOINFO_ROOT'}\/$prj_root/;
				
			}elsif($path =~/\$\(BIOINFO_ROOT\)/ or $path =~ /\$\(BIOINFO_HOST\)/){
				$path =~ s/\$\(BIOINFO_ROOT\)/$ENV{'BIOINFO_ROOT'}/;
				$path =~ s/\$\(BIOINFO_HOST\)/$ENV{'BIOINFO_HOST'}/;
			}else{
				my $relative_path = $makefile_name;
				$relative_path =~ s/\/[^\/]+$/\//;
				$path = $relative_path.$path;
			}
			my $fh2;
			#if(!-e $path){
			#	my $dir = '.';
			#	my $targhet = $path;
			#	if($path =~ m|^(.*)/([^/]+)$|){
			#		$dir = $1;
			#		$targhet = $2;
			#	}
			#	my $cmd = "cd $dir; make -f $makefile_name $targhet";
			#	print "$cmd\n";
			#	print `$cmd`;
			#}
			if(-e $path){
				open $fh2,$path or die("Can't open $path");
				$out.=readmakefile($fh2);
				close $fh2;
			}else{
				if ($raw_path ne '$(BIOINFO_HOST).mkpaths') {
					print STDERR "WANRING: (MakefileParser) $path not found, ignoring.\n";
					$out.="$raw_path: \$(BIOINFO_HOST).mkpaths\n";
				}
				$out.=$_;
			}
		}else{
			$out.=$_;
		}
	}
	return $out;
}

1
