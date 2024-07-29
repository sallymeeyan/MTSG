#!/usr/bin/env perl
use strict;
use Getopt::Long;
use Bundle::Wrapper;
use Data::Dumper;

my $bundle = Bundle::Wrapper->new();
my $opt_h;
my $opt_chr;
my $opt_start;
my $opt_end;
my $opt_db;
my $opt_db_pattern;
my $opt_query;
my $opt_delimiter="\t";
my $opt_silent;
my $opt_outmatch;
my $opt_header;
my $opt_exact;
my $opt_excute;
my $opt_norecord;
my $opt_seperated;
my $opt_chrnum;

GetOptions(
    "help|h" => \$opt_h,
    "chr=i" => \$opt_chr,
    "start=i" => \$opt_start,
    "end=i" => \$opt_end,
    "db=s" => \$opt_db,
    "db_pattern=s" => \$opt_db_pattern,
    "query=s" => \$opt_query,
    "delimiter=s" => \$opt_delimiter,
    "silent!" => \$opt_silent,
    "seperated!" => \$opt_seperated,
    "outmatch=s" =>\$opt_outmatch,
    "header=s" => \$opt_header,
    "norecord!" => \$opt_norecord,
    "exact!" =>\$opt_exact,
    "chrnum!" =>\$opt_chrnum,
    "e|excute!" => \$opt_excute,
);

if(!$opt_excute || $opt_h){
    &printopts2(
	 {
	     "help|h" => \$opt_h,
		 "chr=i" => \$opt_chr,
		 "start=i" => \$opt_start,
		 "end=i" => \$opt_end,
		 "db=s" => \$opt_db,
		 "db_pattern=s" => \$opt_db_pattern,
		 "query=s" => \$opt_query,
		 "delimiter=s" => \$opt_delimiter,
		 "silent!" => \$opt_silent,
		 "seperated!" => \$opt_seperated,
		 "outmatch=s" =>\$opt_outmatch,
		 "header=s" => \$opt_header,
		 "norecord=s" => \$opt_norecord,
		 "exact!" =>\$opt_exact,
		 "e|excute!" => \$opt_excute,
	 }
	);
    print "type -e to run\n";
    exit 0;
}

if ($opt_query){
    open (QUERY,"<$opt_query");
}
else{
    die "Fail open input file";
}

open(OM,">$opt_outmatch"), if $opt_outmatch;

my $head1=`vcfhead $opt_db`;
my $head2=`vcfhead $opt_db 2`; 
if($opt_header==1){ print $head1; }
elsif($opt_header==2){print $head2;}


my $count;
while(my $q=<QUERY>){
    chomp $q;
    my $result;
    
    my @entry = split/\t/,$q;
    
    $entry[$opt_chr-1]=~s/[ \t]//g;
    $entry[$opt_start-1]=~s/[ \t]//g;
    $entry[$opt_end-1]=~s/[ \t]//g;

    #print Data::Dumper @entry;
    ##Skip the non-digital line
    if($entry[$opt_start-1]!~/\D/)
    {
	##tabix core cmd;
	my $tmp= $entry[$opt_chr-1];
	my $db=$opt_db;
	if($opt_seperated){
	    #print STDERR $tmp,"\n";
	    #print STDERR $opt_db_pattern,"\n";
	    $db=~s/$opt_db_pattern/$1$tmp$2/;
	    print STDERR $db,"\n";
	}

	if($opt_chrnum){
	    $entry[$opt_chr-1]=~s/chr//i;
	}
	my $string = "$entry[$opt_chr-1]:$entry[$opt_start-1]-$entry[$opt_end-1]";

	## print STDERR "tabix $opt_db $string\n";
	print STDERR "tabix $db $string\n";
	my $result = `tabix $db $string`;

	#my $result = `wc -l tmpout.tabix`;
	## specify the delimiter you wanted.
	if($opt_delimiter ne "\t"){ $result =~s/\t/$opt_delimiter/g;}
	
	if($result)
	{
	    #Print matched query lines.
	    if($opt_outmatch){ print OM "$q\n"; }
	    
	    open(TMP,"<",\$result);
	    while(<TMP>)
	    {
		#Only output variants which have exact coordination.
		if($opt_exact)
		{
		    if(/$entry[$opt_chr-1]\t$entry[$opt_start-1]/){
			if(! $opt_silent)
			{
			    print"$q\t$_";
			}
			else
			{
			    print $_;
			}
		    }
		}
		#choose to silent $q parts in the result
		else 
		{
		    if(! $opt_silent)
		    {
			print"$q\t$_";
		    }
		    else
		    {
			print $_;
		    }	
		}
		
	    }
	    close TMP;
	}
	else{
	    if(!$result && $opt_norecord){
		$count++;
		#print STDERR "$q\tno found\n";
		if(! $opt_silent)
		{
		    print"$q\tnocomm$count\n";
		}
	    }
	}

	#system("rm tmpout.tabix");
	#To make sure every query has a result.
        #else
        #{
        #print $q,"\tnoprobe\n";
        #}
    }
}


sub printopts2
{
    my $class=shift;
    my $opt=$_[0];
    for (sort keys %$opt)
    {
	my $description = $opt->{$_};
	my $format='%-30s%-20s';
	if(!$$description)
	{
	    $$description='';
	}
	
	print "\t";
	printf $format,$_,$$description;
	print "\n";
    }
}
