#!/bin/env perl
use strict;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Storable;
use Switch;
use Bundle::Wrapper;
use Cwd;
use String::Random;
use List::Util qw /sum/;

####################################################################
##                             OPTIONS
####################################################################

my $Version=Bundle::Wrapper->date;
my $bundle=Bundle::Wrapper->new();
my %opt;
$opt{time}=Bundle::Wrapper->date;
$opt{threads}=5;
$opt{gene}="ENSG00000000419";
$opt{file_genePos}="/nobackup/cgg/yany14/MSG/data/gene.500k.id";
$opt{file_snp}="/nobackup/cgg/yany14/MSG/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz";
$opt{file_sample}="/nobackup/cgg/yany14/MSG/data/samples.filtered.used";
$opt{file_rsid}="/nobackup/cgg/yany14/MSG/data/hg38.vcf";
$opt{file_splicing}="/fs0/yany14/MSG/data/rawSplice/*Whole_Blood*.gz";
$opt{file_ldRef}="/nobackup/cgg/yany14/MSG/data/LDREF/1000G.EUR.";
$opt{dir_script}="/fs0/yany14/yany14/myscript/Rscripts/rmd/cca_multiple";
$opt{sumstats}="/fs0/yany14/MSG/data/sumstats/clozukscz.sumstats";

GetOptions(
    "help|h!" => \$opt{help},
    "excute|e!" => \$opt{e},
    "chkopt!" => \$opt{chkopt},
    "version|v!" => \$opt{version},
    "s|save!" => \$opt{save},
    "l|load!" => \$opt{load},
    "cmd=s" => \@{$opt{cmd}},
    "t|threads=i" => \$opt{threads},
    "temp!" => \$opt{temp},
    "i|in=s" => \$opt{in},
    "o|out=s" => \$opt{out},
    "header!" => \$opt{header},

    "dir_out=s" => \$opt{dir_out},
    "dir_log=s" => \$opt{dir_log},

    "gene=s" => \$opt{gene},
    "tissue=s" => \$opt{tissue},
    "file_splicing=s" => \$opt{file_splicing},
    "file_genePos=s" => \$opt{file_genePos},
    "file_snp=s" => \$opt{file_snp},
    "file_sample=s" => \$opt{file_sample},
    "file_rsid=s" => \$opt{file_rsid},
    "file_xmatrix=s" => \$opt{file_xmatrix},
    "file_ymatrix=s" => \$opt{file_ymatrix},
    "file_cov=s" => \$opt{file_cov},
    "file_ldRef=s" => \$opt{file_ldRef},
    "sumstats=s" => \$opt{sumstats},
    "dir_xmatrix=s" => \$opt{dir_xmatrix},
    "dir_ymatrix=s" => \$opt{dir_ymatrix},
    "dir_cov=s" => \$opt{dir_cov},
    "dir_script=s" => \$opt{dir_script},
    "dir_msg=s" => \$opt{dir_msg},
    
    # "~1=s" => \\$opt{~1},
    # "~1!" => \\$opt{~1},
    # "~1=i" => \\$opt{~1},
    ) or pod2usage(
    verbose => 0,
    exitstatus => 1
    );

if ($opt{help}) {
    pod2usage(
	verbose => 2,
	exitstatus => 0
	);
}

if ($opt{version}) {
    print $Version;
    exit 0;
}

##initiate some opts if not specify
$opt{dir_out}= Cwd::getcwd."/out",if ! $opt{dir_out};
$opt{out}= !$opt{out}?"$opt{dir_out}/out":"$opt{dir_out}/$opt{out}";
$opt{dir_log}="$opt{dir_out}/log",if !$opt{dir_log};

my $bundle=Bundle::Wrapper->new(\%opt);
$bundle->mkdir($opt{dir_out});
$opt{dir_xmatrix}=$opt{dir_out}, if !$opt{dir_xmatrix};
$opt{dir_ymatrix}=$opt{dir_out}, if !$opt{dir_ymatrix};
$opt{dir_cov}=$opt{dir_out}, if !$opt{dir_cov};
$opt{dir_msg}=$opt{dir_out}, if !$opt{dir_msg};
$bundle->mkdir($opt{dir_xmatrix});
$bundle->mkdir($opt{dir_ymatrix});
$bundle->mkdir($opt{dir_cov});
$bundle->mkdir($opt{dir_msg});

if($opt{log}){
    $bundle->mkdir($opt{dir_log});
}


if(!$opt{e}){
    $opt{cmdlist}=$bundle->getSub("main");
    $bundle->opt_print(\%opt);
    print "Add -e to excute\n";
    exit '0';
}


####################################################################
##                               Filehandle
####################################################################
my $fh_in;
open($fh_in,"<$opt{in}");

my $fh_log;
open($fh_log,">",$opt{log}), if $opt{log};

####################################################################
##                               STORAGE
####################################################################
# if($opt{save} && !$opt{load}){
#     store $hash,"$opt{in}.hash";
# }

# if($opt{load} && -e "$opt{in}.hash")
# {
#     $hash=retrieve("$opt{in}.hash");
# }
# elsif($opt{load} && ! -e "$opt{in}.hash")
# {
#     die "need hash.result\n";
# }

# if($opt{load} && -e "$opt{in}.hash")
# {
#     $hash=retrieve("$opt{in}.hash");
# }
# elsif($opt{load} && ! -e "$opt{in}.hash")
# {
#     die "need $opt{in}.hash\n";
# }


####################################################################
##                               MAIN
####################################################################
my %rsid;


#while(<$fh_in>){
#}

foreach(@{$opt{cmd}}){
    switch($_)
    {
	case 'extractSnp' {
	    &extractSnp();
	}
	case 'filterVcfSample' {
	    &filterSample();
	}
	case 'calculateAC' {
	    &calculateAC();
	}
	case 'formatXmatrix' {
	    &formatXmatrix();
	}
	case 'extractSplicing' {
	    &extractSplicing();
	}
	case 'filterSplicingSample' {
	    &filterSplicingSample();
	}
	case 'formatYmatrixForSingleTissue' {
	    &formatYmatrixForSingleTissue();
	}
	case 'cleanIntermediate' {
	    &cleanIntermediate();
	}
	case 'generateXmatrix' {
	    &generateXmatrix();
	}
	case 'generateYmatrix' {
	    &generateYmatrix();
	}
	case 'generateCov' {
	    &generateCov();
	}
	case 'runMSG' {
	    &runMSG();
	}

	else{die "wrong cmd!\n"}
    }
}


####################################################################
##                               SUBS
####################################################################
#sub
sub extractSnp {
    my $cmd;	
    my $bundle=Bundle::Wrapper->new(\%opt);
    $bundle->step_print("extract position of gene");
    $bundle->input("$opt{file_genePos}");
    $bundle->input("$opt{file_snp}");
    
    open(my $fh_gene,"grep $opt{gene} $opt{file_genePos} | ");
    my $geneLine = <$fh_gene>;
    my @gene =split("\t",$geneLine); #"$gene:$chr:$start:$end"

    $opt{file_geneVcf} = $opt{dir_xmatrix}."/".Bundle::File->new($opt{gene})->Filename.".vcf";
    $bundle->output("$opt{file_geneVcf}");
    $cmd = "tabix $opt{file_snp} chr$gene[1]:$gene[2]-$gene[3] -h | grep -v \"^##\" > $opt{file_geneVcf}";
    $bundle->run($cmd);
    &stopWhenLinesLTN($opt{file_geneVcf},2);
}

sub filterVcfSample {
    my $cmd;	
    my $bundle=Bundle::Wrapper->new(\%opt);
    $bundle->step_print("Select samples from vcf");
    $bundle->input("$opt{file_geneVcf}");
    $bundle->input("$opt{file_sample}");
    $opt{file_geneVcfSelect} = $opt{dir_xmatrix}."/".Bundle::File->new($opt{gene})->Filename.".select.rsid.vcf";
    $bundle->output("$opt{file_geneVcfSelect}");
    if($bundle->isrun){
	open(my $fh_file_sample,"<$opt{file_sample}");
	open(my $fh_file_geneVcf,"<$opt{file_geneVcf}");
	chomp(my @sample=<$fh_file_sample>);
	chomp(my @geneVcf=<$fh_file_geneVcf>);
	my @head_geneVcf = split/\t/,$geneVcf[0];
	## print Dumper @head_geneVcf;
	## print Dumper @sample;
	my $head_hash = &arrayIdxInHash(\@head_geneVcf);
	my $select_idx = &getArrayFromHashByName(\@sample,$head_hash);
	## print Dumper $select_idx;
	my $rsid_hash = &rsidToHash();
	## print Dumper $rsid_hash;
	open(my $fh_geneVcfSelect, ">$opt{file_geneVcfSelect}");
	foreach my $geneVcfLine(@geneVcf){
	    my @tmp = split/\t/,$geneVcfLine;
	    if($geneVcfLine=~/^#/){
		print $fh_geneVcfSelect join ("\t",@tmp[(0..8,@$select_idx)])."\n";
	    }
	    if(exists $rsid_hash->{"$tmp[0]:$tmp[1]"}){
		$tmp[2] = $rsid_hash->{"$tmp[0]:$tmp[1]"};
		print $fh_geneVcfSelect join ("\t",@tmp[(0..8,@$select_idx)])."\n";
	    }
	}
    }
    &stopWhenLinesLTN($opt{file_geneVcfSelect},2);
}

sub calculateAC {
    my $cmd;	
    my $bundle=Bundle::Wrapper->new(\%opt);
    $bundle->step_print("calculate allele count");
    $bundle->input("$opt{file_geneVcfSelect}");
    $opt{file_geneVcfSelectAc} = $opt{dir_xmatrix}."/".Bundle::File->new($opt{gene})->Filename.".select.rsid.vcf.ac";
    $bundle->output("$opt{file_geneVcfSelectAc}");
    if($bundle->isrun){
	open(my $fh_geneVcfSelect,$opt{file_geneVcfSelect});
	chomp(my @geneVcfSelect=<$fh_geneVcfSelect>);
	open(my $fh_geneVcfSelectAc, ">$opt{file_geneVcfSelectAc}");
	foreach my $geneVcfSelectLine(@geneVcfSelect){
	    my @F = split/\t/,$geneVcfSelectLine;
	    if($geneVcfSelectLine=~/^#/){
		print $fh_geneVcfSelectAc $geneVcfSelectLine."\n";
	    }elsif(/(\.\|\.)|(\.\/\.)/){
		next
	    }else{
		my @ac = map {&gtToAc($_)} @F[9..$#F];
		my @ac_notna = grep(!/NA/,@ac);
		if(sum(@ac_notna)!=0){
		    ##print $fh_geneVcfSelectAc $geneVcfSelectLine."\n";
		    print $fh_geneVcfSelectAc join ("\t",(@F[0..8],@ac))."\n";
		}
	    }
	}
    }
    &stopWhenLinesLTN($opt{file_geneVcfSelectAc},2);
}

sub formatXmatrix {
    my $cmd;	
    my $bundle=Bundle::Wrapper->new(\%opt);
    $bundle->step_print("Output Xmatrix");
    $bundle->input("$opt{file_geneVcfSelectAc}");
    $opt{file_xmatrix} = $opt{dir_xmatrix}."/".Bundle::File->new($opt{gene})->Filename.".select.rsid.vcf.ac.x_all";
    $bundle->output("$opt{file_xmatrix}");
    if($bundle->isrun){
	open(my $fh_geneVcfSelectAc,$opt{file_geneVcfSelectAc});
	chomp(my @geneVcfSelectAc=<$fh_geneVcfSelectAc>);
	open(my $fh_xmatrix,">$opt{file_xmatrix}");
	foreach my $geneVcfSelectAcLine(@geneVcfSelectAc){
	    my @F = split/\t/,$geneVcfSelectAcLine;
	    if($geneVcfSelectAcLine=~/^#/){
		$F[2]="rsid";
		print $fh_xmatrix join ("\t",(@F[0,2,5,1,3,4,7,9..$#F]))."\n";
		next;
	    }
	    my @info_array=split/;/,$F[7];
	    my %info={};
	    map {my @tmp = split/=/,$_;$info{$tmp[0]}=$tmp[1]} @info_array;
	    $F[7]=$info{"AF"};
	    print $fh_xmatrix join("\t",(@F[0,2,5,1,3,4,7,9..$#F]))."\n";;
	}
    }
    &stopWhenLinesLTN($opt{file_xmatrix},2);
}


sub getSplicing {
    my $cmd;	
    my $bundle=Bundle::Wrapper->new(\%opt);
    $bundle->step_print("Get Splicing from gtex data");
    my @input;
    if($opt{file_splicing}=~/\*/){
	chomp(@input = `ls $opt{file_splicing}`)
    }
    $bundle->input(@input);
    $opt{file_geneSplicing} = $opt{dir_ymatrix}."/".Bundle::File->new($opt{gene})->Filename.".splicing";
    $bundle->output($opt{file_geneSplicing});
    $cmd="parallel -j $opt{threads} --line-buffer -k -q bash -c ' zq --raw '\\''select a.line from index_first a where (a.key like \"ID\" OR a.key like \"\%$opt{gene}\%\")'\\'' {1} | perl -F\"\\t\" -slane '\\''if(\$.==1){print \"tissue\\t\$_\";next}print \$tissue,\"\\t\$_\"'\\'' -- -tissue={=1 s#.*/##g;s#.v\\d+.*##g=}' ::: $opt{file_splicing}  ::: $opt{gene} > $opt{file_geneSplicing}";
    $bundle->run($cmd);
    &stopWhenLinesLTN($opt{file_geneSplicing},3);
}

sub extractSplicing {
    my $cmd;	
    my $bundle=Bundle::Wrapper->new(\%opt);
    $bundle->step_print("extract splicing of gene");
    $bundle->input("$opt{file_genePos}");
    my @input;
    if($opt{file_splicing}=~/\*/){
	chomp(@input = `ls $opt{file_splicing}`)
    }
    $bundle->input(@input);
    
    open(my $fh_gene,"grep $opt{gene} $opt{file_genePos} | ");
    my $geneLine = <$fh_gene>;
    my @gene =split("\t",$geneLine); #"$gene:$chr:$start:$end"

    $opt{file_geneSplicing} = $opt{dir_ymatrix}."/".Bundle::File->new($opt{gene})->Filename.".splicing";
    $bundle->output("$opt{file_geneSplicing}");
    $cmd = "parallel -j $opt{threads} --line-buffer -k -q bash -c 'tabix {1} chr$gene[1]:$gene[4]-$gene[5] -h  | perl -F\"\\t\" -slane '\\''print \$tissue,\"\\t\$_\"'\\'' -- -tissue={=1 s#.*/##g;s#.v\\d+.*##g=}' ::: $opt{file_splicing} > $opt{file_geneSplicing}";
    $bundle->run($cmd);

    $b=Bundle::Wrapper->new(\%opt);
    chomp(my $colnum = `perl -F\"\\t\" -lane 'print scalar(\@F)' $opt{file_geneSplicing} | sort | uniq -c | wc -l`);
    if($colnum!=1){
	$bundle->throw("The gene splicing file has lines with various columns. Make sure input the correct spicing files: a single tissue file or multiple files with the same columns")
    }
    &stopWhenLinesLTN($opt{file_geneSplicing},3);
}


sub filterSplicingSample {
    my $cmd;	
    my $bundle=Bundle::Wrapper->new(\%opt);
    $bundle->step_print("Select samples from geme splicing file");
    $bundle->input("$opt{file_geneSplicing}");
    $bundle->input("$opt{file_sample}");
    $opt{file_geneSplicingSelect} = $opt{dir_ymatrix}."/".Bundle::File->new($opt{gene})->Filename.".select.splicing";
    $bundle->output("$opt{file_geneSplicingSelect}");
    if($bundle->isrun){
	open(my $fh_file_sample,"<$opt{file_sample}");
	open(my $fh_file_geneSplicing,"<$opt{file_geneSplicing}");
	chomp(my @sample=<$fh_file_sample>);
	chomp(my @geneSplicing=<$fh_file_geneSplicing>);
	my @head_geneSplicing = split/\t/,$geneSplicing[0];
	## print Dumper @head_geneSplicing;
	## print Dumper @sample;
	my $head_hash = &arrayIdxInHash(\@head_geneSplicing);
	my $select_idx = &getArrayFromHashByName(\@sample,$head_hash);
	open(my $fh_geneSplicingSelect, ">$opt{file_geneSplicingSelect}");
	foreach my $geneSplicingLine(@geneSplicing){
	    my @tmp = split/\t/,$geneSplicingLine;
	    print $fh_geneSplicingSelect join ("\t",@tmp[(0..4,@$select_idx)])."\n";
	}
    }
    &stopWhenLinesLTN($opt{file_geneSplicingSelect},3);
}


sub formatYmatrixForSingleTissue {
    my $cmd;	
    my $bundle=Bundle::Wrapper->new(\%opt);
    $bundle->step_print("formatYmatrix");
    $bundle->input("$opt{file_geneSplicingSelect}");
    $opt{file_ymatrix} = $opt{dir_ymatrix}."/".Bundle::File->new($opt{gene})->Filename.".all.final.matrix.decomp";
    $bundle->output("$opt{file_ymatrix}");
    $cmd="cat $opt{file_geneSplicingSelect} | cut -f 6- | perl -F\"\\t\" -MList::Util=sum -lane 'if(\$.==1){print join \"\\t\",(\"Expression\",\@F)} if(sum(\@F)==0){next}print \"splice\".(\$.-1).\"\\t\",join \"\\t\",\@F' >  $opt{file_ymatrix}";
    $bundle->run($cmd);
    &stopWhenLinesLTN($opt{file_geneSplicingSelect},3);
}

sub cleanIntermediate {
    my $cmd;	
    my $bundle=Bundle::Wrapper->new(\%opt);
    $bundle->step_print("Clean intermediate files");
    $cmd="rm -rf $opt{file_geneVcf} $opt{file_geneVcfselect} $opt{file_geneVcfSelectAc} $opt{file_geneSplicing} $opt{file_geneSplicingSelect}";
    $bundle->run($cmd);
}

sub generateCov {
    my $cmd;	
    my $bundle=Bundle::Wrapper->new(\%opt);
    $bundle->step_print("generate cov for pair-wise snps for each gene from X matrix");
    chomp(my @input = `ls $opt{file_ldRef}*`);
    $bundle->input(@input);
    $bundle->input("$opt{file_xmatrix}");
    $opt{file_cov} = $opt{dir_cov}."/".Bundle::File->new($opt{file_xmatrix},".x_all")->Prefix.".cov.RData";
    $bundle->output($opt{file_cov});
    $cmd="Rscript $opt{dir_script}/generate_db_and_cov.R --input $opt{file_xmatrix} --dir_out $opt{dir_cov}  --ref_ld_chr $opt{file_ldRef}";
    $bundle->run($cmd);
}

sub runMSG {
    my $cmd;	
    my $bundle=Bundle::Wrapper->new(\%opt);
    $bundle->step_print("run MSG using x matrix, y matrix and cov file");
    $bundle->input($opt{file_xmatrix});
    $bundle->input($opt{file_ymatrix});
    $bundle->input($opt{file_cov});
    $bundle->input($opt{sumstats});
    $opt{rlt} = $opt{dir_msg}."/results/all/".Bundle::File->new($opt{file_xmatrix},".x_all")->Prefix.".results.MSG_GBJ_ACAT.txt";
    $bundle->output($opt{rlt});
    $cmd="Rscript $opt{dir_script}/gtex_comp_MSG_ACAT_GBJ_120522.R --x $opt{file_xmatrix} --y $opt{file_ymatrix} --model_training --save_model --cov $opt{file_cov} --sumstats $opt{sumstats} --dir_out $opt{dir_out} --verbose TRUE";
    $bundle->run($cmd);
}


sub generateXmatrix {
    &extractSnp();
    &filterVcfSample();
    &calculateAC();
    &formatXmatrix();
}

sub generateYmatrix {
    &extractSplicing();
    &filterSplicingSample();
    &formatYmatrixForSingleTissue;
}

####################################################################
##                internal subroutine
####################################################################

sub gtToAc{
    my $gt = $_[0];
    if($gt eq "0|0" || $gt eq "0/0"){
	return 0
    }elsif($gt eq "0|1" || $gt eq "1|0" || $gt eq "0/1"){
	return 1
    }elsif($gt eq "1|1" || $gt eq "1/1"){
	return 2
    }else{
	return "NA"
    }
}

sub arrayIdxInHash{
    my $array = $_[0];
    my $hash;
    for(my $i=0;$i<@$array;$i++){
	$hash->{"$$array[$i]"}=$i;
    }
    return $hash;
}

sub getArrayFromHashByName{
    my $name = $_[0];
    my $hash = $_[1];
    my $array;
    foreach my $key(@$name){
	if(exists $hash->{$key}){
	    push @$array,$hash->{$key};
	}
    }
    return $array;
}

sub rsidToHash{
    my $cmd;	
    my $bundle=Bundle::Wrapper->new(\%opt);
    $bundle->step_print("read rsid to hash...");
    $bundle->input($opt{file_rsid});

    if($bundle->isrun){
	open(my $fh_rsid,$opt{file_rsid});
	chomp(my @rsid = <$fh_rsid>);
	my $hash;
	foreach my $rsidLine(@rsid){
	    my @tmp = split/\t/,$rsidLine;
	    $hash->{"$tmp[0]:$tmp[1]"} = $tmp[2];
	}
	return $hash;
    }else{
	$bundle->throw("filed to read rsid to hash");
    }
}

sub stopWhenLinesLTN {
    my $cmd;	
    my $bundle=Bundle::Wrapper->new(\%opt);

    my $file = $_[0];
    my $threshold = $_[1];
    my $n = `wc -l $file | perl -ane 'chomp(\@F);print \$F[0]'`;
    if($n<$threshold){
	print STDERR "STOP: $file doesn't have enough number of lines (n=$n,threshold=$threshold)...\n";
	exit()
    }
}


####################################################################
##                             TEMP FILE SAVE
####################################################################
my $fh_temp;
if($opt{temp}){
    open($fh_temp,">$opt{out}.temp"),if $opt{out} || open($fh_temp,">out.temp");
}



__END__
####################################################################
##                             Now Docs...
####################################################################
=head1 NAME

$0  - DESCRIBE ME

=head1 SYNOPSIS

$0  [-h] [-v]

=head1 EXAMPLES Step by Step


Generate X matrix

perl ~/script/msg.pl -cmd generateXmatrix -gene ENSG00000000457 -file_genePos /nobackup/cgg/yany14/MSG/data/gene.500k.id -file_snp /nobackup/cgg/yany14/MSG/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz  -file_sample /nobackup/cgg/yany14/MSG/data/samples.filtered.used -e


Generate Y matrix

perl ~/script/msg.pl -cmd generateYmatrix -gene ENSG00000000457 -file_genePos /nobackup/cgg/yany14/MSG/data/gene.500k.id -file_splicing /fs0/yany14/MSG/data/rawSplice/Whole_Blood.v8.leafcutter_phenotypes.bed.gz -file_sample /nobackup/cgg/yany14/MSG/data/samples.filtered.used -e


Generate cov

perl ~/script/msg.pl -cmd generateCov -file_xmatrix /fs0/chenr6/other/yany/MSG/out/ENSG00000000419.select.rsid.vcf.ac.x_all --file_ldRef /nobackup/cgg/yany14/MSG/data/LDREF/1000G.EUR. -e

Run MSG

perl ~/script/msg.pl -cmd runMSG -file_xmatrix /fs0/chenr6/other/yany/MSG/out/ENSG00000000419.select.rsid.vcf.ac.x_all --file_ymatrix /fs0/chenr6/other/yany/MSG/out/ENSG00000000419.all.final.matrix.decomp --cov /fs0/chenr6/other/yany/MSG/out/ENSG00000000419.select.rsid.vcf.ac.cov.RData -sumstats /fs0/yany14/MSG/data/sumstats/clozukscz.sumstats


=head1 EXAMPLES In one step

perl ~/script/msg.pl -cmd generateXmatrix -cmd generateYmatrix -cmd generateCov -cmd runMSG -cmd cleanIntermediate -gene ENSG00000000457 file_genePos data/gene.500k.id -file_snp data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz  -file_sample data/samples.filtered.used -file_splicing data/rawSplice/Whole_Blood.v8.leafcutter_phenotypes.bed.gz --file_ldRef data/LDREF/1000G.EUR. -sumstats data/sumstats/clozukscz.sumstats -dir_script scripts/ -e


=head1 OPTIONS

=over 1
		
=item B<-h|--help>

Print help message and exit successfully.

=item B<-v|--version>

Print version information and exit successfully.

=item B<-excute|e>
Run cmd

=item B<-dir_out>

Folder of output.
Default: pwd/out


=item B<-cmd>

Command build in this script, including:

=over 1

=item B<extractSnp>

Purpose : Extract SNPs around gene.

Input: --gene, --file_genePos, --file_snp

Output: file_geneVcf


=item B<filterVcfSample>

Purpose: Extract only samples in a sample file.

Input: --file_geneVcf, --file_sample. --file_geneVcf can be specified or generated from "--cmd extractSnp"

Output: file_geneVcfSelect


=item B<calculateAC>

Purpose: Calculate Allele Count for each indivudial.

Input: --file_geneVcfSelect

Output: file_geneVcfSelectAc


=item B<formatXmatrix>

Purpose: generate x matrix

Input: --file_geneVcfSelectAc

Output: --file_xmatrix


=item B<extractSplicing>

Purpose: extract splicing for gene 

Input: -gene, --file_genePos, --file_splicing

Output: file_geneSplicing


=item B<filterSplicingSample>

Purpose: select samples from gene splicing file	

Input: --file_geneSplicing

Output: file_geneSplicingSelect


=item B<formatYmatrixForSingleTissue>

Purpose: format Y matrix

Input: --file_geneSplicingSelect

Output: file_ymatrix


=item B<generateXmatrix>

Purpose: generate X matrix in one step

Input:  --gene, --file_genePos, --file_snp, --file_sample

Output: file_xmatrix


=item B<generateYmatrix>

Purpose: generate Y matrix in one step. Including "-cmd extractSplicing -cmd filterSplicingSample, -cmd formatYmatrixForSingleTissue";

Input: --gene, --file_genePos, --file_splicing, --file_sample

Output: file_ymatrix


=item B<generateCov>

Purpose: generate cov for pair-wise snps for each gene from X matrix

Input: --file_ldRef, --file_xmatrix

Output: file_cov


=item B<runMSG>

Purpose: run MSG using x matrix, y matrix and cov file

Input: --file_xmatrix, --file_ymatrix, --file_cov, --sumstats

Output: rlt


=item B<cleanIntermediate>

Purpose: Clean intermediate files

=back

=back


=head1 INPUT OPTIONS

=item B<--gene>

Which gene to run

=item B<--file_genePos>

File including the upstreaming and downstreaming 500k for each gene

=item B<--file_snp>

File including the genotypes

=item B<--file_rsid>

File including snp and rsid mapping

=item B<--file_sample>

File including samples that should be used 

=item B<--file_splicing>

File including splicing expression. 

=item B<--file_ldRef>

File including the LD info for snps


=item B<--sumstats>

File including summary statistics


=head1 DIRECTORY OPTIONS

=item B<--dir_script>

Dir including all scripts

=item B<--dir_xmatrix>

Dir including x matrix. Default is --dir_out.

=item B<--dir_ymatrix>

Dir including y matrix. Default is --dir_out.

=item B<--dir_cov>

Dir including cov files. Default is --dir_out.

=item B<--dir_msg>

Dir including final results. Default is --dir_out.


=item B<-header>

Parse header 

=item B<-out>

Prefix of output files. Default: out

=item B<-t|threads>
      
Threads used 
      
=item B<-save>
      
Save hash if there exists.

=item B<-load>

Load hash if there exists.

=item B<chkopt>

Check auto-parsed opts

=item B<-version>

Version of this script.
		
=back

=cut
