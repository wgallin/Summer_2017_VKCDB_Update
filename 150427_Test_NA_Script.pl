#!/usr/bin/perl 
use strict;
use warnings;
#140527 script to test for failures of the nucleic acid record recovery in Main_Update.pl
use lib "/System/Library/Perl/Extras/5.12/XML/:";
#use lib "/opt/local/lib/perl5/site_perl/5.12.3/:";

use LWP::Simple;
use Bio::DB::GenPept;
use Bio::DB::Taxonomy;
use Bio::Seq;
use Bio::Seq::SeqBuilder;
use Bio::DB::EUtilities;
use Bio::DB::Taxonomy::entrez;

require common_functions; 
require DB_connection;
require DB_Queries;
require NCBI_Retrival;
require SQL_output;


print "\n";
print "Starting to recover the sequence records based on the accession numbers.\n";

#Block of code to replace the upstream material from the main scrript

my @genbank_na_acc = ("NM_010595", "X17621", "NM_004770");
my $file3 = "NA_TEST_150427.txt";
my $adminEmail = "wgallin\@ualberta.ca";

my $number = scalar(@genbank_na_acc);
print "Initial number of records of array before sub dividing ".$number."\n";

my $min = 0;
my $max = 20;
my $step = $max - $min;

my @subarray;

open (WRITE, '>', $file3) || die "Can't open file!\n";

while ($min < $number){

	if ($max < $number){ 
	print "start index positions: $min     $max\n";
	@subarray = @genbank_na_acc[$min...$max];
	}
	else{
	@subarray = @genbank_na_acc[$min...$number-1];
	print "last post index positions $min    ".($number-1)."\n";
	}
	print "subarray had this many elements ".scalar(@subarray)."\n";
	
	my $retry = 0;
	my $factory;
	my $record_count;
	
	REDO_ESEARCH_FETCH:
	eval{
	$factory = Bio::DB::EUtilities-> new(-eutil => 'esearch',
										-email => $adminEmail,
										-db => 'nucleotide',
										-term => join(',', @subarray),
										-field => '[ACC]',
										-usehistory => 'y',
										);



	$record_count = $factory -> get_count;
	print "esearch only posted this many to the history server ".$record_count."\n";
	
	my $hist = $factory-> next_History || die 'No History Data returned';
	print "History Returned\n";


	$factory -> set_parameters(-eutil => 'efetch',
						-rettype => 'gbwithparts',
						-retmode => 'text',
						-history => $hist
						);
	};
	if ($@) {
	die "Server error on esearch and fetch: $@.  Try again later" if $retry == 5;
    print STDERR "$@\n";
    print STDERR "Server error on esearch and fetch, redo #$retry\n";
    $retry++;
    sleep(5);
    goto REDO_ESEARCH_FETCH;
	}

	$min = ($max +1);
	$max += $step;

	my $retry1 = 0;						
	my $retmax1 = 500;
	my $retstart1 = 0;



	RETRIEVE_SEQS:
	while ($retstart1 < $record_count){
	$factory -> set_parameters(-retmax => $retmax1,
								-retstart => $retstart1,
								);
		eval{
		$factory->get_Response(-file => sub {my($data)=@_; print WRITE $data});
		};
		if ($@){
		die "Server error at nucleic accid seq recovery: $@. Try again later" if $retry1 == 5;
		print STDERR "Server error at nucleic accid seq recovery\n$@\n, redo #$retry1\n";
		$retry1++ && redo RETRIEVE_SEQS;
		}
			
	print "Retrieved on attempt number $retstart1.\n";
	$retstart1 += $retmax1;
	}
}
close WRITE;

print "Finished Eutilities nucleic acid record recovery.\n";

