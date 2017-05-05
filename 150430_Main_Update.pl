#!/usr/bin/perl 
use strict;
use warnings;

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

#Rewriting

##opens a text file to write any errors to "problems.txt" and associates it with the file handler "PROBLEM"
open (PROBLEM, ">>problems.txt") || die "Can not open file for messages about problems during update.\n";

##Set change flag to one to force the initial cycle through the process of updating
##the relevant gi numbers.
my $change_flag = 1;

##create a new Common functions object to have access to its methods
my $cmObj = new common_functions();

#get admin email
my $adminEmail;
$cmObj->getEmail(\$adminEmail);


##make a time variable via calling the get current time method from the common funtions class
my $newDate;
$newDate = $cmObj->get_current_time();


##create a new DBConnection object
my $myConnObj = new DB_connection();

##call the open db connection method
my $dbh = $myConnObj->open_db_connection();

##create a new db queries object
my $myDBqObj = new DB_Queries();

##create a variable new_refindex and call the method to Get the last ref_id value from VKCDB, since new reference records may be created
my $new_refindex;
$myDBqObj->ref_id(\$dbh, \$new_refindex);


#declare some global variables for holding data, an array of gi numbers and a set of hashes keyed on the gi numbers for other parameters
my @ginumber_search = ();
my %vkcid;
my %lastmod;
my %code_taxon;
my %vkcnt_id;

my $count;

#call the vkcID_gi_lastMod method to set the passed in variables
$myDBqObj->vkcID_gi_lastMod(\$dbh, \@ginumber_search, \%vkcid, \%lastmod, \%vkcnt_id, \%code_taxon, \$count);

#Block of alternative code for single inputs for testing - Comment out when single runs are working.
#print "Input gi number for single entry test:\n";

#$ginumber_search[0]=<STDIN>;
#chomp $ginumber_search[0];

#print "Input VKCID number for entry:\n";

#$vkcid{$ginumber_search[0]}=<STDIN>;
#chomp $vkcid{$ginumber_search[0]};

#print "Input last date modifed (yyyy-m-d):\n";

#$lastmod{$ginumber_search[0]}=<STDIN>;
#chomp $lastmod{$ginumber_search[0]};

#print "Input VKCNT number:\n";

#$vkcnt_id{$ginumber_search[0]}=<STDIN>;
#chomp $vkcnt_id{$ginumber_search[0]};

#print "Input taxon code\n";

#$code_taxon{$ginumber_search[0]}=<STDIN>;
#chomp $code_taxon{$ginumber_search[0]};

#$count =1;

#End of test block
#test line
print "$count records for evaluation.\n";

#declare some more  global variables for holding data
my $retry;
my $file = 'temp_hold.gb';
my $gpeptfactory;
my $history;
my $history2;
my $retr_seq;
my $retmax;
my $retstart;
my $temp_file = 'buffer_file';
my $rec_rem = 0;

	my $ncbiObj = new NCBI_Retrival();

	$ncbiObj->eutilities_getData(\$file, \@ginumber_search, \%vkcid, \%lastmod, \%vkcnt_id, \$adminEmail);

#Open the output file for the SQL commands that will update VKCDB

open( OUT, '>','revise_VKCDB.sql' )
  || die "Can't create the SQL output file.\n";


#last of the global variable used mostly by the sql_output subroutine global variables
my %protein_gi_w_prot_acc;
my %aaseq;
my @pro_acc;
my %flag;
my $pubmed_list = "";
my %pubmed_hash;
my $GenBank_flag  = 0;
my $REFSEQ_flag  = 0;
my $EMBL_flag  = 0;
my @genbank_na_acc = ();
my @refseq_na_acc = ();
my @embl_na_acc = ();
my $j = 0;
my $k = 0;
my $l = 0;
my $start_frame;
my %codon_table;	
my $na_eutil;
my $genbank_count= 0;
my $file3 = 'na_gi_numbers.txt';

# This file was created by the retrieval subroutine

my $file4 = 'removed_gis.txt';

open (REMOVED, $file4)||die "Can not open file of dead entries, $file4.\n";

while (<REMOVED>){
	if (($_=~/(\d+), Status, dead/)||($_=~/(\d+), Status, suppressed/)){
	
		print OUT "UPDATE protein SET flag = \"N\" WHERE ginumber = \"$1\";\n";
	
	}
}

#This reopens the file that holds all of the GENBANK entries and starts the process of
#creating the necessary SQL update statements

$retr_seq = Bio::SeqIO->new(
    -file   => $file,
    -format => 'genbank'
);

#Cycle through the retrieved records
while ( my $record = $retr_seq->next_seq ) {
#Hide this line when 500 sample test run is desired


#Testblock limiting the search to the first 500, for testing only

#for(my $cycle=0; $cycle<500; $cycle++){
#	my $record = $retr_seq->next_seq;
#Hide these two lines  and matching bracket at bottom when full run is desired

#TEST LINE
print "cycling through records\n";
	#Start of comparing the dates with the previously obtained LastMod dates from VKCDB
	    
    #Get the sequence
    my $gproseq = $record->seq();

    #		print "$gproseq\n";
    my $ginumber = $record->primary_id();
    my $pro_acc  = $record->accession_number();
    $pro_acc =~ s/\.\d+//;
    $protein_gi_w_prot_acc{$pro_acc} = $ginumber;
    $aaseq{$protein_gi_w_prot_acc{$pro_acc}} = $gproseq;

    #    	print "$ginumber\n";
    my $anno_collection = $record->annotation;

#This goes through the date annotations in the record, finding the most recent date and saving it
#in the form yyyymmdd for numerical comparison with the last_mod date from the VKCDB entry

    my $startdate = 0;
    my @dates     = $record->get_dates;
    my $comp_newdate;
    foreach my $nextdate (@dates) {

        #	   		print "Next date is $nextdate\n" ;
        $nextdate =~ /(\d{2})-(\w{3})-(\d{4})/;
        
        #call the month_converte method from the common functions class
        my $num_month = $cmObj->month_converter($2);
        $comp_newdate = "$3" . $num_month . "$1";

        #    		print "Reformatted for comparison is $comp_newdate\n";
        if ( $comp_newdate > $startdate ) {
            $startdate = $comp_newdate;
        }
    }

    my $acctemp      = $ginumber;
    my $comp_olddate = $lastmod{$acctemp};

    #    	print "Date for comparison is $comp_olddate\n";
    $comp_olddate =~ s/-//g;

	#End of date comparision part of the while loop
#___________________________________________________________________________________
#Beginning of comparision if dates are different here then the records is parsed and a sql update
#statemt is made for it.

    if ( $comp_newdate > $comp_olddate ) {
#TEST LINE
print "date match\n";
        #Identify the corresponding vkc_id value

        my $vkcid = $vkcid{$acctemp};

    # return some annotations to "", so they don't add up for different entries.

        my $gpfam     = "";
        my $gsmart    = "";
        my $glocus    = "";
        my $gmim      = "";
        my $gmgd      = "";
        my $gprodef   = "";
        my $gAcc      = "";
        my $gorganism = "";
        my $gntAcc    = "";
        my $genNT_vers;
        my $cds_range;
        my $gVers;
        my $gCDD = "";

        my $gtaxon = "";
        my ( @gtitle, @gauthor, @grefpos, @gjournal, @gmedline, @gpubmed );

        #parse definition from GenBank.
        $gprodef = $record->description;

        #parse the primary gAcc.
        $gAcc = $record->accession_number;

        #parse the version number

        $gVers = $record->seq_version;




#this section needs to be in an if block to test if the $record is initialized by species->species call
#if (exists take a hash and its key, probablu use
#$record->species->species; if its defined do the following line else send a warning. )
#$record is a sequence object, 

        #parse organism.
        $gorganism = $record->species->binomial();
        #TEST LINE
        print "$gorganism\n";
		$code_taxon{$ginumber} = $gorganism;
        #parse references
        #cycle through references, collecting relevant information
        my @references = $anno_collection->get_Annotations('reference');

        foreach my $temp_ref (@references) {

 #If the PUBMED reference number already exists, then write an SQL
 #statement that just adds a link between the protein entry and the reference
 #entry in the protein_reference table.  Otherwise parse the information for
 #the new reference, assign it the next ref_id value and write a new SQL command
 #to add a new line to the reference table and write a new SQL command to add
 #a new link to the protein_reference table.

         #recover the hash reference that contains all the reference information
            my $hash_ref = $temp_ref->hash_tree;

            #parse reference title.
            my $title = $hash_ref->{'title'};
            if ( defined($title) ) {
                $title =~ s/"//g;
            }

            my $authors = $hash_ref->{'authors'};
            if ( defined($authors) ) {
                $authors =~ s/"//g;
            }

            #			print "$authors\n";
            my $biblio = $hash_ref->{'location'};
            $biblio =~ s/"//g;

            #			print "$biblio\n";
            #			my $medline = $hash_ref->{'medline'};
            my $pubmed = $hash_ref->{'pubmed'};

            #			print "$pubmed\n";

           #Do not enter a reference if it is a Direct Submission or Unpublished

            if ( !defined($title) ) {
                $title = 'ERRATUM';
                print "Title ($title) not defined for an entry in $ginumber\n";
            }
            if ( !defined($biblio) ) {
                print
"No bibliographic information defined for an entry in $ginumber\n";
                $biblio = 'No bibliographic information available';
            }
            if (   ( $title eq "Direct Submission" )
                || ( $title eq "" )
                || ( $biblio =~ /Submitted.*EMBL/ )
                || ( $biblio =~ /Unpublished/ ) )
            {
                next;
            }
            if ( !defined($pubmed) ) {
                print "No PUBMEDID for entry in $ginumber\n";
                next;
            }

            my $sth2 = $dbh->prepare(
                "SELECT COUNT(*) FROM reference
                         WHERE pubmed_id = \"$pubmed\""
            ) || die "Can't prepare statement:$dbh->errstr\n";

            my $rv2 = $sth2->execute
              || die "Can't execute the query: $sth2->errstr\n";

            my $rv22 = $sth2->fetchrow_array;

            if ( $rv22 == 0 ) {
                if ( $pubmed_list =~ /$pubmed/ ) {
                    print "PUBMEDID $pubmed already found earlier in update\n";

                    print OUT
"INSERT INTO protein_reference VALUES (\"VKC$vkcid\", \"Ref$pubmed_hash{$pubmed}\");\n";

                }
                else {

                    print OUT
"INSERT INTO reference VALUES (\"Ref$new_refindex\", \"$authors\", \"$title\", \"$biblio\", \"\", \"\", \"\", \"$pubmed\");\n";

                    print OUT
"INSERT INTO protein_reference VALUES (\"VKC$vkcid\", \"Ref$new_refindex\");\n";

                    $pubmed_hash{$pubmed} = $new_refindex;

                    $pubmed_list = $pubmed_list . ", " . $pubmed;

                    $new_refindex++;
                }
            }
            elsif ( $rv22 == 1 ) {
                my $sth3 =
                  $dbh->prepare(
                    "SELECT  ref_id FROM reference WHERE pubmed_id = $pubmed")
                  || die
"Could not prepare the statement for retrieving Ref ID given PubMed value.\n";

                my $rh        = $sth3->execute;
                my @temparray = $sth3->fetchrow_array;
                my $refid     = $temparray[0];

#If there is no matched entry for the current VKCID and Reference number, then insert a new entry into
#the protein_reference table.
                my $sth4 =
                  $dbh->prepare(
"SELECT vkc_id, ref_id FROM protein_reference WHERE vkc_id = \"VKC$vkcid\" AND ref_id = \"$refid\""
                  )
                  || die
                  "Could not prepare statement for retrieving VKC RefID pair\n";
                my $rv4 = $sth4->execute
                  || die
"Could not execute statement for retrieving VKC RefID pair.\n";
                my @temp_hold = $sth4->fetchrow_array;
                if ( !( defined( $temp_hold[0] ) ) ) {
                    print OUT
"INSERT INTO protein_reference VALUES (\"VKC$vkcid\", \"$refid\");\n";
                }
            }
            else {
                print
"There is more than one entry for PUBMED = $pubmed, requires manual fix.\n";
            }

        }#end of going throug heach of the references info for the record

        #parse crossreferenced database information.
        #Start with finding ORF nucleic acid sequence
        #Set default gennnt as unannotated, changed if CDS is found

        $gntAcc     = "Unannotated";
        $genNT_vers = "";
        $cds_range  = "";

        foreach my $seq_feature ( $record->get_SeqFeatures ) {

            if ( $seq_feature->primary_tag eq "CDS" ) {

                if ( $seq_feature->has_tag('coded_by') ) {

                    my @cds_annot = $seq_feature->get_tag_values('coded_by');
                    foreach my $cds_annot (@cds_annot) {

                        if ( $cds_annot =~
                            /^([0-9A-Za-z_]+)\.(\d*)\:<?.*?(\d+\.\..*\d+)/ )
                        {
                            $gntAcc     = $1;
                            $genNT_vers = $2;
                            $cds_range  = $3;
                            $genbank_na_acc[$genbank_count] = $gntAcc;
							$genbank_count++;
                        }
                        elsif ( $cds_annot =~
/^complement\(join\(([0-9A-Za-z_]+)\.(\d*)\:.*(\d+\.\..*\d+)\)\)/
                          )
                        {
                            $gntAcc     = $1;
                            $genNT_vers = $2;
                            $cds_range  = $cds_annot;
                            $genbank_na_acc[$genbank_count] = $gntAcc;
							$genbank_count++;
                        }
                        elsif ( $cds_annot =~
                            /^join\(([0-9A-Za-z_]+)\.(\d*)\:.*(\d+\.\..*\d+)\)/
                          )
                        {
                            $gntAcc     = $1;
                            $genNT_vers = $2;
                            $cds_range  = $cds_annot;
                            $genbank_na_acc[$genbank_count] = $gntAcc;
							$genbank_count++;

                        }
                        else {
                            $gntAcc     = "Unannotated";
                            $genNT_vers = "";
                            $cds_range  = "";
                        }
                    }
                }
                print "Protein $ginumber_search[0] has CDS in $gntAcc\n";
            }

            if ( $seq_feature->has_tag('db_xref') ) {

                foreach my $temp_tag ( $seq_feature->get_tag_values('db_xref') )
                {

                    #		print "$temp_tag\n";

                    #pares taxon id number

                    if ( $temp_tag =~ /taxon\:(\d+)/ ) {

                        $gtaxon = $1;
                    }

                    #parse pfam.

                    if ( $temp_tag =~ /CDD\:(pfam\d*)/ ) {
                        $gpfam .= " $1";
                        $gpfam =~ s/^\s//;    #Take the front space off.
                    }

                    #parse smart.

                    if ( $temp_tag =~ /CDD\:(smart\d*)/ ) {
                        $gsmart .= " $1";
                        $gsmart =~ s/^\s//;
                    }

                    #parse straight CDD references

                    if ( $temp_tag =~ /CDD\:(\d*)/ ) {
                        $gCDD .= " $1";
                        $gCDD =~ s/^\s//;
                    }

                    #parse Locus ID.

                    if ( $temp_tag =~ /LocusID\:(\d*)/ ) {
                        $glocus .= " $1";
                        $glocus =~ s/^\s//;
                    }

                    #parse OMIM ID.

                    if ( $temp_tag =~ /MIM\:(\d*)/ ) {
                        $gmim .= " $1";
                        $gmim =~ s/^\s//;
                    }

                    #parse MGD ID.

                    if ( $temp_tag =~ /MGD\:(\d*)/ ) {
                        $gmgd .= " $1";
                        $gmgd =~ s/^\s//;
                    }

                }

            }
        }#end of going through all the seq features objects for the record

        #Create the update statements for VKCDB and print to the OUT file

        #update protein table.
        if ( !defined($gAcc) ) {
            print "No accession found for gi=$ginumber\n";
        }
        if ( !defined($gVers) ) {
            print "No accession version found for gi=$ginumber\n";
            $gVers = 0;
        }
        if ( !defined($gprodef) ) {
            print "No protein definition found for gi=$ginumber\n";
        }
        if ( !defined($gproseq) ) {
            print "No protein sequence found for gi=$ginumber\n";
        }
        if ( !defined($newDate) ) {
            print "No new modification date found for gi=$ginumber\n";
        }

        #			$gprodef =~ s/\[|\]|\.|\-/ /g;
        print OUT
"UPDATE protein SET ginumber = \"$ginumber\", pro_acc = \"$gAcc\", pro_acc_ver = \"$gVers\", pro_def = \"$gprodef\", aaseq = \"$gproseq\", lastmod = \"$newDate\" WHERE vkc_id = \"VKC$vkcid\";\n";

        #update DBref table.

        print OUT
"UPDATE dbref SET mim_id = \"$gmim\", locus_id = \"$glocus\", mgd_id = \"$gmgd\", pfam_id = \"$gpfam\", smart_id = \"$gsmart\", cdd_ids = \"$gCDD\" WHERE dbref_id = \"DBref$vkcid\";\n";

#update GenNT table.
#-----------------------------------------------------------------------------------------------------------------------------------
#090511 - Inserted and modified code from 090121_Update_Empty_CDS.pl to do full ORF update within this script rather than having todo it as a separate operation
#140530 - Removed dblinks annotation recovery - was pointing the nucleic acid accession search to wrong sources.        
    }#end of if statement of checking the date 

}#end of while loop going through each record in temp_hold.gb


#First use Eutilities to retrieve the gi numbers corresponding to the accession numbers that have
#been collected above.

print "\n";
print "Starting to recover the sequence records based on the accession numbers.\n";


my $number = scalar(@genbank_na_acc);
print "Initial number of records of array before sub dividing ".$number."\n";

my $min = 0;
my $max = 5;
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
	print "last post index positions $min    ".($number-1);
	}
	print "subarray had this many elements ".scalar(@subarray)."\n";
	for (my $i=0; $i<@subarray; $i++){
		$subarray[$i].="[accn]";
	}
	
	$retry = 0;
	my $factory;
	my $record_count;
	my $term_search = join(' OR ',@subarray);
	print "term search string is $term_search.\n";
	REDO_ESEARCH:
	eval{
	$factory = Bio::DB::EUtilities-> new(-eutil => 'esearch',
										-email => $adminEmail,
										-db => 'nucleotide',
										-term => $term_search,
										-usehistory => 'y',
										);
	};
	if ($@) {
	die "Server error on esearch: $@.  Try again later" if $retry == 5;
    print STDERR "$@\n";
    print STDERR "Server error on esearch and fetch, redo #$retry\n";
    $retry++;
    sleep(5);
    goto REDO_ESEARCH;
	}

	REDO_ECOUNT:
	eval{
	$record_count = $factory -> get_count;
	};
		if ($@) {
	die "Server error on esearch: $@.  Try again later" if $retry == 5;
    print STDERR "$@\n";
    print STDERR "Server error on esearch and fetch, redo #$retry\n";
    $retry++;
    sleep(5);
    goto REDO_ECOUNT;
	}

	
	print "esearch only posted this many to the history server ".$record_count."\n";
	
	my $hist = $factory-> next_History || die 'No History Data returned';
	print "History Returned\n";
	REDO_FETCH:
	eval{
	$factory -> set_parameters(-eutil => 'efetch',
						-rettype => 'gbwithparts',
						-retmode => 'text',
						-file => $file3,
						-history => $hist
						);
	};
	if ($@) {
	die "Server error on fetch: $@.  Try again later" if $retry == 5;
    print STDERR "$@\n";
    print STDERR "Server error on esearch and fetch, redo #$retry\n";
    $retry++;
    sleep(5);
    goto REDO_FETCH;
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
			
	print "Retrieved $retstart1.\n";
	$retstart1 += $retmax1;
	}
}
close WRITE;

print "Finished Eutilities nucleic acid record recovery.\n";


#--------------------------------section it off----------------------------------------
my $target_prot_acc;
my $temp_acc;
my $temp_acc_vers;
my $test_prot_obj;

#Go through each of the retrieved records and find the CDS feature corresponding to the protein
#First create a filehandle to the Entrez Taxonomy DB; need to query for correct codon table info

my $entrez_dbh = Bio::DB::Taxonomy::entrez->new();

my $retr_seq_na = Bio::SeqIO -> new(-file   => $file3,-format => 'genbank');



while (my $temp_na_record = $retr_seq_na->next_seq) {

    my $rec_acc = $temp_na_record->accession_number;
    print "Retrieved na record w/ acc no $rec_acc.\n";
    for my $feat_object ( $temp_na_record->get_SeqFeatures ) {
        if ( $feat_object->primary_tag eq "CDS" ) {

            #	   		print "Found a CDS in $rec_acc.\n";
            if ( $feat_object->has_tag('protein_id') ) {
                my @target_prot_acc = $feat_object->get_tag_values('protein_id');
                foreach $target_prot_acc (@target_prot_acc) {
                #Get rid of the version number, it isn't on the protein accession number value
                	$target_prot_acc =~ s/\.\d+//;
                	#print "Protein w/ Acc No $target_prot_acc found.\n";
                    if ( defined( $protein_gi_w_prot_acc{$target_prot_acc} ) ) {

                        #Need to find the start position for translation
                        if ( $feat_object->has_tag('codon_start') ) {
                            my @start_posn =
                              $feat_object->get_tag_values('codon_start');
                            foreach my $place (@start_posn) {
                                $start_frame = $place - 1;
                            }
                        }
                        print
"CDS encodes protein $target_prot_acc, gi number is $protein_gi_w_prot_acc{$target_prot_acc}, frame is $start_frame.\n";

#If the protein accession number of the CDS feature matches one from the records being updated
#then get the sequence from the CDS feature, translate and compare to the matched aa record
                        my $db_prot =
                          $aaseq{ $protein_gi_w_prot_acc{$target_prot_acc} };
                        my $temp_na_seq = $feat_object->spliced_seq->seq;
                        my $temp_cds_obj = Bio::Seq->new( -seq => $temp_na_seq );
                        my $taxon_info = $entrez_dbh->get_taxon(-gi => $protein_gi_w_prot_acc{$target_prot_acc}, -db => 'protein'  );
                        my $code_table_num = $taxon_info->genetic_code();
                        print
"Translating with genetic code $code_table_num from taxon $code_taxon{$protein_gi_w_prot_acc{$target_prot_acc}}\n";
                        my $test_prot_obj = $temp_cds_obj->translate(
                            -codontable_id => $code_table_num,
                            -complete      => 1,
                            -frame         => $start_frame
                        );
                        my $trans_prot = $test_prot_obj->seq;

                        #Trim off any translation of the stop codon to *
                        $trans_prot =~ s/\*$//;
                        if($db_prot eq ""){
                        	$db_prot = $trans_prot;
                        	 $aaseq{ $protein_gi_w_prot_acc{$target_prot_acc} } = $trans_prot;
                        	 print OUT "UPDATE protein SET aaseq =\"$db_prot\" WHERE pro_acc = \"$target_prot_acc\";\n";
                        }

#Make sure that the first amino acid is M, this avoids rejecting bacterial CDSs with alternative start codons
                        if ( !( $db_prot eq $trans_prot ) ) {
                            print PROBLEM
"Sequence mismatch between gi $protein_gi_w_prot_acc{$target_prot_acc} and $rec_acc.\n";
                            print PROBLEM "$db_prot\n\n$trans_prot\n\n";
                            next;
                        }
                        my $temp_range = "";
                        my $temp_gene_name;
                        my @temp_gene_name;
                        if ( $trans_prot =
                            $aaseq{ $protein_gi_w_prot_acc{$target_prot_acc} } )
                        {
                            if ( $rec_acc =~ m/(\w+)\.(\w+)/ ) {
                                $temp_acc      = $1;
                                $temp_acc_vers = $2;
                            }
                            else {
                                $temp_acc      = $rec_acc;
                                $temp_acc_vers = 0;
                            }
                            if (
                                $feat_object->location->isa('Bio::Location::SplitLocationI')
                              )
                            {
                                foreach my $location (
                                    $feat_object->location->sub_Location )
                                {
                                    $temp_range =
                                        $temp_range
                                      . $location->start . ".."
                                      . $location->end . ", ";
                                }
                            }
                            elsif ( defined( $feat_object->location ) ) {
                                $temp_range =
                                    $feat_object->location->start . ".."
                                  . $feat_object->location->end;
                            }
                            else {
                                $temp_range = "";
                            }
                            if ( !( $feat_object->has_tag('gene') ) ) {
                                $temp_gene_name = "";
                            }
                            elsif (
                                !(
                                    @temp_gene_name =
                                    $feat_object->get_tag_values('gene')
                                )
                              )
                            {
                                $temp_gene_name = "";
                            }
                            else {
                                $temp_gene_name = $temp_gene_name[0];

                            }
                            print OUT
"UPDATE gennt SET accnum = \"$temp_acc\", acc_vers = \"$temp_acc_vers\", code_range = \"$temp_range\", gene_name = \"$temp_gene_name\", cds = \"$temp_na_seq\" WHERE vkcnt_id = \"$vkcnt_id{$protein_gi_w_prot_acc{$target_prot_acc}}\";\n";
                        }
                    }
                }
            }
        }
    }
}
$retr_seq->close;


  
#----------------finish code call PorE----------------------
#create a sql output object
my $sqlObj = new SQL_output();

#call the PorE method in the sql output class
$sqlObj-> P_or_E(\$dbh);

my $rc2 = $dbh->disconnect
  || die "Could not close the DB handle at end of script run.\n";

close OUT;
close PROBLEM;
exit;

