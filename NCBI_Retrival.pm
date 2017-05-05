#!/usr/bin/perl 
package NCBI_Retrival;
use lib "/opt/local/lib/perl5/site_perl/5.12.3/";

use strict;
use warnings;
use LWP::Simple;
use DBI;
use Bio::DB::GenPept;
use Bio::DB::Taxonomy;
use Bio::Seq;
use Bio::DB::EUtilities;
use Bio::DB::Taxonomy::entrez;
use Bio::SeqIO;
use Bio::Seq::SeqBuilder;
use Bio::DB::SeqVersion;


sub new {
    my $class = shift;
    my $self = {};
    
    bless $self, $class;
    return $self;
}

##################################################################################
##################################################################################
#Create new EUTILS  history object for submitting a large array of gi numbers and
#attempt to retrieve it - try 5 times and if it fails then quit - there is either
#something wrong with the server or the script is not working correctly - check for whether
#updates to eUtilities at NCBI or to the Eutilities module have changed the behavior if the
#server is functioning as it should.
##################################################################################
##################################################################################
sub eutilities_post {
	my $this = shift;
	my $retry = shift;
	my $gpeptfactory = shift;
	my $ginumber_search = shift;
	my $count = shift;
#	my $count = 500;
	my $history = shift;
	my $file = shift;
	my $adminEmail = shift;

  eval {
        $$gpeptfactory = Bio::DB::EUtilities->new(
            -eutil          => 'epost',
            -db             => 'protein',
            -rettype        => 'gp',
            -retmode        => 'text',
            -tool           => 'VKCDB_Update',
            -email          => $$adminEmail,
            -id             => $ginumber_search,
            -keep_histories => 1
        );
    };
    if ($@) {
        die "Server error on post: $@.  Try again later, post not working" if $$retry == 5;
        print STDERR "$@\n";
        print STDERR "Server error on post, redo #$$retry\n";
        $$retry++;
        sleep(5);
        goto RETRIEVE_HIST;
    }
    else {
        print "Posted successfully\n"; 
    }

    print "Count is $$count\n";

    $$history = $$gpeptfactory->next_History;

    print "The file for appending the retrieved results is $$file.\n";
}
##################################################################################
##################################################################################
# fetch data from ncbi
##################################################################################
##################################################################################
sub eutilities_fetch {
	my $this = shift;
	my $gpeptfactory = shift;
	my $history = shift;
	my $retmax = shift;
	my $retstart = shift;
	my $retry = shift;
	my $count = shift;
#	my $count = 500;
	my $temp_file = shift;
	my $file = shift;
	my $adminEmail = shift;
	
	
#Now set up to retrieve the records - this code was updated 090507 to work with
#the eUtilities that was updated some time in the past few months

    $$gpeptfactory = Bio::DB::EUtilities->new(
        -eutil   => 'efetch',
        -db      => 'protein',
        -rettype => 'gp',
        -retmode => 'text',
        -tool    => 'VKCDB_Update',
        -email   => $$adminEmail,
        -history => $$history,

    );
    ( $$retmax, $$retstart ) = ( 100, 0 );
    $$retry = 0;

 #retrieve the results in blocks of 100 records to avoid overloading the servers

  DOWNLOAD: while ( $$retstart < $$count ) {

        $$gpeptfactory->set_parameters(
            -retmax   => $$retmax,
            -retstart => $$retstart,
            -rettype => 'gp',
        	-retmode => 'text',

        );

#This evaluation is set up to allow the script to try to get the results 5X since the
#server or the transfer can fail easily and otherwise the script stops with an error

      RETRIEVE_SEQS:
        eval { $$gpeptfactory->get_Response( -file => $$temp_file ) };

        #        print "Output file is $$file.\n";
        if ($@) {
            die "Server error on download: $@.  Try again later, retrieval not working" if $$retry == 5;
            print STDERR "$@\n";
            print STDERR "Server error on download, redo #$$retry\n";
            $$retry++;
            system ("rm $$temp_file\n");
            sleep(5);
            goto DOWNLOAD;
        }

#Once the results are successfully retrieved the contents of the temporary holding file
#$temp_file are read in and written to the file that holds the complete results of
#the search, $$file.

        else {
            open( TEMP, "$$temp_file" ) || die "Can not open $$temp_file.\n";
            while (<TEMP>) {
                if (/temporarily unavailable/si) {
                    close TEMP;
                    redo DOWNLOAD;
                }
            }
            my $retend = $$retstart + $$retmax - 1;
            open( TEMP, "$$temp_file" ) || die "Can not open $$temp_file\n";
            open( DATA, ">>$$file" )    || die "Can not open $$file.\n";
            while (<TEMP>) {
                print DATA $_;
            }
            close TEMP;
            close DATA;
            print
              "Loaded entries $$retstart through $retend on retry #$$retry.\n";

            $$retstart += $$retmax;
            $$retry = 0;
        }
    }
}
##################################################################################
##################################################################################
#this subroutine takes an array of gi's and downloads all there summary info and tests each of those summary info's if there status is suspended if it is don't add it to the new clean array. After clean array is populated download all the genbank records for these files.
##################################################################################
##################################################################################
sub eutilities_getData{

	my $this = shift;

	my $file = shift;

	my $ginumber_search = shift;
	my $vkcid = shift;
	my $lastmod = shift;
	my $vkcnt_id = shift;
	my $adminEmail = shift;
	my $retry;	



#array to hold all records that are not suppressed
	my @clean_array = ();
	#array to hold all gi's from the file to clean up
	my @original_array = @$ginumber_search;
	#file to print all suppressed gi's to
	open (my $myfile, '>', 'removed_gis.txt') || die "Can't open file:$!";
	#file to print all not repressed records to
	open (my $out, '>', $$file) || die "Can't open file:$!";
	print "Finished first part, ready to start retrieving data from GB.\n";
	#initial number of records to post
	my $number = scalar(@original_array);
#	my $number = 600;
	print "Initial number of records of array before sub dividing ".$number."\n";
	#starting and ending index of which records to post from the original array
	my $min = 0;
my $max = 100;
	
my $step = $max - $min;
	#subarray to hold the segment of records specified by min and max
	my @subarray;
	#while loop to post and retreive chunks of the original array
	while ($min < $number){
		if ($max < $number){
			print "start and end subarray index positions: $min     $max\n";
			@subarray = @original_array[$min...$max];
		}
		else{
			@subarray = @original_array[$min...($number-1)];
			print "end of array last post, index positions: $min    ".($number-1)."\n";
		}
		print "subarray had this many elements ".scalar(@subarray)."\n";

	#retry counter to handle the eval block of post and fetch
	my $retry = 0;
	REDO_POST_FETCH:
	eval{
	my $factory = Bio::DB::EUtilities->new(-eutil => 'epost',
													-email => 'wgallin@ualberta.ca',
			
                                       -db => 'protein',
			
                                       -id => \@subarray,
                                       				-verbose => 2,
                                       				
);
                print "finished second part, post of all gis in subarray of original array\n";
                


###check the status of the record if suppressed dont add it to the new gi array###

	
                # get history from queue
                
	my $hist = $factory->next_History || die 'No history data returned';
                
	print "History returned from subarray post\n";
                

	print "beginning of third part, summary retrevial\n";
                

	$factory->set_parameters(-eutil => 'esummary',
                
        	                 -history => $hist,
        	                				 -verbose => 2,);
        	        print "beginning of fourth part, going through each summary record and testing\n";
        	        while (my $ds = $factory->next_DocSum) {
    
		
        	        	my $gi_id = $ds->get_id;
		
						if($ds->get_contents_by_name('Status')){
							my $contents = $ds->to_string;
							my ($status) = $ds->get_contents_by_name('Status');
							#print "Status is $status\n";

		if($status eq "live"){
								push (@clean_array, $gi_id);
						
						
	}
							elsif(($status eq "suppressed") or ($status eq "dead")){

								print $myfile "$gi_id, Status, $status\n";
								#print  "$gi_id, Status, $status\n";
							}
							elsif($status eq "replaced"){
								print  "$gi_id, Status, $status\n";
								my ($replaced) = $ds->get_contents_by_name('ReplacedBy');
								print "ReplacedBy returns $replaced\n";
								my $temp_factory = Bio::DB::EUtilities->new(-eutil	=>	'efetch',
																			-db		=>	'protein',
																			-id		=>	$replaced,
																			-email	=>	'wgallin@ualberta.ca',
																			-rettype	=>	'gi');
								my $latest_gi = $temp_factory -> get_Response -> content;
								chomp $latest_gi;
								#print "The New gi is returned as $latest_gi\n";											
					print $myfile "$gi_id, Status, $status with $latest_gi\n";
						push (@clean_array, $latest_gi);
					
						$$vkcid{ $latest_gi }   = $$vkcid{$gi_id};
                		$$lastmod{ $latest_gi } = $$lastmod{$gi_id};
                		$$vkcnt_id{ $latest_gi} = $$vkcnt_id{$gi_id};
					}
				}	

		}
		#end of docsum while loop
print "finished fourth part\n";
	};
	if ($@) {
	die "Server error on post and fetch: $@.  Try again later, trouble with testing of records\n" if $retry == 5;
    print STDERR "$@\n";
    print STDERR "Server error on post and fetch, redo #$retry\n";
    $retry++;
    sleep(5);
    goto REDO_POST_FETCH;
	}
	
#increment the sub index pointers
$min = ($max +1);
$max += $step;
}
#end of while loop of $min < $number
close $myfile;
@$ginumber_search = @clean_array;
print "Beginning of fifth part, epost of cleaned array\n";
my $number_clean = scalar(@clean_array);
#my $number_clean = 600;
print "Initial number of records of array before sub dividing ".$number_clean."\n";
#starting and ending index of which records to post from the original array
my $min_clean = 0;
my $max_clean = 100;
my $step_clean = $max_clean - $min_clean;
#subarray to hold the segment of records specified by min and max
my @subarray_clean;
#cumulaitve count to keep track of if you are at the end of the original amount of records.
#while loop to post and retreive chunks of the original array
while ($min_clean < $number_clean){
	if ($max_clean < $number_clean){
		print "start and end subarray index positions: $min_clean     $max_clean\n";
		@subarray_clean = @clean_array[$min_clean...$max_clean];
	}
	else{
		@subarray_clean = @clean_array[$min_clean...($number_clean-1)];
		print "end of array last post, index positions: $min_clean    ".($number_clean-1)."\n";
	}
	print "subarray had this many elements ".scalar(@subarray_clean)."\n";


$retry = 0;
my $count = scalar(@subarray_clean);
print "number of records in clean sub array: $count\n";
my $factory;
	
	REDO_POST_FETCH_OF_CLEAN_RECORDS:
	eval{
	$factory = Bio::DB::EUtilities->new(-eutil => 'epost',
		
                                       -email => $$adminEmail,
		
                                       -db => 'protein',
		
                                       -id => \@subarray_clean,
                                       			-verbose => 2,
        
                                       );
    	# get history from queue
    	my $hist = $factory->next_History || die 'No history data returned';
    	
print "History returned from clean array post\n";
    	

	print "beginning sixth part, efetch of clean subarray post\n";
    	
	# note db carries over from above
    	
	$factory->set_parameters(-eutil => 'efetch',
    	
                         -rettype => 'gb',
    	
                         -history => $hist,
                         -verbose => 2,);
		};
	if($@){
	die "Server error on post and fetch of clean records: $@.  Try again later" if $retry == 5;
    print STDERR "$@\n";
    print STDERR "Server error on post and fetch of clean records, redo #$retry\n";
    $retry++;
    sleep(5);
    goto REDO_POST_FETCH_OF_CLEAN_RECORDS;
	}
 
#increment the sub index pointers
$min_clean = ($max_clean +1);
$max_clean += $step_clean;

 
$retry = 0;

my ($retmax, $retstart) = (500,0);

 
	RETRIEVE_SEQS:
	while ($retstart < $count) {
		
    $factory->set_parameters(-retmax => $retmax,
		
                            -retstart => $retstart,
                            		-verbose => 2,
                            		);
                eval{
                
   $factory->get_Response(-cb => sub {my ($data) = @_; print $out $data});
                 };
                if ($@) {
                
        die "Server error at RETRIVE SEQS: $@.  Try again later" if $retry == 5;
        print STDERR "$@\n";
		print STDERR "Server error on retrive seqs of clean records, redo #$retry\n";
        $retry++;
		sleep(5);	
		redo RETRIEVE_SEQS;
    }
    print "Retrieved on attempt $retstart\n";
    $retstart += $retmax;
	}

		
		}
		close $out;}


##################################################################################
##################################################################################

##################################################################################
##################################################################################
sub  bio_seqIO {
    my $this = shift;
    my $retr_seq = shift;
    my $file = shift;
    my $change_flag = shift;
    my $vkcid = shift;
    my $lastmod = shift;
    my $rec_rem = shift;
	my $count = shift;
#	my $count = 500;
	my $ginumber_search = shift;
	my $retry = shift;
	
	
#A Bio::SeqIO object is created to take the entries form the holding file one entry at a
#time and parse the results to check for records that have been replaced by new records.
#$change_flag is set to 0, and will only be reset to 1 if ther is a need to download
#a replacement record.

    $$retr_seq = Bio::SeqIO->new(
        -file   => $$file,
        -format => 'genbank',
    );

    my $incr = 0;
    $$change_flag = 0;
    my @newgi = ();
  RECORD_TEST: while ( my $record = $$retr_seq->next_seq ) {
        my $current_gi      = $record->primary_id;
        my $anno_collection = $record->annotation;
        my @comments        = $anno_collection->get_Annotations('comment');
        foreach my $comment (@comments) {
            my $contents = $comment->as_text;

            #print "Parsing\n$contents\n";

#This is statement identifies any records that have been replaced by a more current
#record attaches the new gi number and date of last modification to the correct
#vkc_id number.  The $change_flag value is reset to 1 to force a new round of record
#retrieval from NCBI

            if ( $contents =~
                /\[WARNING\].*?this sequence was replaced by gi\:(\d+)/s )
            {
                $newgi[$incr] = $1;

                print "$current_gi replaced by $newgi[$incr]\n";

                $$vkcid{ $newgi[$incr] }   = $$vkcid{$current_gi};
                $$lastmod{ $newgi[$incr] } = $$lastmod{$current_gi};
                goto ENDING;
            }

#If the record was removed, then it needs to be manually annotated in VKCDB, but there
#is no longer a valid gi number for that vkc_id number, so no gi number is stored in the
#updated array of gi numbers

            if ( $contents =~ /this sequence was removed/ ) {
                print "Record gi$current_gi was removed, check manually.\nNotprocessing this entry any further.\n";
                $$rec_rem++;
                goto RECORD_TEST;

            }

        }

        #If none of the comments triggers the two conditions above
        #then the gi number is not changed

        $newgi[$incr] = $current_gi;
ENDING:        $incr++;
    }    
    my $new_count = scalar (@newgi);
    if ( ($$count - $$rec_rem) != $new_count ) {
        print
"Mismatch between initial number of records, $$count,\nand number of records finally retrieved, $new_count, so repeating.\n";
	for (my $n = 0; $n<$$count; $n++){
		print "$n     $$ginumber_search[$n]     $newgi[$n]\n";
	}
		@newgi = ();
        $$change_flag = 1;
        $incr = 0;
        $$retry = 0;
        
#TEST LINE
print "this should not happen \n";
        exit;
        
########
        goto RETRIEVE_HIST;
    }

    @$ginumber_search = @newgi;
   #test line
    print
      "Found a total of $incr records when searching for replaced entries.\n";
      $incr = 0;
    $$retr_seq->close;
}
    

1;