#!/usr/bin/perl 
package DB_Queries;
use lib "/opt/local/lib/perl5/site_perl/5.12.3/";

use strict;
use warnings;
use DBI;

sub new {
    my $class = shift;
    my $self = {};
    
    bless $self, $class;
    return $self;
}

#Get the last ref_id value from the DB, since new reference records may be created.
#takes in a reference parameter of a scalar variable and a connection object
#sets the new_refindex passed, dereferences it where needed and modifies it
#to be the new value.
sub ref_id {
	my $this = shift;
	my $dbh = shift;#can use $_ as well, makes a local variable equal to the parameter passed to this method
	my $new_refindex = shift;
	
	my $statement;
	my $sth = $$dbh->prepare("SELECT ref_id FROM reference")#gets a statement ready for execution and returns a statement handle OBJECT
  		|| die "Can't prepare $statement:$$dbh->errstr\n";

	my $rv = $sth->execute || die "Can't execute the query: $sth->errstr\n";#executes a given statement or gives an error, $rv is just 
																			#to hold the return true/false/undef value of the execute method.
	my $oldref_id = 0;

	while ( my @row = $sth->fetchrow_array ) {#declares a local variable of @row and makes it equal to the 
    	$row[0] =~ /Ref(\d+)/;					#array that is returned 
    	if ( $1 > $oldref_id ) { #$1 is a special variable that is the first part of the string returned from the =~ evaluation above
        $oldref_id = $1;
    }
}
#calls the finish method, $rc is just to hold the return value of the method. optional? as its 
#an explicit call?
	my $rc           = $sth->finish;
	$$new_refindex = $oldref_id + 1;
}


#Collect the vkc_id, gi number and lastmod date for each reference,
#these will be used for accessing the latest GenPept entries and tracking changes that are made.
sub vkcID_gi_lastMod {
	my $this = shift;
	my $dbh = shift;#can use $_ as well, makes a local variable equal to the parameter passed to this method
	my $ginumber_search = shift;
	my $vkcid = shift;
	my $lastmod = shift;
	my $vkcnt_id= shift;
	my $code_taxon = shift;
	my $count = shift;

	my $sth = $$dbh->prepare("SELECT protein.vkc_id, protein.ginumber, protein.lastmod, organism.taxon, gennt.vkcnt_id  
	FROM gennt 
	INNER JOIN gennt_protein ON gennt.vkcnt_id = gennt_protein.vkcnt_id 
	INNER JOIN protein ON gennt_protein.vkc_id = protein.vkc_id 
	INNER JOIN protein_organism ON protein.vkc_id = protein_organism.vkc_id 
	INNER JOIN organism ON protein_organism.org_id = organism.org_id 
	WHERE flag = \"Y\" AND ginumber!=\"\"")#LIMIT 0 , 2 
	|| die "Can't prepare statement:$$dbh->errstr\n";

	my $rv = $sth->execute
  		|| die "Can't execute the query: $sth->errstr\n";

	my $i = 0;

	while ( my @row = $sth->fetchrow_array ) {
		
	
    $$ginumber_search[$i] = $row[1];
    $row[0] =~ /VKC(\d+)/;
    $$vkcid{ $$ginumber_search[$i] }   = $1;
    $$lastmod{ $$ginumber_search[$i] } = $row[2];
    $$code_taxon{$$ginumber_search[$i]} =$row[3];
	$$vkcnt_id{$$ginumber_search[$i]} = $row[4];
	
	#print "$$ginumber_search[$i]     $$vkcid{ $$ginumber_search[$i] }     $$lastmod{ $$ginumber_search[$i]}     $$code_taxon{$$ginumber_search[$i]}     $$vkcnt_id{$$ginumber_search[$i]} \n";
    #print "$$ginumber_search[$i] \n";
    
    $i++;
}
#$$count = 500;
$$count = scalar (@$ginumber_search) +1;
my $rc = $sth->finish;
}






1;
