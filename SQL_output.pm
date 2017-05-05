#!/usr/bin/perl 
package SQL_output;

use strict;
use warnings;
use lib "/opt/local/lib/perl5/site_perl/5.12.3/";

use LWP::Simple;
use DBI;
use Bio::DB::GenPept;
use Bio::DB::Taxonomy;
use Bio::Seq;
use Bio::DB::EUtilities;
use Bio::DB::Taxonomy::entrez;



sub new {
    my $class = shift;
    my $self = {};
    
    bless $self, $class;
    return $self;
}

sub P_or_E {
	my $this = shift;
	my $dbh = shift;
    
    my $sth = $$dbh->prepare("SELECT  org_id, taxon, P_or_E FROM organism")
      || die "Can't prepare. $$dbh->errstr\n";
    my $rv = $sth->execute || die "Can't execute the query: $sth->errstr\n";

    #Open connection to Entrez Taxonomy DB

    my $tdbh = new Bio::DB::Taxonomy( -source => 'entrez' );
    my $ancestor;

    #Open output file for the SQL update commands

    open( OUT, ">>division_update" );

    #Retrieve all the taxon IDs from VKCDB
    while ( my @row = $sth->fetchrow_array ) {
        if ( defined( $row[2] ) ) {
            next;
        }
        my $org_id         = $row[0];
        my $taxon_id       = $row[1];
        my $lineage_object = $tdbh->get_taxon( -taxonid => $taxon_id );

        while ( defined( $ancestor = $tdbh->ancestor($lineage_object) ) ) {

            #Retrieve the full division from the taxon object
            my $name = $ancestor->node_name;
            print "$org_id is $name\n";

            if ( $name eq "Bacteria" ) {

                print OUT
"UPDATE organism SET P_or_E = \"B\" WHERE org_id = \"$org_id\";\n";
                print "$org_id is a Prokaryote\n";

            }
            elsif ( $name eq "Eukaryota" ) {

                print OUT
"UPDATE organism SET P_or_E = \"E\" WHERE org_id = \"$org_id\";\n";
                print "$org_id is a Eukaryote\n";

            }
            elsif ( $name eq "Archaea" ) {

                print OUT
"UPDATE organism SET P_or_E = \"A\" WHERE org_id = \"$org_id\";\n";
                print "$org_id is an Archaea\n";

            }
            $lineage_object = $ancestor;

        }
    }

    close OUT;
   
}
1;