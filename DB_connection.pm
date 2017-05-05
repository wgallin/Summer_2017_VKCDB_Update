#!/usr/bin/perl 
package DB_connection;
use lib "/opt/local/lib/perl5/site_perl/5.12.3/";

use strict;
use warnings;
use DBI;

#Create database handle to VKCDB
my $dbh = DBI->connect( "dbi:mysql:database=vkcdb_feb12_2013;host=iolaus.biology.ualberta.ca",
    "vkcdb", "blivkcdb" );

sub new {
    my $class = shift;
    my $self = {};
    
    bless $self, $class;
    return $self;
}

#return connection
sub open_db_connection {

    return $dbh;
}

#close connection
sub close_db_connectin{
	
	$dbh->disconnect();
}

1;
