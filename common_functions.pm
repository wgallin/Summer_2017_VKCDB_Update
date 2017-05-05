#!/usr/bin/perl 
package common_functions;

use strict;
use warnings;
use lib "/opt/local/lib/perl5/site_perl/5.12.3/";


sub new {
    my $class = shift;
    my $self = {};
    
    bless $self, $class;
    return $self;
}

#Capture current date for insertion into protein.lastmod
#then returns that date
  sub get_current_time {
  	my $this = shift;
  	my @timedata = localtime(time);
    my $newdate =
  	( $timedata[5] + 1900 ) . "-" . ( $timedata[4] + 1 ) . "-" . $timedata[3];
  return $newdate; 
}

sub getEmail{
	my $this = shift;
	my $adminEmail = shift;
	print "Please type admin email then hit enter to continue:\n";
	$$adminEmail = <>;
}

sub month_converter {
	my $this = shift;
    my $mon = shift;

    my $digitmon;

    if ( $mon eq "JAN" ) { $digitmon = "01"; }
    if ( $mon eq "FEB" ) { $digitmon = "02"; }
    if ( $mon eq "MAR" ) { $digitmon = "03"; }
    if ( $mon eq "APR" ) { $digitmon = "04"; }
    if ( $mon eq "MAY" ) { $digitmon = "05"; }
    if ( $mon eq "JUN" ) { $digitmon = "06"; }
    if ( $mon eq "JUL" ) { $digitmon = "07"; }
    if ( $mon eq "AUG" ) { $digitmon = "08"; }
    if ( $mon eq "SEP" ) { $digitmon = "09"; }
    if ( $mon eq "OCT" ) { $digitmon = "10"; }
    if ( $mon eq "NOV" ) { $digitmon = "11"; }
    if ( $mon eq "DEC" ) { $digitmon = "12"; }

    return $digitmon;
}



1;