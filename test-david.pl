use strict;
use warnings FATAL => qw( all );

use SOAP::Lite;
use HTTP::Cookies;

my $soap = SOAP::Lite                             
	-> uri('http://service.session.sample')                
	-> proxy('http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService',
			cookie_jar => HTTP::Cookies->new(ignore_discard=>1));

 #user authentication by email address
 #For new user registration, go to http://david.abcc.ncifcrf.gov/webservice/register.htm
 my $check = $soap->authenticate('christian.frech@ccri.at')->result;
  	print STDERR "User authentication: $check\n";