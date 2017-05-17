#!/usr/bin/perl -w
#  print, in readable form, the archive block from a Gaussian output file
# KKI 3/28/06
#
if ( $#ARGV != 0 ) {
	print "Usage: perl arch.pl <g03.out>\n";
	exit(0);
}
open( G03, $ARGV[0] ) or die "Failure reading file $ARGV[0]\n";
$arch = "";
$inarch = 0;
while ( <G03> ) {
	$inarch = 1 if ( /1\\1\\/ or /1\|1\|/ );
	if ( $inarch ) {
		chomp();
		s/^\s+//;
		s/\s+$//;
		$arch .= $_;
	}
	$inarch = 0 if ( /\\\\\@/ or /\|\|\@/ or /^\s*$/ );
}
# split at the backslashes
if ( $arch =~ /\\\\/ ) {
	@arch = split( /\\+/, $arch );
} else {
	@arch = split( /\|+/, $arch );
}
foreach ( @arch ) {
	# next line improves readability
	s/=/ = /g;
	print "$_\n";
}
