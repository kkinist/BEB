#!/usr/bin/perl -w
#  return value of Nimag from Gaussian output file
#  KKI 2/16/2012
#
@files = glob "$ARGV[0]";
foreach $file ( @files ) {
	print "$file: ";
	open( NIM, "arch.pl $file |" ) or die "failure opening pipe";
	while( <NIM> ) {
		if ( /nimag/i ) {
			@line = split();
			$nimag = $line[-1];
			print "nimag = $nimag\n";
		}
	}
}
