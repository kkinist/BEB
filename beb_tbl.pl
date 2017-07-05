#!/usr/bin/perl
#
# read a tab-delimited table of molecular data, compute BEB cross sections
# each row: MO label, B (eV), U (eV), N, Q (default = 1),
#   'dbl' ("yes" for dbl ionization, default is no), 'special' (for special treatment),
#   'remarks' (for comments)
# comment lines in data file begin with '#'
# Karl Irikura, NIST
#
# First version:  6/9/2000 KKI
# Add optional 2nd cmd-line arg for single incident KE (6/15/00)
# add numerical derivative, turned on by the flag $deriv inside the script (7/15/03)
# add peak value (10/15/04)
# send (energy, xsec) data to a file; add file suffix when missing (7/5/17)
#
$i = 0;
$tmin = 9999.0;
$comment = "";
$deriv = 0;	# flag to compute and print derivative
$dt = 0.01;	# for numerical derivative
( $fname, $just_one, $details ) = @ARGV;
if (! -f $fname) {
    # try adding suffix ".bun"
    print "*** File \"$fname\" not found: adding suffix \".bun\" ***\n";
    $fname .= ".bun";
}
open( DATA_FILE, "<$fname" ) or die "Can't open data file $fname!\n";
# create name for output file (energy, xsec), CSV
$fcsv = $fname;
$fcsv =~ s/\..*$//;
$fcsv .= "_BEBxsec.csv";
open( OUTPUT_FILE, ">$fcsv" ) or die "Failure creating BEB output file $fcsv\n";
while ( <DATA_FILE> ) {
	if ( /^#/ ) {
		# comment line; append to $comment
		$comment .= $_;
		next;
	}
	chomp;
	@word = split( /\t/, $_ );
	next if ( @word < 4 );		# skip bad or blank lines (require at least lbl, b, u, n)
	( $lbl[$i], $b[$i], $u[$i], $n[$i], $q[$i], $dbl[$i], $special[$i], $remarks[$i] ) = @word;
	$q[$i] = 1 unless( $q[$i] );	# default q=1
	if ( $dbl[$i] =~ /^[yY]/ ) {
		$dbl[$i] = 'Y';	# upper-case 'Y'
	} else {
		$dbl[$i] = 'n';	# lower-case 'n'
	}
	$tmin = $b[$i] if ( $tmin > $b[$i] );
	$i++;
}
close( DATA_FILE );

# echo comment and orbital info
print "$comment\n";
print "MO\tB\tU\tN\tQ\tdbl\tspecial\tremarks\n";
for ( $i = 0; $i < @b; $i++ ) {
	print "$lbl[$i]\t$b[$i]\t$u[$i]\t$n[$i]\t$q[$i]\t$dbl[$i]\t$special[$i]\t$remarks[$i]\n";
}
print "\n";

$ntot = 0;
foreach ( @n ) { $ntot += $_; }
print "Total number of electrons = $ntot\n\n";

# compute BEB/BEQ cross section
$tpeak = $xpeak = 0;	# energy and cross-section at peak
$tmax = 5000;
if ( $just_one ) {
	@t = ( $just_one );
} else {
	@t = t_grid( $tmin, $tmax );	# generate list of kinetic energies to compute
}
$i = 0;		# flag for choosing to print: "BEB" or "BEQ"
foreach ( @q ) { $i = 1 if ( $_ != 1 ); }
unless ($details) {
    print " Energy (eV)\t ", ( $i ? "BEQ" : "BEB" ), " (A**2)";
}
print ( $deriv ? "\t gradient (A**2/eV)\n" : "\n" );
for ( $i = 0; $i < @t; $i++ ) {			# loop over incident energies
	$sigma[$i] = fsigma( $t[$i] );
	if ( $deriv ) {
		# numerical computation of derivative wrt incident energy
		$dsigma[$i] = fsigma( $t[$i] + $dt ) - fsigma( $t[$i] - $dt );
		$dsigma[$i] /= 2 * $dt;
	}
    if ($details) {
        print "\n Energy (eV)\t ", ( $i ? "BEQ" : "BEB" ), " (A**2)\n";
    }
	printf " %9.2f \t %7.3f", $t[$i], $sigma[$i];
	printf OUTPUT_FILE "%.4f,%.4f\n", $t[$i], $sigma[$i];
	if ( $deriv ) {
		printf "\t%10.6f\n", $dsigma[$i];
	} else {
		print "\n";
	}
	if ( $sigma[$i] > $xpeak ) {
		$xpeak = $sigma[$i];
		$tpeak = $t[$i];
	}
}
close OUTPUT_FILE;
print "Numerical output written to file $fcsv\n";
unless ( $just_one ) {
	printf "\nPeak cross-section is %.3f at %.0f eV.\n", $xpeak, $tpeak;
}

exit(0);
# print in multiple columns 24 lines high
print "\n";
for ( $i = 0; $i < 24; $i++ ) {
	for ( $j = $i; $j < @t; $j += 24 ) {
		printf "%9.2f %7.3f ", $t[$j], $sigma[$j];
	}
	print "\n";
}

sub fsigma {
	# arg is incident electron energy (in eV)
	# return value is total cross section
	my ( $T ) = @_;
	my( $j, $sigma, $tsigma );
    my $oflag = 1;
	$tsigma = 0;
	for ( $j = 0; $j < @b; $j++ ) {		# loop over orbitals
		next if ( $T < $b[$j] );
		$sigma = beq( $T, $b[$j], $u[$j], $n[$j], $q[$j], $special[$j] );
		$sigma *= 2 if ( $dbl[$j] =~ /^[yY]/ );
		if ( $details ) {
			# print contribution from each orbital
            if ($oflag) {
                print "Orbital contributions:\n";
                print "MO\tB\txsec\n";
                $oflag = 0;
            }
			printf "%d\t%.2f\t%f\n", $lbl[$j], $b[$j], $sigma;
		}
		$tsigma += $sigma;
	}
	return $tsigma;
}

sub t_grid {
	# args are lower and upper limits
	# generate list of incident energies
	# interval is in %interv (key = incident eV lower end, value = interval)
	my( $e, $emax ) = @_;
	my( $de, @list, @sorted );
	%interv = (
		24	=>	2,
		40	=>	5,
		150	=>	10,
		250	=>	50,
		1000 =>	100,
		3000 =>	200,
		5000 =>	500,
	);
	@sorted = sort { $a <=> $b } keys %interv;
	while ( $e <= $emax ) {
		push( @list, $e );
		$de = 0.5;	# interval for very low energies
		foreach ( @sorted ) {
			$de = $interv{ $_ } if ( $e >= $_ );
		}
		$e += $de;
		unless ( $#list ) {
			# this is the first point above the lower limit
			# round $e upward from initial value
			$e -= $de;
			$e = $de * int( 1.1 + $e / $de );
		}
	}
	return @list;
}

sub beq {
	# args are: T (eV), B (eV), U (eV), N, Q, 'special'
	# valid 'special' values:
	#	'ion'						treat as +1 ion: (u+1)/2
	#	nl (e.g., '3s' or '3p')		treat as n>2 orbital: use (u+1)/n
	#	'heavy_core'				treat as n<3 heavy core orbital: multiply by [1+(t+u+1)/(t)]/2
	#
	# compute BEQ cross section
	my( $pi, $a0, $r ) = ( 3.141592654, 0.5291772, 13.6057 );
	my( $T, $B, $U, $n, $q, $special ) = @_;
	my( $t, $u, $s, $x );
	$t = $T / $B;
	$u = $U / $B;
	$s = $a0 * $r / $B;
	$s = 4 * $pi * $n * $s * $s;
	$x = 1 - 1 / $t - log( $t ) / ( $t + 1 );
	$x *= 2 - $q;
	$x += $q * log( $t ) * ( 1 - 1 / ( $t * $t ) ) / 2;
	$x *= $s;
	if ( $special eq "ion" ) {
		# treat as +1 ion
		$x /= $t + ( $u + 1 ) / 2;
	} elsif ( $special =~ /^(\d)[spdf]{1,4}\b/i ) {
		# treat as atom-localized valence orbital with n = $special
		$x /= $t + ( $u + 1 ) / $1 if ( $1 > 2 );
	} else {
		# ordinary and 'heavy_core' cases
		$x /= $t + $u + 1;
	}
	if ( $special eq "heavy_core" ) {
		# core orbital (K or L shell) for heavy atom 
		$x *= 1 + ( $u + 1 ) / ( 2 * $t );
	}
	return $x;
}
