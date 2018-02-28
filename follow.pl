#!/usr/bin/perl -w
#  Follow a G09 normal mode (probably imaginary freq.)
#  KKI 2/16/2012
#
# element symbols 
@element = qw(
	X  H  He
	Li Be B  C  N  O  F  Ne
	Na Mg Al Si P  S  Cl Ar
	K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
	Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe
	Cs Ba La                                              Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu
			 Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
	Fr Ra Ac                                              Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No Lr
);
#
if ( $#ARGV < 0 ) {
	die "Usage:  perl follow.pl <g09.out> [mode] [displacement] [MOL] [-v]\n";
}
open( G09, "<$ARGV[0]" ) or die "Failure reading file $ARGV[0]\n";
$target_mode = 0;	# default is to follow the first (lowest-frequency) mode
$displ = 0.5;	# default displacement along the target mode
$mol = 0;	# flag to print differently
$verbose = 0;	# flag for more output
shift( @ARGV );
foreach ( @ARGV ) {
	if ( /^\d+$/ ) {
		# an integer; set as target mode
		$target_mode = $_ - 1;
	}
	if ( /^[-]?(\d+)?\.\d+$/ ) {
		# decimal point; set as displacement
		$displ = $_;
	}
	if ( /MOL/i ) {
		$mol = 1;
	}
	if ( /-v/i ) {
		$verbose = 1;
	}
}
printf "target mode = %d, displacement = %f\n", ( $target_mode + 1, $displ ) if ( $verbose );
$instd = $inmode = 0;	# flags
while ( <G09> ) {
	$instd = 0 if ( /Rotational constants|Distance matrix / );
	$inmode = 0 if ( /Thermochemistry/ );
	if ( /Harmonic frequencies / ) {
		@freq = @nmx = @nmy = @nmz = ();
	}
	if ( $instd && /^\s+\d+\s+/ ) {
		@line = split();
		$i = $line[0] - 1;
		$nuc[$i] = $line[1];
		$x[$i] = $line[3];
		$y[$i] = $line[4];
		$z[$i] = $line[5];
		$natom = $line[0];
	}
	if ( $inmode && /\d+\.\d+/ && not /[a-zA-Z]/ ) {
		@line = split();
		$iat = shift( @line ) - 1;
		shift( @line );	# discard atomic number
		$imode = $#freq - int($#line / 3);
		for ( $j = $imode; $j <= $#freq; $j++ ) {
			$nmx[$j]->[$iat] = shift( @line );
			$nmy[$j]->[$iat] = shift( @line );
			$nmz[$j]->[$iat] = shift( @line );
		}
	}
	if ( /(Standard|Input|Z-Matrix) orientation:/ ) {
		$instd = 1;
		@x = @y = @z = @nuc = ();
	}
	if ( /Frequencies --/ ) {
		$inmode = 1;
		@line = split();
		shift( @line );
		shift( @line );
		push( @freq, @line );
	}
}
close( G09 );
if ( $verbose ) {
	print "Original coordinates:\n";
	for ( $i = 0; $i <= $#nuc; $i++ ) {
		printf "%-2s%12.6f%12.6f%12.6f\n", ( $element[$nuc[$i]], $x[$i], $y[$i], $z[$i] );
	}
	if ( $mol ) {
		print "Displaced coordinates (MOL format):\n";
	} else {
		print "Displaced coordinates:\n";
	}
}
for ( $i = 0; $i <= $#nuc; $i++ ) {
	$dx = $displ * $nmx[$target_mode]->[$i];
	$dy = $displ * $nmy[$target_mode]->[$i];
	$dz = $displ * $nmz[$target_mode]->[$i];
	$dx += $x[$i];
	$dy += $y[$i];
	$dz += $z[$i];
	if ( $mol ) {
		# MOL-file format
		printf "%10.5f%10.5f%10.5f %-3s", ( $dx, $dy, $dz, $element[$nuc[$i]] );
		print "  0  0  0  0  0  0  0  0  0  0  0  0\n";
	} else {
		# ordinary format
		printf "%-2s%12.6f%12.6f%12.6f\n", ( $element[$nuc[$i]], $dx, $dy, $dz );
	}
}
exit(0);
# print the frequencies and normal mode vectors
for ( $i = 0; $i <= $#freq; $i++ ) {
	printf "%3d%12.4f\n", ( $i+1, $freq[$i] );
	for ( $j = 0; $j < $natom; $j++ ) {
		printf "\t%8.2f%8.2f%8.2f\n", ( $nmx[$i]->[$j], $nmy[$i]->[$j], $nmz[$i]->[$j] );
	}
}
