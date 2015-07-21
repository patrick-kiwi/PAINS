#!/usr/bin/perl
use strict;
use warnings;
use Chemistry::File::SMARTS;
use Chemistry::File::SDF;

my $sdfFile = $ARGV[0];
my @molecules = Chemistry::Mol->read("$sdfFile");

my @painsFiles = qw/PAINS_FilterFamilyA.sieve 
		    PAINS_FilterFamilyB.sieve
		    PAINS_FilterFamilyC.sieve
		   /;

my $pains;  	#Hash ref 'readable warning' => 'pattern object'
my $painHits; 	#HoA_ref 'molecule_ID' => [ warnings ] 



loadPAINS(\@painsFiles); #load the PAINS to memory 
searchMols(\@molecules); #run the PAINS search


sub searchMols {
	foreach my $sdfMol ( @{$_[0]} )  {
	my $snid = $sdfMol->attr("sdf/data")->{SNID}; 	#Some unique SDF identifier
	next if ( $painHits->{$snid} );  #Skip if it's already been red-flagged
		foreach my $warning ( keys %$pains ) {
		if ( $pains->{$warning}->match($sdfMol) ) {
		push ( @{$painHits->{$snid}}, $warning ); #red flat it!
		print "Hit: ", $snid, " Warning: ", $warning, "\n";
	}
}
}
}


sub loadPAINS {
foreach my $infile ( @{$_[0]} ) {
print "*** Loading filter from file: $infile\n";
open (my $inFH, "<", $infile) 
	|| die "Cannot open PAINS filter file: $!\n";
while (<$inFH>) {
next unless /^FRAGMENT\s+regId=(\S+)\s+(\S+)\s+/;
my $readableWarning = $1;
my $patternObject = Chemistry::Pattern->parse("$2", format => 'smarts' );
$pains->{$readableWarning} = $patternObject;
}
close $inFH;
}
print "*** Done loading PAIN filters ***\n";
}
