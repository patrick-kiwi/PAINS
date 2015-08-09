#!/usr/bin/perl
use strict;
use warnings;
use Chemistry::File::SMARTS;
use Chemistry::File::SDF;
use Chemistry::Ring 'aromatize_mol';

my $sdfFile = $ARGV[0];
my @molecules = Chemistry::Mol->read("$sdfFile");


my @painsFiles = qw/PAINS.sieve
		   /;


my $pains;  	#Hash ref 'readable warning' => 'pattern object'
my $painHits; 	#HoA_ref 'molecule_ID' => [ warnings ] 



loadPAINS(\@painsFiles); #load the PAINS to memory 
searchMols(\@molecules); #run the PAINS search
writeOutput(\@molecules); #write a new sdf file with the pains hits




sub writeOutput {
my $outfile = $sdfFile;
$outfile =~ s/(.*)\.sdf/$1_PAINS\.sdf/;
	foreach my $sdfMol ( @{$_[0]} ) {
	my $snid = $sdfMol->attr("sdf/data")->{SNID};

if (ref($sdfMol->attr("sdf/data")->{'comments'}) eq 'ARRAY' && $painHits->{$snid} ) {
push (@{$sdfMol->attr("sdf/data")->{'comments'}}, @{$painHits->{$snid}})
} else {
$sdfMol->attr("sdf/data")->{'comments'} .= "\n$_" for (@{$painHits->{$snid}});
}	
					
Chemistry::Mol->write( $outfile, mols => \@{$_[0]} );
}
}

	 


sub searchMols {
	foreach my $sdfMol ( @{$_[0]} )  {
	my $snid = $sdfMol->attr("sdf/data")->{SNID}; 	#Some unique SDF identifier.  In my case this is a field called SNID
	#next if ( $painHits->{$snid} );  #Skip if it's already been red-flagged
		foreach my $smartsString ( keys %$pains ) {
			if ( $smartsString =~ /\p{L}/ ) {
			aromatize_mol($sdfMol);
				if ( $pains->{$smartsString}->[1]->match($sdfMol) ) {
				push ( @{$painHits->{$snid}}, $pains->{$smartsString}->[0] ); #red flat it!
				print "Hit: ", $snid, " Warning: ", $pains->{$smartsString}->[0], "\n";
				}
			} else {
			if ( $pains->{$smartsString}->[1]->match($sdfMol) ) {
				push ( @{$painHits->{$snid}}, $pains->{$smartsString}->[0] ); #red flat it!
				print "Hit: ", $snid, " Warning: ", $pains->{$smartsString}->[0], "\n";
			}
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
my $smartsString = $2;
my $patternObject = Chemistry::Pattern->parse("$smartsString", format => 'smarts' );
$pains->{$smartsString}->[0] = $readableWarning;
$pains->{$smartsString}->[1] = $patternObject;
}
close $inFH;
}
print "*** Done loading PAIN filters ***\n\n";
}
