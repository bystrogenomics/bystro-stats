use 5.10.0;
use strict;
use warnings;
use DDP;
use Cpanel::JSON::XS;

open(my $fh, '<', $ARGV[0]);

open(my $writeFh, '|-', "seqant_statistics -referenceColumnIdx 7 -alleleColumnIdx 6 "
  . "-heterozygotesColumnIdx 3 -homozygotesColumnIdx 4 --outJsonPath test.out");

my $lines = '';
my $count = 0;
while(<$fh>) {
  $lines .= $_;

  $count++;

  if($count > 1e4) {
    print $writeFh $lines;

    $lines = '';
  }
}

if($lines) {
  print $writeFh $lines;
}

close $writeFh;

open(my $jsonFh, '<', 'test.out');
my $jsonString = <$jsonFh>;

say "json is";

my $jsonHref = decode_json($jsonString);

p $jsonHref;

close $jsonFh;