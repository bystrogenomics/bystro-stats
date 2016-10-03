use 5.10.0;
use strict;
use warnings;
use DDP;
use Cpanel::JSON::XS;

open(my $fh, '<', $ARGV[0]);

open(my $writeFh, '|-', 'seqant-statistics -outputJSONPath ExAC.hg19.snpOnly-10klines.statistics.json '
  .' -outputTabPath ExAC.hg19.snpOnly-10klines.statistics.tab -outputQcTabPath ExAC.hg19.snpOnly-10klines.statistics.qc.tab '
  .' -referenceColumnName hg38 -alleleColumnName minorAlleles -homozygotesColumnName homozygotes '
  . '-heterozygotesColumnName heterozygotes -siteTypeColumnName refSeq.siteType '
  . ' -dbSNPnameColumnName dbSNP146.name -exonicAlleleFunctionColumnName refSeq.exonicAlleleFunction '
  . " -numberInputHeaderLines 1 -primaryDelimiter \$\";\" -secondaryDelimiter \$\"|\" -fieldSeparator \$\"\t\"");

my $lines = '';
my $count = 0;
while(<$fh>) {
  chomp;
  $lines .= $_;

  
  
  $count++;

  if($count > 1e4) {
    say $writeFh $lines;

    $lines = '';
    sleep(1);
  }
}

if($lines) {
  say $writeFh $lines;
}

close $writeFh;

open(my $jsonFh, '<', 'test.out');
my $jsonString = <$jsonFh>;

say "json is";

my $jsonHref = decode_json($jsonString);

p $jsonHref;

close $jsonFh;