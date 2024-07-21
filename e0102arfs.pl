#! /usr/bin/perl -w
use strict;

use PDL;
use PDL::Graphics::PGPLOT;
use Chandra::Tools::Common qw( read_bintbl_cols );
use Astro::FITS::CFITSIO;

my ($elo, $ehi) = (0.31,2);
$_ = log $_ for $elo, $ehi;
my $egrid = exp(sequence(1000) / 1000 * ($ehi-$elo) + $elo);

# I3 arfs

my %specresp;
my @obsids = (440, 3526, 6751, 10649 );

for my $o (@obsids) {

  my %f = (
	   '04' => "I3arfs/obs${o}_ciao412_cir_j2000_specextract.warf",
	   '89' => "I3arfs/obs${o}_ciao4_EAvF_ContamN9989.warf",
	   '31' => "I3arfs/obs${o}_ciao4_EAvF_ContamN9931.warf",
	   );

  for my $k (keys %f) {
    -f $f{$k} or die $f{$k};
    my ($el, $eh, $specresp) = read_bintbl_cols($f{$k},
						qw( energ_lo energ_hi specresp ), { extname => 'specresp' },
					       );
    $specresp{$k}{$o} = interpol($egrid, 0.5 * ($el+$eh), $specresp);
    print $el->min,"\n";
  }

  my $hdr = Astro::FITS::CFITSIO::fits_read_header($f{'04'});
  my $year = $hdr->{TSTART}/86400/365.24 + 1998;
}

my $dev = 'e0102_arf_ratios.ps/ps';
dev $dev, 2,1;
autolog(1);
for my $k ('89', '31') {
  my $ls = 1;
  for my $o (@obsids) {
    line $egrid, $specresp{'04'}{$o} / $specresp{$k}{$o},
      {
       linestyle => $ls,
       axis => 'logx',
       yrange => [0.35,2],
       title => "E0102 I3 arfs : 04 / $k",
       };
    hold;
    ++$ls;
  }
  release;
}
