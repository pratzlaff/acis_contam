#! /usr/bin/perl

use strict;
use warnings;

use PDL;
use PDL::Graphics::PGPLOT;
use Chandra::Tools::Common 'read_bintbl_cols';

my $xrange = [0, 70];

my $dev = 'trans_s3_2015.ps/vcps';
#$dev = '/xs';

dev $dev, 1, 3, { hardlw => 2, hardch => 1.25 };
autolog(1);

my ($wav_lo, $wav_hi, $specresp) = read_bintbl_cols('/data/legs/rpete/flight/acis_letg_cedge/data/mkn421/17682/tg_reprocess_osort1/acisf17682N002LEG_1_garf.fits', qw/ bin_lo bin_hi specresp /, { extname => 'specresp' });

line
  0.5*($wav_lo + $wav_hi), $specresp,
  {
#   axis => 'LOGY',
   border => 1,
   xrange => $xrange,
   ytitle => 'EA (cm\u2\d)',
   title => 'Mkn 421 SPECRESP, Obsid 17682, TG_M=+1',
   };


my %data;

my @chipy = (180, 300, 400, 512);

my %ci = ( 180 => 1, 300 => 2, 400 => 3, 512 => 5 );

for my $chipy (@chipy) {
  $data{$chipy} = { };

  my ($e, $t) = rcols("trans_s3_chipy_$chipy.txt");

  $data{$chipy}{energy} = $e;
  $data{$chipy}{trans} = $t;
  $data{$chipy}{wavelength} = 12.39854 / $data{$chipy}{energy};

  line $data{$chipy}{wavelength}, $data{$chipy}{trans},
      {
       border => 1,
       xrange => $xrange,
       yrange => [0.01, 1],
       axis =>'LOGY',
       color => $ci{$chipy},
       ytitle => "Contamination Transmittance",
       title => "ACIS-S3, 2015.0, CHIPX=512",
      };
  hold;
}

legend
  [map { "CHIPY=$_" } @chipy],
  30, log10(.9),
  {
   color=>[@ci{@chipy}],
  };

release;


my $black = 0;
for my $chipy (300, 400, 512) {

  my $ci;

  if ($chipy == 300 and not $black) {
    $ci = 1;
  }

  else {
    $ci = $ci{$chipy}
  }

  line $data{$chipy}{wavelength}, $data{$chipy}{trans} / $data{180}{trans},
      {
       border => 1,
       xrange => $xrange,
       yrange => [1, 2],
       color => $ci,
       xtitle => "\\gl (\\A))",
       ytitle => "Transmittance Ratio vs CHIPY=180",
       title => "ACIS-S3, 2015.0, CHIPX=512",
      };
  hold;

  if ($chipy == 300 and not $black) {
    $black = 1;
    redo;
  }

}

legend
  [map { "CHIPY=$_" } (300, 400, 512)],
  50, 1.7,
  {
   color=>[@ci{300,400,512}],
  };

release;
