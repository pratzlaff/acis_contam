package Contam;

use 5.006;
use strict;
use warnings;

use Exporter 'import';
our @EXPORT_OK = qw(
		     image
		  );

use Carp;
use Log::Any '$log';
use PDL;

sub year2time {
  my $year = shift;
  return 325.2425*86400*($year-1998);
}

sub image {
  my ($fname, $ccd_id, $energy, $year) = @_;
  my $h = rfits($fname."[@{[$ccd_id+1]}]");
  my $n = +($h->{'component'}->dims)[0];
  my $taus = zeros($h->{'fxy     '});
  my $time = year2time($year);

  for my $i (0..$n-1) {
    my $mu = interpol($energy, $h->{'energy  '}->slice("($i)"), $h->{'mu      '}->slice("($i)"));
    my $tau0 = interpol($time, $h->{'time    '}->slice("($i)"), $h->{'tau0    '}->slice("($i)"));
    my $tau1 = interpol($time, $h->{'time    '}->slice("($i)"), $h->{'tau1    '}->slice("($i)"));
    my $tau = $tau0 + $tau1 * $h->{'fxy     '}->slice("($i)");
    $taus->slice("($i)") .= $mu * $tau;
  }
  return $taus;
}

sub energy {
  my ($fname, $ccd_id, $x, $y, $year) = @_;
  my $energy;
  my $h = rfits($fname."[@{[$ccd_id+1]}]");
  my $n = +($h->{'component'}->dims)[0];
  my $taus = zeros($h->{'energy  '});
  my $time = year2time($year);

  for my $i (0..$n-1) {
    if (not defined $energy) {
      $energy = $h->{'energy  '}->slice("($i)");
    }
    my $tau0 = interpol($time, $h->{'time    '}->slice("($i)"), $h->{'tau0    '}->slice("($i)"));
    my $tau1 = interpol($time, $h->{'time    '}->slice("($i)"), $h->{'tau1    '}->slice("($i)"));
    my $tau = $tau0 + $tau1 * $h->{'fxy     '}->slice("$i,$x,$y");
    $taus->slice("($i)") .= $h->{'mu      '}->slice("($i)") * $tau;
  }
  return $energy, $taus;
}

sub year {
  my ($fname, $ccd_id, $x, $y, $energy) = @_;
  my $times;
  my $ntimes = 101;
  my $h = rfits($fname."[@{[$ccd_id+1]}]");
  my $n = +($h->{'component'}->dims)[0];
  my $taus = zeros($n, $ntimes);

  for my $i (0..$n-1) {
    if (not defined $times) {
      my ($tmin, $tmax) = $h->{'time    '}->slice("($i)")->minmax;
      $times = sequence($ntimes)*($tmax-$tmin)/($ntimes-1)+$tmin;
    }
    my $mu = interpol($energy, $h->{'energy  '}->slice("($i)"), $h->{'mu      '}->slice("($i)"));
    my $tau0 = interpol($times, $h->{'time    '}->slice("($i)"), $h->{'tau0    '}->slice("($i)"));
    my $tau1 = interpol($times, $h->{'time    '}->slice("($i)"), $h->{'tau1    '}->slice("($i)"));
    my $tau = $tau0 + $tau1 * $h->{'fxy     '}->slice("$i,$x,$y");
    $taus->slice("($i)") .= $mu * $tau;
  }
  return 1998+$times/86400/365.2425, $taus;
}

1; # End of Contam.pm
