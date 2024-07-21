#! /usr/bin/perl -w
use strict;

use PDL;
use PDL::Graphics::PGPLOT;
use PGPLOT;

dev 'foo.ps/cps', 2, 4, { aspectratio => 4/2 };

for (1..8) {
  pgpage();

  # viewport size of each chip
  my $d = 0.14;
  my $ics = (1-(6*$d))/7; # inter-chip spacing
  my $vm = 0.5 * (1- 3*$d + 2*$ics); # vertical margins

  my ($x1, $x2, $y1, $y2);

  my $box = sub {
    pgsvp($x1, $x2, $y1, $y2);
    pgswin(0.5, 5.5, 0.5, 5.5);
#    pgbox('BC', 0, 0, 'BC', 0, 0);
    my $img = sequence(5,5);
    pggray($img->float->get_dataref,
	   $img->getdim(0),
	   $img->getdim(1),
	   1, $img->getdim(0),
	   1, $img->getdim(1),
	   0, 25,
	   [
	    0.5, 1, 0,
	    0.5, 0, 1,
	   ],
	   );
  };

  # I0
  $y2 = 1-$vm;
  $y1 = $y2-$d;
  $x1 = 3*$ics + 2*$d;
  $x2 = $x1+$d;
  $box->();

  # I1
  $_ += $ics + $d for $x1, $x2;
  $box->();

  # I2
  $_ -= $ics + $d for $y1, $y2, $x1, $x2;
  $box->();

  # I3
  $_ += $ics + $d for $x1, $x2;
  $box->();

  # S0
  $_ -= $ics + $d for $y1, $y2;
  $x1 = $ics;
  $x2 = $x1 + $d;
  $box->();

  # S1
  $_ += $ics + $d for $x1, $x2;
  $box->();

  # S2
  $_ += $ics + $d for $x1, $x2;
  $box->();

  # S3
  $_ += $ics + $d for $x1, $x2;
  $box->();

  # S4
  $_ += $ics + $d for $x1, $x2;
  $box->();

  # S5
  $_ += $ics + $d for $x1, $x2;
  $box->();

  pgsave();

  pgsch(3);

  pgsvp(0,1,0,1);
  pgswin(0,1,0,1);
  pgptxt(0.5, 0.8, 0, 0.5, '\fr89 - 2008, 0.4 keV');

  pgunsa();

}
