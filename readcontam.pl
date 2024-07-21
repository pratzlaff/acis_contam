#! /usr/bin/perl -w
use strict;

use PDL;
use PDL::Graphics::PGPLOT;
use PDL::Graphics::LUT;
use Astro::FITS::CFITSIO;
use PGPLOT;
use Chandra::Tools::Common qw( read_bintbl_cols check_status );

my %opts = (
	    separateimages => 1,
	    combinedimages => 1,
	    imagesamp => 8,
	    fakeimages => 0,
	    odimage => 1, # images of optical depth instead of transmission
	    ctable => 'heat',
	    cidev => '89_31_cdiff.png/png',
	    images_energy => 0.3, # energy at which to compute above
	    itf => 0, # 0=linear, 1=log, 2=sqrt
	   );

my ($dev, $win);

#my $f04 = '/usr/local/caldb-4.1.2/data/chandra/acis/contam/acisD1999-08-13contamN0004.fits';

my $f04 = '/data/legs/rpete/contam/acisD1999-08-13contamN0004.fits';
my $f89 = '/data/legs/rpete/contam/acisD1999-08-13contamN9989.fits';
my $f31 = '/data/legs/rpete/contam/acisD1999-08-13contamN9931.fits';

my %p04 = read_contamfile($f04);
my %p89 = read_contamfile($f89);
my %p31 = read_contamfile($f31);

my $year = sequence(1000)/1000 * 10 + 2000;
my ($x, $y, $e) = (512, 512, 0.5);
my $depth = depth_vs_time($p31{'ACIS-7'}, $x, $y, $e, $year);
$dev = '31_depth_vs_time.ps/ps';
dev $dev;
line $year, $depth,
  {
   xtitle => 'Year',
   ytitle => 'optical depth',
   title => "31 : ACIS-S3, CHIPX = $x, CHIPY = $y, $e keV",
   };
exit;

#i3aimcheck(\%p89, \%p31);
#exit;

=begin comment

$dev = '/xs';
$win = PDL::Graphics::PGPLOT::Window->new(Device => $dev,
					       nxpanel => 2,
					       nypanel => 2,
					       aspectratio => 2/2,
					      );
my $m = transmission_map($p89{'ACIS-3'}, 1, 1024, 1, 1024, 2008, 0.4);
my $max = $m->max;
my $tr = $win->transform( $m->dims,
			  { RefPos => [[0,0], [1, 1]],
			    PixInc => [1,1],
			  });
#my $title = "$ps : $s $year, $e keV";
$win->imag($m, { Transform => $tr, range => [0,$max], drawwedge => 1,
#		 title => $title,
		 ytitle => 'CHIPY',
		 xtitle => 'CHIPX',
	       }
	  );
exit;

=cut

# print  transmission($p89{'ACIS-7'}, 1, 1, 2008, pdl(0.3,0.4)),"\n";
# print  transmission($p89{'ACIS-7'}, 512, 512, 2008, pdl(0.3, 0.4)),"\n";
# #print  transmission($p89{'ACIS-7'}, long(1,512,1,512), long(1,512,1,512), 2008, pdl(0.3,0.4)),"\n";
# #exit;

# my ($x, $y) = (sequence(1e1) % 512 +1, sequence(1e1)% 512+1);
# my @t1;
# for (0..$x->nelem-1) {
#   push @t1, transmission($p89{'ACIS-7'}, $x->at($_), $y->at($_), 2008, 0.3)->at(0);
# }
# my $t1 = pdl \@t1;
# my $t2 = transmission($p89{'ACIS-7'}, $x, $y, 2008, pdl(0.3,0.3));
# print join(', ', stats($t2 - $t1)),"\n";
# print $t1,"\n";
# print $t2,"\n";
# exit;


my ($elo, $ehi) = (0.31,2);
$_ = log $_ for $elo, $ehi;
my $egrid = exp(sequence(1000) / 1000 * ($ehi-$elo) + $elo);
$_ = exp $_ for $elo, $ehi;

my %chips = (
	     'ACIS-0' => 'ACIS-I0',
	     'ACIS-1' => 'ACIS-I1',
	     'ACIS-2' => 'ACIS-I2',
	     'ACIS-3' => 'ACIS-I3',
	     'ACIS-4' => 'ACIS-S0',
	     'ACIS-5' => 'ACIS-S1',
	     'ACIS-6' => 'ACIS-S2',
	     'ACIS-7' => 'ACIS-S3',
	     'ACIS-8' => 'ACIS-S4',
	     'ACIS-9' => 'ACIS-S5',
	     );

#my @years = (2002, 2004, 2006, 2008);
# E0102 observation dates
my @years = (
2000.25727842936,
2003.59340731671,
2006.19863593764,
2009.0474530666,
);
# aimpoints from http://cxc.harvard.edu/proposer/POG/html/chap6.html#tth_sEc6.10
my @config = (
	      ['ACIS-3', 970, 967], # I3 aimpoint
#	      ['ACIS-7', 221, 520], # S3 aimpoint
#	      ['ACIS-5', 512, 150], # center of S1
	      );

$dev = 'e0102dates_04_89_31.ps/ps';
$dev = '/xs';
$win = PDL::Graphics::PGPLOT::Window->new(device => $dev,
					     nxpanel => 3,
					     nypanel => scalar(@config),
#					     aspectratio => 3/2,
					     );
$win->autolog(1);

for (@config) {
  my ($key, $x, $y) = @$_;
  my $s = $chips{$key};

  for ([\%p04, '04'], [\%p89, '89'], [\%p31, '31']) {
    my ($p, $ps) = @$_;

    my $ls = 1;
    for my $year (@years) {
      $win->line($egrid, transmission($p->{$key}, $x, $y, $year, $egrid),
	{ title => "$ps $s: CHIPX = $x, CHIPY = $y",
	  axis => 'logx',
	  yrange => [0,1],
	  linestyle => $ls,
	});
      $win->hold;
      ++$ls;
    }
    $win->release;
  }
}
$win->close;

$dev = 'e0102dates_04_89_31_ratio.ps/ps';
$win = PDL::Graphics::PGPLOT::Window->new(device => $dev,
					     nxpanel => 2,
					     nypanel => scalar(@config),
#					     aspectratio => 2/3,
					     );
$win->autolog(1);

for (@config) {
  my ($key, $x, $y) = @$_;
  my $s = $chips{$key};

  my %t1 = map { $_, transmission($p04{$key}, $x, $y, $_, $egrid) } @years;

  for ([\%p89, '89'], [\%p31, '31']) {

    my ($p, $ps) = @$_;

    my $ls = 1;

    for my $year (@years) {

      my $t2 = transmission($p->{$key}, $x, $y, $year, $egrid);

      $win->line($egrid, $t1{$year} / $t2,
	{ title => "04 / $ps: $s CHIPX = $x, CHIPY = $y",
	  axis => 'logx',
	  yrange => [0.35,2],
	  linestyle => $ls,
	});
      $win->hold;
      ++$ls;
    }
    $win->release;
  }
}
$win->close;

my %images;
my $e = $opts{images_energy};

if ($opts{separateimages} or $opts{combinedimages}) {

  my ($max, $min);
  for ([\%p89, '89'], [\%p31, '31']) {
    my ($p, $ps) = @$_;
    for my $year (@years) {
      for my $key (keys %$p) {
	my $m;
	if ($opts{fakeimages}) {
	  $m = sequence(10,10);
	}
	else {
	  my ($t,$d) = transmission_map_multiplepos($p->{$key}, 1, 1024, 1, 1024, $year, $e);
	  $d = (transmission_map_multiplepos($p->{$key}, 1, 1024, 1, 1024, $year, 0.3))[1] - (transmission_map_multiplepos($p->{$key}, 1, 1024, 1, 1024, $year, 0.27))[1];
	  $m = $opts{odimage} ? $d : $t;
	}
	$images{$ps}{$key}{$year} = $m;
	my $mmax = $m->max;
	$max = $mmax unless defined $max and $max > $mmax;
	my $mmin = $m->min;
	$min = $mmin unless defined $min and $min < $mmin;
      }
    }
  }
  $images{max} = $max;
  $images{min} = $min;
}

if ($opts{separateimages}) {

  for my $key (keys %p89) {
    my $s = $chips{$key};

    my $dev = $s.'_image.ps/ps';
    my $win = PDL::Graphics::PGPLOT::Window->new(Device => $dev,
						 nxpanel => 2,
						 nypanel => 2,
						 aspectratio => 2/2,
						);
    for ([\%p89, '89'], [\%p31, '31']) {
      my ($p, $ps) = @$_;
      my ($max, $min);

      for my $year (@years) {
	my $m = $images{$ps}{$key}{$year};
#	my $m = transmission_map_multiplepos($p->{$key}, 1, 1024, 1, 1024, $year, $e);
	$max = $m->max unless defined $max;
	$min = $m->min unless defined $min;
	my $tr = $win->transform( $m->dims,
				  { RefPos => [[0,0], [1, 1]],
				    PixInc => [1,1],
				  });
	my $title = "$ps : $s $year, $e keV";
	$win->imag($m, { Transform => $tr, range => [$min,$max], drawwedge => 1,
			 title => $title,
			 ytitle => 'CHIPY',
			 xtitle => 'CHIPX',
		       }
		  );
      }
    }
    $win->close;
  }
}

if ($opts{combinedimages}) {

  dev $opts{cidev}, 2, 4, { aspectratio => 4/2 };

  if ($opts{cidev} =~ /png$/) {
    pgpap(10, 2);
    pgsubp(2, 4);
  }
  set_ctab($opts{ctable});
  pgsitf($opts{itf});

  for my $y (@years) {
    for my $p ('89', '31') {
      pgpage();
      acis_layout($images{min}, $images{max}, "$p : $y \\gt\\d0.30\\u - \\gt\\d0.27\\u", #"$p : $y $e keV",
		  map { $images{$p}{'ACIS-'.$_}{$y} } (0..9)
		  );
    }
  }

}

#line $egrid, transmission($p{'ACIS-7'}, 512, 512, 2009, $egrid), { title => 'ACIS-S3: CHIPX = CHIPY = 512', axis => 'logx', };

exit;

sub transmission_map {
  my ($p, $xlo, $xhi, $ylo, $yhi, $year, $e) = @_;

  my $nx = $xhi-$xlo+1;
  my $ny = $yhi-$ylo+1;

  $_ /= $opts{imagesamp} for $nx, $ny;

  my $map = zeroes($nx, $ny);

  my $xvals = xvals($map)*$opts{imagesamp} + $opts{imagesamp}/2 + $xlo;
  my $yvals = yvals($map)*$opts{imagesamp} + $opts{imagesamp}/2 + $ylo;

  $_ = $_->flat for $xvals, $yvals;

  for my $i (0..$xvals->nelem-1) {
    my ($x, $y) = ($xvals->at($i), $yvals->at($i));
    $map->set($x-$xlo, $y-$ylo,
	      transmission($p, $x, $y, $year, $e)
	      );
  }
  return wantarray ? ($map, -log($map)) : $map;
}

sub transmission_map_multiplepos {
  my ($p, $xlo, $xhi, $ylo, $yhi, $year, $e) = @_;

  my $nx = $xhi-$xlo+1;
  my $ny = $yhi-$ylo+1;

  $_ /= $opts{imagesamp} for $nx, $ny;

  my $map = zeroes($nx, $ny);

  my $xvals = xvals($map)*$opts{imagesamp} + $opts{imagesamp}/2 + $xlo;
  my $yvals = yvals($map)*$opts{imagesamp} + $opts{imagesamp}/2 + $ylo;

#  for my $i (0..$xvals->getdim(0)-1) {
#    my ($x, $y) = ($xvals->slice("($i)"), $yvals->slice("($i)"));
#    (my $tmp = $map->slice("($i)")) .= 
#	      transmission_multiplepos($p, $x, $y, $year, $e);
  my $n = 4 * $opts{imagesamp};
  for my $i (0..$xvals->getdim(0)/$n-1) {
    my $ss = ($n*$i) . ':' . ($n*$i+$n-1);
    my ($x, $y) = ($xvals->slice($ss)->flat, $yvals->slice($ss)->flat);
    (my $tmp = $map->slice($ss)->flat) .=
	      transmission_multiplepos($p, $x, $y, $year, $e);
  }

  return wantarray ? ($map, -log($map)) : $map;

}

sub transmission {
  my ($p, $x, $y, $year, $e) = @_;
  $e = PDL::Core::topdl($e);

  my $tau = tdep_tau($p, $x, $y, $year);

  # dims 2 x edims
  my $mu = interpol($e->dummy(0,$p->{energy}->getdim(1)), $p->{energy}, $p->{mu});

  # FIXME : not dealing with distinct j
  return exp(-sumover($mu * $tau));
}

sub transmission_multiplepos {
  my ($p, $x, $y, $year, $e) = @_;
  $e = PDL::Core::topdl($e);

  my $tau = tdep_tau_multiplepos($p, $x, $y, $year);

  # dims 2 x edims
  my $mu = interpol($e->dummy(0,$p->{energy}->getdim(1)), $p->{energy}, $p->{mu});

  # FIXME : not dealing with distinct j
  return exp(-sumover($mu * $tau));
}

sub tdep_tau_multiplepos {
  my ($p, $x, $y, $year) = @_;

  my $time = ($year - 1998) * 365.24 * 86400;

  my $i = 1024 * ($y-1) + ($x-1);

  return
    interpol($time, $p->{time}, $p->{tau0}) +
    interpol($time, $p->{time}, $p->{tau1}) * $p->{fxy}->dice($i)->xchg(0,1);
}

sub tdep_tau {
  my ($p, $x, $y, $year) = @_;

  my $time = ($year - 1998) * 365.24 * 86400;

  my $i = 1024 * ($y-1) + ($x-1);

  return
    interpol($time, $p->{time}, $p->{tau0}) +
    interpol($time, $p->{time}, $p->{tau1}) * $p->{fxy}->index($i);
}


sub read_contamfile {

  my $f = shift;

  my $s=0;
  my $fptr = Astro::FITS::CFITSIO::open_file($f,Astro::FITS::CFITSIO::READONLY(),$s);
  check_status($s) or die;

  my %p;

  for my $n (2..11) {
    $fptr->movabs_hdu($n, undef, $s);
    check_status($s) or die "could not move to HDU number $n";

    my $hdr = Astro::FITS::CFITSIO::fits_read_header($fptr);
    exists $hdr->{DETNAM} or die "no DETNAM keyword in HDU number $n";

    my $detnam = $hdr->{DETNAM};
    ($detnam) =  $detnam =~ /^'(\S+)/;
    my @cols = qw( component weight kappa n_energy energy mu n_time time tau0 tau1 fxy );
    my @d = read_bintbl_cols($fptr, @cols);

    @{$p{$detnam}}{@cols} = @d;

    # FIXME : only uses the first row values of n_energy and n_time

    $p{$detnam}{$_} = $p{$detnam}{$_}->slice('0:'.($p{$detnam}{n_energy}->at(0)-1).',')->sever for qw( energy mu );
    die if which($p{$detnam}{n_energy} != $p{$detnam}{n_energy}->at(0))->nelem;

    $p{$detnam}{$_} = $p{$detnam}{$_}->slice('0:'.($p{$detnam}{n_time}->at(0)-1).',')->sever for qw( tau0 tau1 time );
    die if which($p{$detnam}{n_time} != $p{$detnam}{n_time}->at(0))->nelem;

    # FIXME : only dealing with a single value of j
    die if which($p{$detnam}{component} != $p{$detnam}{component}->at(0))->nelem;
  }

  return %p;
}

sub rotccw {
  my $img = shift;
  return $img->xchg(0,1)->slice('-1:0');
}

sub rotcw {
  my $img = shift;
  return $img->xchg(0,1)->slice(',-1:0');
}

sub acis_layout {
  my ($l, $h, $text, @img) = @_;

    # viewport size of each chip
  my $d = 0.14;
  my $ics = (1-(6*$d))/7; # inter-chip spacing
  my $vm = 0.5 * (1- 3*$d + 2*$ics); # vertical margins

  my ($x1, $x2, $y1, $y2);

  my $box = sub {
    my $img = shift;
    pgsvp($x1, $x2, $y1, $y2);
    pgswin(0.5, 1024.5, 0.5, 1024.5);
    pgimag($img->float->get_dataref,
	   $img->getdim(0),
	   $img->getdim(1),
	   1, $img->getdim(0),
	   1, $img->getdim(1),
	   $l, $h,
	   [
	    $opts{imagesamp}/2 + 0.5, $opts{imagesamp}, 0,
	    $opts{imagesamp}/2 + 0.5, 0, $opts{imagesamp},
	   ],
	   );
  };

  # I0
  $y2 = 1-$vm;
  $y1 = $y2-$d;
  $x1 = 3*$ics + 2*$d;
  $x2 = $x1+$d;
  $box->(rotcw($img[0]));

  # I1
  $_ += $ics + $d for $x1, $x2;
  $box->(rotccw($img[1]));

  # I2
  $_ -= $ics + $d for $y1, $y2, $x1, $x2;
  $box->(rotcw($img[2]));

  # I3
  $_ += $ics + $d for $x1, $x2;
  $box->(rotccw($img[3]));

  # S0
  $_ -= $ics + $d for $y1, $y2;
  $x1 = $ics;
  $x2 = $x1 + $d;
  $box->($img[4]);

  # S1
  $_ += $ics + $d for $x1, $x2;
  $box->($img[5]);

  # S2
  $_ += $ics + $d for $x1, $x2;
  $box->($img[6]);

  # S3
  $_ += $ics + $d for $x1, $x2;
  $box->($img[7]);

  # S4
  $_ += $ics + $d for $x1, $x2;
  $box->($img[8]);

  # S5
  $_ += $ics + $d for $x1, $x2;
  $box->($img[9]);

  pgsave();

  pgsch(3);

  pgsvp(0,1,0,1);
  pgswin(0,1,0,1);
  pgptxt(0.5, 0.8, 0, 0.5, $text);

  pgsvp(0.25, 0.75, 0.4, 0.75);
  pgsch(7);
  pgwedg('BI', 1.5, 1, $l, $h, '');

  pgunsa();
}

sub set_ctab {
  my $name = shift;
  my ($l, $r, $g, $b) = lut_data($name);
  pgctab(
         $l->float->get_dataref,
         $r->float->get_dataref,
         $g->float->get_dataref,
         $b->float->get_dataref,
         $l->nelem,
         1, 0.5);
}

sub i3aimcheck {
  my ($p89, $p31) = @_;
  my $e = 0.5;

  $_ = $_->{'ACIS-3'} for $p89, $p31;

  my $p = $p31;

  my ($k, $x, $y) =  ('ACIS-3', 970, 967); # I3 aimpoint

  my $fxyi = 1024 * ($y-1) + ($x-1);

  print "components = ". $p->{component}."\n";
  print "kappa = ". $p->{kappa}."\n";


  my @y = (2002, 2004, 2006, 2008);
  for my $y (@y) {
    my $exp = 0;
    my $time = ($y - 1998) * 365.24 * 86400;
    for my $i (0..$p->{component}->nelem-1) {
      my $c = $p->{component}->at($i);
      my $ss = ",($i)";
      my $tau0 = interpol($time, $p->{time}->slice($ss), $p->{tau0}->slice($ss));
      my $tau1 = interpol($time, $p->{time}->slice($ss), $p->{tau1}->slice($ss));
      my $fxy = $p->{fxy}->slice($ss)->index($fxyi);
      my $mu = interpol($e, $p->{energy}->slice($ss), $p->{mu}->slice($ss));

      printf "$y: component $c - tau0 = %.2f, tau1 = %.2f, fxy = %.2f, mu = %.2f\n", $tau0, $tau1, $fxy, $mu;
      $exp += $mu * ( $tau0 + $tau1 * $fxy );
    }
    my $trans = exp(-$exp);
    printf "$y: tau = %.2f, trans = %.2f\n", $exp, $trans;
  }
}

sub depth_vs_time {
  my ($p, $x, $y, $e, $year) = @_;

  my @depth;

  for (0..$year->nelem-1) {
    push @depth, -log(transmission($p, $x, $y, $year->at($_), $e)->at(0));
  }
  return pdl \@depth;

}
