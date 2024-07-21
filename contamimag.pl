#! /usr/bin/perl -w
use strict;

#
# In addition to the CPAN modules PDL, PGPLOT and
# Astro::FITS::CFITSIO, this program also requires the custom module
# Chandra::Tools::Common which is available from Pete Ratzlaff
# <pratzlaff@cfa.harvard.edu>
#
# on the HEAD network, the messy details of the modules can usually be
# ignored on Linux machines by running
#
#   env PGPLOT_DIR=/data/legs/rpete/linux-x86_64/pgplot /proj/axaf/bin/perl /path/to/contamimag
#


=head1 NAME

contamimag - display images of ACIS contamination layer optical depth and transmittance

=head1 SYNOPSIS

contamimag [options]

=head1 DESCRIPTION

Display images of ACIS contamination layer optical depth and
transmittance. The contamination calibration file used by default is
the latest found in the directory name stored in the CALDB environment
variable. Thus, typically this program would be run from within a CIAO
environment.

=head1 OPTIONS

=over 4

=item --help

Show help and exit.

=item --version

Show version and exit.

=item --energy=f

The energy, in keV, for which the optical depth is calculated. The default is 1 keV.

=item --year=f

The time at which the optical depth is calculated. The default is 2012.

=item --trans

Plot the contamination layer transmittance rather than optical depth.

=item --det=i0|i1|i2|i3|s0|s1|s2|s3|s4|s5

Plot only a single detector, rather than all of them.

=item --file=s

Specify a calibration contamination file, rather than choosing the
latest from the CALDB directory.

=item --table=s

Specify the color table name. The default is I<heat>.

=item --itf=i

The image transfer function. Possible values are 0 (default), 1
and 2 for linear, log and sqrt, respectively.

=item --dev=s

PGPLOT device. The default is I</xs>. Examples are I<file.ps/cps>
(color Postscript file) and I<file.png/png> (PNG file).

==item --cinvert, --nocinvert

Explicitly set foreground (text) and background colors to black and white,
respectively.  This is enabled by default if the plotting device is
either the PNG or GIF format.

=back

=head1 AUTHOR

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> March 2012

=head1 SEE ALSO

perl(1).

=cut

use PDL;
use PDL::Graphics::PGPLOT;
use PDL::Graphics::LUT;
use PGPLOT;
use Astro::FITS::CFITSIO;
use Chandra::Tools::Common qw( read_bintbl_cols check_status );

my $version = '0.1';

use FindBin;
use Config;
use Carp;
use Getopt::Long;

my %default_opts = (
		    table => 'heat',
		    itf => 0,		  # 0=linear, 1=log, 2=sqrt
		    dev => '/xs',
		    year => 2012,
		    energy => 1.0,
#		    file => '/usr/local/ciao/CALDB/data/chandra/acis/contam/acisD1999-08-13contamN0006.fits',
		   );
my %opts = %default_opts;

GetOptions(\%opts,
	   'help!', 'version!', 'debug!', 'trans!',
	   'dev=s', 'year=f', 'n=i', 'det=s', 'energy=f',
	   'table=s', 'itf=i', 'file=s', 'cinvert!',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

@ARGV and die "Usage: $0 [options]\nTry --help for more information\n";

$opts{cinvert} = 1 if ($opts{dev} =~ /(png|gif)$/i and ! exists $opts{cinvert});

my $e = $opts{energy};

$e > 0 && $e < 20 or die "unsupported energy = $e keV\n";

for ($opts{itf}) {
  last if $_==0 or $_==1 or $_==2;
  die "unsupported ITF=$opts{itf}\n";
}

my ($dev, $win);

if (! exists $opts{file}) {
  if (! exists $ENV{CALDB}) {
    die "--file unspecified and there is no CALDB environment variable\n";
  }
  my $dir = $ENV{CALDB}.'/data/chandra/acis/contam';
  opendir DIR, $dir or die "--file unspecified and could not opendir $dir: $!\n";
  my @f = sort grep /^acis.*contamN\d+\.fits$/, readdir DIR;
  closedir DIR;
  @f or die "--file unspecified and no contamination were files found in $dir\n";
  $opts{file} = $dir. '/' . $f[-1];
}

my $f = $opts{file};

my %p = read_contamfile($f);
my @det = qw/ i0 i1 i2 i3 s0 s1 s2 s3 s4 s5 /;
my %det;
@det{@det} = 0..$#det;

my $desc = $opts{trans} ? 'contam layer transmittance' : 'contam layer optical depth';
$desc .= " at $e keV";
$desc .= sprintf(', year=%.10g', $opts{year});

if ($opts{det}) {

  my $det = lc($opts{det});
  exists $det{$det} or die "detector '$det' is unrecognized\n";

  my $dev = '/xs';
  my $win = PDL::Graphics::PGPLOT::Window->new(device => $opts{dev},
					       justify => 1,
					      );
  $win->ctab(lut_data($opts{table}));
  #my $m = transmittance_map($p{'ACIS-'.$det{$det}}, 1, 1024, 1, 1024, 2008, $e);

  my ($trans,$tau) = transmittance_map_multiplepos($p{'ACIS-'.$det{$det}}, 1, 1024, 1, 1024, $opts{year}, $e);
  my $m = $opts{trans} ? $trans : $tau;

  my $title = uc($det).' - '.$desc;

  my $tr = $win->transform( $m->dims,
			    { RefPos => [[0,0], [1, 1]],
			      PixInc => [1, 1],
			    });
  if ($opts{cinvert}) {
    pgscr(0, 1, 1, 1);
    pgscr(1, 0, 0, 0);
  }

  $win->imag($m, { Transform => $tr, range => [$m->minmax], drawwedge => 1,
		   title => $title,
		   ytitle => 'CHIPY',
		   xtitle => 'CHIPX',
		 }
	    );
}
else {

  dev $opts{dev};
  set_ctab($opts{table});

  my @images;
  my ($min, $max);
  for my $i (0..$#det) {
    my ($trans,$tau) = transmittance_map_multiplepos($p{'ACIS-'.$i}, 1, 1024, 1, 1024, $opts{year}, $e);
    my $m = $opts{trans} ? $trans : $tau;
    push @images, $m;

    my $mmax = $m->max;
    $max = $mmax unless defined $max and $max > $mmax;
    my $mmin = $m->min;
    $min = $mmin unless defined $min and $min < $mmin;

  }

  pgpap(8.5, 1);

  if ($opts{cinvert}) {
    pgscr(0, 1, 1, 1);
    pgscr(1, 0, 0, 0);
  }

  pgsitf($opts{itf});

  acis_layout($min, $max, $desc, @images);
}

exit;

sub transmittance_map {
  die "transmittance_map() is too slow, transmittance_map_multipos() is what you want";
  my ($p, $xlo, $xhi, $ylo, $yhi, $year, $e) = @_;

  my $nx = $xhi-$xlo+1;
  my $ny = $yhi-$ylo+1;

  my $map = zeroes($nx, $ny);

  my $xvals = xvals($map) + $xlo;
  my $yvals = yvals($map) + $ylo;

  $_ = $_->flat for $xvals, $yvals;

  for my $i (0..$xvals->nelem-1) {
    my ($x, $y) = ($xvals->at($i), $yvals->at($i));
    $map->set($x-$xlo, $y-$ylo,
	      transmittance($p, $x, $y, $year, $e)
	      );
  }
  return wantarray ? ($map, -log($map)) : $map;
}

sub transmittance_map_multiplepos {
  my ($p, $xlo, $xhi, $ylo, $yhi, $year, $e) = @_;

  my $nx = $xhi-$xlo+1;
  my $ny = $yhi-$ylo+1;

  my $map = zeroes($nx, $ny);

  my $xvals = xvals($map) + $xlo;
  my $yvals = yvals($map) + $ylo;

#  for my $i (0..$xvals->getdim(0)-1) {
#    my ($x, $y) = ($xvals->slice("($i)"), $yvals->slice("($i)"));
#    (my $tmp = $map->slice("($i)")) .= 
#	      transmittance_multiplepos($p, $x, $y, $year, $e);
  my $n = 1;
  for my $i (0..$xvals->getdim(0)/$n-1) {
    my $ss = ($n*$i) . ':' . ($n*$i+$n-1);
    my ($x, $y) = ($xvals->slice($ss)->flat, $yvals->slice($ss)->flat);
    (my $tmp = $map->slice($ss)->flat) .=
	      transmittance_multiplepos($p, $x, $y, $year, $e);
  }

  return wantarray ? ($map, -log($map)) : $map;
}

sub transmittance {
  my ($p, $x, $y, $year, $e) = @_;
  $e = PDL::Core::topdl($e);

  my $tau = tdep_tau($p, $x, $y, $year);

  # dims 2 x edims
  my $mu = interpol($e->dummy(0,$p->{energy}->getdim(1)), $p->{energy}, $p->{mu});

  # FIXME : not dealing with distinct j
  return exp(-sumover($mu * $tau));
}

sub transmittance_multiplepos {
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

  print STDERR "Reading $f\n";

  my $s=0;
  my $fptr = Astro::FITS::CFITSIO::open_file($f,Astro::FITS::CFITSIO::READONLY(),$s);
  check_status($s) or die;

  my %p;

  for my $n (2..11) {
    $fptr->movabs_hdu($n, undef, $s);
    check_status($s) or die "could not move to HDU number $n";

    my $hdr = Astro::FITS::CFITSIO::fits_read_header($fptr);
    exists $hdr->{CCD_ID} or die "no CCD_ID keyword in HDU number $n";

    my $id = $hdr->{CCD_ID};
    my $detnam = 'ACIS-' . $id;

# this was for when detnam was taken from the header DETNAM entry
#    ($detnam) =  $detnam =~ /^'(\S+)/;

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
	   [0.5, 1, 0, 0.5, 0, 1],
#	   [
#	    $opts{sample}/2 + 0.5, $opts{sample}, 0,
#	    $opts{sample}/2 + 0.5, 0, $opts{sample},
#	   ],
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

  pgsch(1.2);

  pgsvp(0,1,0,1);
  pgswin(0,1,0,1);
  pgptxt(0.5, 0.73, 0, 0.5, $text);

  pgsvp(0.25, 0.75, 0.4, 0.75);
  pgsch(3);
  pgwedg('BI', 2.7, 1, $l, $h, '');

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

sub set_ctab_pdl {
  my $name = shift;
  my ($l, $r, $g, $b) = lut_data($name);
  ctab(
       $l, $r, $g, $b,
#       1, 0.5,
       );
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
    push @depth, -log(transmittance($p, $x, $y, $year->at($_), $e)->at(0));
  }
  return pdl \@depth;

}

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
}
