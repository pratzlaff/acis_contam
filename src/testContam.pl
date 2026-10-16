use PDL;
use PDL::Graphics::PGPLOT::Window;
use Contam;

my ($f, $ccd_id, $x, $y, $e, $year) = (
    "/data/legs/rpete/flight/acis_contam/data/acisD1999-08-13contamN0016.fits",
    7, 220, 520, 0.3, 2012,
  );

my $sample=8;
my $taus = Contam::image($f, $ccd_id, $e, $year, $sample);
my $tau = $taus->sumover();
printf "min=%.4g, max=%.4g, mean=%.4g\n", $tau->min, $tau->max, $tau->avg;
my $win1 = PDL::Graphics::PGPLOT::Window->new();
my $tr = $win1->transform( $tau->dims,
			  { RefPos => [[0,0], [$sample/2+0.5, $sample/2+0.5]],
			    PixInc => [$sample, $sample],
			  });
$win1->imag($tau, { transform => $tr, drawwedge => 1 });
$win1->label_axes('CHIPX', 'CHIPY', "S3 Contamination Layer Optical Depth at $e keV, $year");

my  ($energy, $taus) = Contam::energy($f, $ccd_id, $x, $y, $year);
$tau = $taus->sumover();
printf "min=%.4g, max=%.4g, mean=%.4g\n", $tau->min, $tau->max, $tau->avg;
my $win2 = PDL::Graphics::PGPLOT::Window->new();
$win2->line($energy->log10, $tau);
$win2->label_axes("log\\d10\\u Energy (keV)", '\gt', "S3 Contamination Layer Optical Depth at CHIPXY=($x,$y), $year");

my  ($dates, $taus) = Contam::year($f, $ccd_id, $x, $y, $e);
$tau = $taus->sumover();
printf "min=%.4g, max=%.4g, mean=%.4g\n", $tau->min, $tau->max, $tau->avg;
my $win3 = PDL::Graphics::PGPLOT::Window->new();
$win3->line($dates, $tau);
$win3->label_axes("Year", '\gt', "S3 Contamination Layer Optical Depth at CHIPXY=($x,$y), $e keV");

