use PDL;
use PDL::Graphics::PGPLOT::Window;
use Contam;

my $f = "/usr/local/ciao/CALDB/data/chandra/acis/contam/acisD1999-08-13contamN0016.fits";

my ($ccd_id, $e, $x, $y, $year)  = (7, 0.3, 220, 520, 2012);

my $taus = Contam::image($f, $ccd_id, $e, $year);
my $tau = $taus->sumover();
print $tau->min, " ", $tau->max, "\n";
my $win1 = PDL::Graphics::PGPLOT::Window->new();
$win1->imag($tau);
$win1->label_axes('CHIPX', 'CHIPY', "S3 Contamination Layer Optical Depth at $e keV, $year");

my  ($energy, $taus) = Contam::energy("/usr/local/ciao/CALDB/data/chandra/acis/contam/acisD1999-08-13contamN0016.fits", $ccd_id, $x, $y, $year);
$tau = $taus->sumover();
print $tau->min, " ", $tau->max, "\n";
my $win2 = PDL::Graphics::PGPLOT::Window->new();
$win2->line($energy->log10, $tau);
$win2->label_axes("log\\d10\\u Energy (keV)", '\gt', "S3 Contamination Layer Optical Depth at CHIPXY=($x,$y), $year");

my  ($dates, $taus) = Contam::year("/usr/local/ciao/CALDB/data/chandra/acis/contam/acisD1999-08-13contamN0016.fits", $ccd_id, $x, $y, $e);
$tau = $taus->sumover();
print $tau->min, " ", $tau->max, "\n";
my $win3 = PDL::Graphics::PGPLOT::Window->new();
$win3->line($dates, $tau);
$win3->label_axes("Year", '\gt', "S3 Contamination Layer Optical Depth at CHIPXY=($x,$y), $e keV");

