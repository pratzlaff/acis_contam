#! /usr/bin/perl -w
use strict;

# on solaris, run with PERL5LIB=/home/rpete/local/perlmods

my $version = '0.1';

use FindBin;
use Config;
use Carp;
use PDL;
use PDL::Graphics::PGPLOT;
use PGPLOT;

use Getopt::Long;
my %default_opts = (
		    exec => 1,
		    obsfile => 'acis_evt2.fits',
		    detsubsys => 'ACIS-I3;QE=1;UNIFORM',
		    xpos => 970,
		    ypos => 967,
		    emin => 0.31,
		    emax => 2,
		    n => 100,
		    dev => '/xs',
		    );
# for each of 31, 89
#   for each year, get the transmission at I3 970,967
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'exec!', 'dev=s', 'obsfile=s', 'detsubsys=s',
	   'xpos=i', 'ypos=i', 'emin=f', 'emax=f', 'n=i',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();


my ($detsubsys, $xpos, $ypos, $emin, $emax, $n) =
  @opts{ qw( detsubsys xpos ypos emin emax n ) };

my @years =  (
	      2000.25727842936,
	      2003.59340731671,
	      2006.19863593764,
	      2009.0474530666,
	     );

my @t = map { ($_ - 1998) * 86400 * 365.24 } @years;
my @mjd = map { 50814 + $_/86400 } @t;


my %f = (
	 '04' => 'acisD1999-08-13contamN0004.fits',
	 '89' => 'acisD1999-08-13contamN9989.fits',
	 '31' => 'acisD1999-08-13contamN9931.fits',
	 '05' => 'acisD1999-08-13contamN0005.fits',
	 '32' => 'acisD1999-08-13contamN9932.fits',
	 );

my %trans;

for my $k (keys %f) {

  run_command('punlearn', 'ardlib');
  run_command('pset', 'ardlib',
	      map {
		'AXAF_ACIS' . $_ . '_CONTAM_FILE=' .
		$f{$k} . '[AXAF_CONTAM'.($_+1).']'
	      } 0..9
	     );

  for my $i (0..$#t) {

    my $out = `./ardlib_qe '$detsubsys;TIME=$t[$i]' $opts{obsfile} $xpos $ypos $emin $emax $n`;
    my (@e, @trans);
    for (split "\n", $out) {
      my ($e, $trans) = split ' ', $_;
      push @e, $e;
      push @trans, $trans;
    }
    $trans{$k}{$years[$i]}{e} = pdl \@e;
    $trans{$k}{$years[$i]}{trans} = pdl \@trans;

  }
}

my @k = ('89', '31');
@k = ('32');
#my $dev = 'e0101dates_ardlib_ratios.ps/ps';
#$dev = '/xs';
dev $opts{dev}, scalar(@k), 1;
autolog(1);
my $num = '05';
for my $k (@k) {
  my @ls = (1..@years);
  my $y = $years[0];
  line($trans{$num}{$y}{e},
       $trans{$num}{$y}{trans} / $trans{$k}{$y}{trans},
       {
	border => 1,
	yrange => [0.35, 2],
	axis => 'logx',
	title => "ardlib I3 : (x,y)=($xpos,$ypos) : $num / $k",
	xtitle =>'Energy',
	ytitle => 'transmittance ratio',
       }
      );
  hold;
  for my $i (1..$#years) {
    my $y = $years[$i];
    line($trans{$num}{$y}{e},
	 $trans{$num}{$y}{trans} / $trans{$k}{$y}{trans},
	 {
	  linestyle => $ls[$i],
	 }
	);
  }

  legend(
	 [ map { sprintf "%.1f", $years[$_] } 0..$#years ],
	 0, 1.6, { LineStyle => [ 1..@years ] },
	);
  release;
}

sub run_command {
  my ($cmd, @args) = @_;
  print STDERR join(' ', $cmd, map ( "'$_'", @args)),"\n\n";
  if ($opts{'exec'}) {
    system($cmd, @args);
    die unless ($opts{errignore} or $? == 0);
  }
}

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
}

