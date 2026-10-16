import argparse
import astropy.io.fits
import matplotlib.pyplot as plt
import numpy as np

# https://cxc.cfa.harvard.edu/ciao/download/papers/spatial_contam_effects.pdf

def year2time(year):
    return 365.2425*86400*(year-1998)

def testfxy(fxy):
    n = fxy.shape[0]
    for i in range(1, n):
        diff = fxy[0] - fxy[i]
        min, max = diff.min(), diff.max()
        print(f'fxy[0] - fxy[{i}]: {min=:.4g}, max={max=:.4g}')

def image(fname, ccd_id, energy, year, sample=1):
    with astropy.io.fits.open(fname) as hdulist:
        hdu = hdulist[ccd_id+1]
        data = hdu.data
        n = data['component'].shape[0]
        taus = np.zeros_like(data['fxy'][:,::sample,::sample])
        time = year2time(year)
        for i in range(n):
            mu = np.interp(energy, data['energy'][i], data['mu'][i])
            tau0 = np.interp(time, data['time'][i], data['tau0'][i])
            tau1 = np.interp(time, data['time'][i], data['tau1'][i])
            tau = tau0 + tau1 * data['fxy'][i,::sample,::sample]
            taus[i] = mu * tau
    return taus

def energy(fname, ccd_id, x, y, year):
    energy = None
    with astropy.io.fits.open(fname) as hdulist:
        hdu = hdulist[ccd_id+1]
        data = hdu.data
        n = data['component'].shape[0]
        taus = np.zeros((n, data['energy'].shape[1]))

        time = year2time(year)
        for i in range(n):
            if energy is None:
                energy = data['energy'][i]
            tau0 = np.interp(time, data['time'][i], data['tau0'][i])
            tau1 = np.interp(time, data['time'][i], data['tau1'][i])
            tau = tau0 + tau1 * data['fxy'][i][y][x]
            taus[i] = data['mu'][i] * tau

    return energy, taus

def year(fname, ccd_id, x, y, energy):
    times = None
    ntimes = 101
    with astropy.io.fits.open(fname) as hdulist:
        hdu = hdulist[ccd_id+1]
        data = hdu.data
        n = data['component'].shape[0]
        taus = np.zeros((n, ntimes))

        for i in range(n):
            if times is None:
                tmin = data['time'][i].min()
                tmax = data['time'][i].max()
                times = np.arange(ntimes)*(tmax-tmin)/(ntimes-1)+tmin
            mu = np.interp(energy, data['energy'][i], data['mu'][i])
            tau0 = np.interp(times, data['time'][i], data['tau0'][i])
            tau1 = np.interp(times, data['time'][i], data['tau1'][i])
            tau = tau0 + tau1 * data['fxy'][i][y][x]
            taus[i] = mu * tau

    return 1998+times/86400/365.2425, taus

def plot_image(args):
    sample=8
    taus = image(args.fname, args.ccd_id, args.energy, args.year, sample=sample)
    tau = taus.sum(axis=0)
    min, max, mean = tau.min(), tau.max(), tau.mean()
    print(f'{min=:.4g}, {max=:.4g}, {mean=:.4g}');
    fig = plt.figure()
    plt.imshow(tau, origin='lower', extent=(0.5, 1024.5, 0.5, 1024.5))
    plt.xlabel('CHIPX')
    plt.ylabel('CHIPY')
    plt.colorbar()
    plt.tight_layout()

def plot_energy(args):
    e, taus = energy(args.fname, args.ccd_id, args.x, args.y, args.year)
    tau = taus.sum(axis=0)
    min, max, mean = tau.min(), tau.max(), tau.mean()
    print(f'{min=:.4g}, {max=:.4g}, {mean=:.4g}');
    fig = plt.figure()
    plt.plot(e, tau)
    plt.xlabel('Energy')
    plt.ylabel(r'$\tau$')
    plt.xscale('log')
    plt.tight_layout()

def plot_year(args):
    y, taus = year(args.fname, args.ccd_id, args.x, args.y, args.energy)
    tau = taus.sum(axis=0)
    min, max, mean = tau.min(), tau.max(), tau.mean()
    print(f'{min=:.4g}, {max=:.4g}, {mean=:.4g}');
    fig = plt.figure()
    ii = y < 2026
    plt.plot(y[ii], tau[ii])
    plt.xlabel('Year')
    plt.ylabel(r'$\tau$')
    plt.xscale('log')
    plt.tight_layout()

def main():
    parser = argparse.ArgumentParser(
        description='Show a FITS HDU keyword.'
    )
    parser.add_argument('fname', help='FITS file name')
    parser.add_argument('ccd_id', type=int)
    parser.add_argument('x', type=int)
    parser.add_argument('y', type=int)
    parser.add_argument('energy', type=float)
    parser.add_argument('year', type=float, help='Year (e.g, 2026.5)')
    parser.add_argument('--image', action=argparse.BooleanOptionalAction, default=True, help='Display an image of optical depth over the full detector.')
    parser.add_argument('--energy', action=argparse.BooleanOptionalAction, help='Display a plot of optical depth vs energy.')
    parser.add_argument('--year', action=argparse.BooleanOptionalAction, help='Display an plot of optical depth vs time.')
    args = parser.parse_args()

    if args.image:
        plot_image(args)
    if args.energy:
        plot_energy(args)
    if args.year:
        plot_year(args)
    plt.show()


if __name__ == '__main__':
    main()
