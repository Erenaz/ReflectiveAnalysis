import os
import pickle
import argparse
import numpy as np
import matplotlib.colors
import matplotlib.pyplot as plt
import astrotools.auger as auger
import NuRadioReco.utilities.cr_flux as cr_flux
from NuRadioReco.utilities import units

import itertools
marker = itertools.cycle((',', '+', '.', 'o', '*', 'x', '^', 'v', 'd'))
color = itertools.cycle(('black', 'blue', 'green', 'orange', 'brown'))
plt.style.use('plotsStyle.mplstyle')

def dBfactor(dB):
    return 10**(-dB/20)

def linefit(x, min=17):
    vals = .8/3 * (x-min) + 0.2
    vals[x<min] = 0
    vals[vals > 1] = 1
    return vals

def piecefit(x):
    vals = linefit(x)
    vals[x < 17.2] = 0
    return vals

#def logfit(x, min=17, a=0.14, b=0.16, c=1):    #Hand Fit
def logfit(x, min=17, a=0.49, b=0.22, c=0.78):    #Scipy Fit
    vals = np.log(a*(x-min-.3))*b + c
    vals[vals < 0] = 0
    vals[np.isnan(vals)] = 0
    return vals


def singleBinERate(eLow, eHigh, zLow, zHigh, Aeff):
#    print(f'check elow {eLow} ehigh {eHigh} zlow {zLow} zhigh {zHigh} aeff {Aeff}')
#    high = auger.event_rate(eLow, eHigh, zHigh, Aeff)
#    low = auger.event_rate(eLow, eHigh, zLow, Aeff)
#Need geometric exposure?
    high = auger.event_rate(eLow, eHigh, zHigh, Aeff * 0.5*(1 + np.cos(np.deg2rad(zHigh))))
    low = auger.event_rate(eLow, eHigh, zLow, Aeff * 0.5*(1 + np.cos(np.deg2rad(zLow))))

#What about using NuRadioReco cr_flux?
#    high = cr_flux.get_cr_event_rate(log10energy=eLow, zenith=zHigh*units.deg, a_eff=Aeff*units.km) * 10**eLow*units.eV
#    low = cr_flux.get_cr_event_rate(log10energy=eHigh, zenith=zLow*units.deg, a_eff=Aeff*units.km) * 10**eHigh*units.eV
#    print(f'for eLow-eHigh {eLow}-{eHigh} and zLow-zHigh {zLow}-{zHigh}')
#    print(f'high of {high/units.year} and low {low/units.year}')

    return high - low

def binnedAeff(energies, zeniths, nTrig, nThrown, trigMask, xLocs, yLocs, maxR, nRadBins=10, binWidth=.1, logScale=False):
    energy = energies
    if not logScale:
        energy = np.log10(energies)

    rads = np.sqrt(xLocs**2 + yLocs**2)

    engBins = np.arange(np.min(energy), np.max(energy), binWidth)
    cosZenBins = np.arange(1, 0.5, -0.05)
    zenBins = np.arccos(cosZenBins)
    zenBins[np.isnan(zenBins)] = 0

    radBins = np.linspace(0, maxR, nRadBins)
    nTrigBins = np.zeros((len(engBins)-1, len(zenBins)-1))
    nThrowBins = np.zeros_like(nTrigBins)

    rDigitize = np.digitize(rads, radBins)-1
    zDigitize = np.digitize(zeniths, zenBins)-1


    for iE in range(len(engBins)-1):
        lowMask = engBins[iE] <= energy
        highMask = energy < engBins[iE+1]
        engMask = lowMask & highMask
        for iZ in range(len(zenBins)-1):
            lowMask = zenBins[iZ] <= zeniths
            highMask = zeniths < zenBins[iZ+1]
            zenMask = lowMask & highMask

#            fMask = trigMask & engMask & zenMask
            fMask = engMask & zenMask
            tMask = trigMask[fMask]
            print(f'shape tmask {np.shape(tMask)}')
            print(f'shape before trigged {np.shape(nTrig)}')
            trigged = nTrig[fMask]
#            trigged = trigged[tMask]
            print(f'shape of trigged {np.shape(trigged)}')
            thrown = nThrown[fMask]
#            thrown = thrown[tMask]
            print(f'shape of rDigitize {np.shape(rDigitize)}')
            rdigits = rDigitize[fMask]
            print(f'shape of rdigits {np.shape(rdigits)}')
#            rdigits = rdigits[tMask]
            zdigits = zDigitize[fMask]
#            zdigits = zdigits[tMask]

            print(f'doing thrown {thrown} and trigged {trigged} and rdig {rdigits} and zdig {zdigits}')
            for iT in range(len(thrown)):
                rs = rdigits.flatten()
                for r in rs:
                    nThrowBins[r][zdigits[0]] += 1

            rdigits = rdigits[tMask]
#            if trigged.size == 0:
            for r in rdigits:
#                for iT in range(len(trigged)):
                nTrigBins[r][zdigits[0]] += 1

    print(f'final ntrig {nTrigBins}')
    print(f'final throw {nThrowBins}')

    nThrowBins[nThrowBins == 0] = 1
    nTrigBins = nTrigBins / nThrowBins
    aeffBins = np.zeros_like(nTrigBins)

    aeffBins = nTrigBins
#    radBins *= 10**-3   #convert from m to km
#    for iE in range(len(radBins)-1):
#        areaHigh = np.pi * radBins[iR+1]**2
#        areaLow = np.pi * radBins[iR]**2
#        aeffBins[iR] = nTrigBins[iR] * (areaHigh - areaLow)

    return aeffBins, nTrigBins, engBins, zenBins
"""
dd_rate = auger.event_rate(np.log10(np.min(energy[dd_mask])), np.log10(np.max(energy[dd_mask])), np.rad2deg(np.max(zenith[dd_mask])), Aeff_dd)
def returnAugerRate(aeffBins, radBins, zenBins):
    eRate = 0
    rCenter = (radBins[1:] + radBins[:-1])/2
    zCenter = (zenBins[1:] + zenBins[:-1])/2
    zCenter = np.rad2deg(zCenter)
    for iR in rCenter:
        for iZ in zCenter:
            
"""
def getERate(energy, zenith, xs, ys, trigMask, ntrig, nthrow, length, title='', plot=True):
    AREA = length**2
    maxR = np.sqrt( 2 * (length/2)**2 ) * 1000  #do in m
    radBins = np.linspace(0, maxR*1.1, 10)
#    radBins = np.array([0, maxR*1.1])
#    engBins = np.arange(np.log10(np.min(energy)), 20.1, .2)
#    cosZenBins = np.arange(1, -0.001, -0.05)
#    engBins = np.arange(15.9, 20.5, 0.2)
#    engBins = np.arange(15.9, 20.5, 0.5)
    engBins = np.arange(16, 20.1, 0.2)
    cosZenBins = np.arange(1, -0.001, -0.2)
#    cosZenBins = np.array([1, 0])
    zenBins = np.arccos(cosZenBins)
    zenBins[np.isnan(zenBins)] = 0

    rads = np.sqrt(xs**2 + ys**2)
    e_dig = np.digitize(np.log10(energy), engBins) - 1
    z_dig = np.digitize(zenith, zenBins) - 1
    r_dig = np.digitize(rads, radBins)-1

    trig = np.zeros( (len(engBins)-1, len(zenBins)-1, len(radBins)-1) )
    throw = np.zeros( (len(engBins)-1, len(zenBins)-1, len(radBins)-1) )

    trigBins = []
    throwBins = []
    engs = []
    for iE, eng in enumerate(energy):
        nE = e_dig[iE]
        zE = z_dig[iE]
        rE = r_dig[iE]
        for iR, r in enumerate(rE):
            if r == len(radBins)-1:
#                print(f'r of {r} in bins {radBins}, maxR {maxR}')
                continue
            if zE == len(zenBins)-1:
                print(f'ze of {zE} in bins {zenBins} from zen {zenith[iE]}')
#            throw[nE,zE,r] += 1
#            trig[nE,zE,r] += int(trigMask[iE,iR])
#            print(f'r dig {r} corresponding to {rads[iE][iR]}')
            throw[nE][zE][r] += 1
            trig[nE][zE][r] += int(trigMask[iE][iR])

#        if np.sum(trig[nE]) > 0:
#            print(f'for energy {eng} got a trig number of {np.sum(trig[nE])}')
#            trigBins.append(np.sum(trig[nE]))
#            engs.append(eng)
    for iE in range(len(engBins)-1):
        for iZ in range(len(zenBins)-1):
            trigBins.append(np.sum(trig[iE][iZ]))
            throwBins.append(np.sum(throw[iE][iZ]))
            engs.append(engBins[iE])

    plt.scatter(engs, trigBins, label='Trig')
    plt.scatter(engs, throwBins, label='Throw')
    plt.legend()
    plt.xlabel('Energy')
    plt.ylabel('Number Triggers')
    plt.yscale('log')
    plt.title(title)
    plt.show()



    throw[throw == 0] = 1
    Aeff = trig/throw

    for iR in range(len(radBins)-1):
        aHigh = np.pi * radBins[iR+1]*10**-3
        aLow = np.pi * radBins[iR]*10**-3
        a = aHigh - aLow
#        Aeff[:][:][iR] *= a
        Aeff[:,:,iR] *= a

    eRate = np.zeros_like(Aeff)
    print(f'len eng {len(engBins)} zen {len(zenBins)}')
    print(f'shape {np.shape(np.sum(eRate, axis=0))}')
    print(f'shape {np.shape(np.sum(eRate, axis=1))}')
    print(f'shape {np.shape(np.sum(eRate, axis=2))}')


        #Test starts here
    minE = 0
    for iE in range(len(engBins)-1):
        if np.any(Aeff[iE]):
            minE = engBins[iE]
            break
    eCenter = (engBins[1:] + engBins[:-1])/2
#    chiCut = linefit(eCenter, min=minE)
#    chiCut = piecefit(eCenter)
#    minE = 16.7
    chiCut = logfit(eCenter, min=minE)
    print(f'minE is {minE}')


    for iE in range(len(engBins)-1):
#        if engBins[iE+1] < 17.6:
#            continue
        for iZ in range(len(zenBins)-1):
            for iR in range(len(radBins)-1):
#                eRate[iE][iZ][iR] = singleBinERate(engBins[iE], engBins[iE+1], np.rad2deg(zenBins[iZ]), np.rad2deg(zenBins[iZ+1]), 1) * (10**((engBins[iE+1] + engBins[iE])/2))**2
                eRate[iE][iZ][iR] = singleBinERate(engBins[iE], engBins[iE+1], np.rad2deg(zenBins[iZ]), np.rad2deg(zenBins[iZ+1]), Aeff[iE][iZ][iR])
#                eRate[iE][iZ][iR] = singleBinERate(engBins[iE], engBins[iE+1], np.rad2deg(zenBins[iZ]), np.rad2deg(zenBins[iZ+1]), Aeff[iE][iZ][iR]) * chiCut[iE]


    asp = 4
    if plot:
#    while not asp == -1:
#        print(f'give aspect')
#        asp = float(input())


        rate = np.sum(eRate, axis=2)
        zenBins = np.rad2deg(zenBins)
        plt.imshow(rate.T, extent=[min(engBins), max(engBins), max(cosZenBins), min(cosZenBins)], aspect=asp, interpolation='none', norm=matplotlib.colors.LogNorm())
#        plt.imshow(rate.T, extent=[min(engBins), max(engBins), min(zenBins), max(zenBins)], aspect=1/25, interpolation='none', norm=matplotlib.colors.LogNorm())
        """
        eCenter = (engBins[1:] + engBins[:-1])/2
        zCenter = (zenBins[1:] + zenBins[:-1])/2
        zCenter = np.rad2deg(zCenter)
        for iE in range(len(eCenter)):
            for iZ in range(len(zCenter)):
#                rate = np.sum(eRate[iE][iZ])
                plt.scatter(eCenter[iE], zCenter[iZ], c=rate[iE][iZ])
        plt.clim(rate.min(), rate.max())
        """
        ax_labels = []
        for zen in zenBins:
            ax_labels.append('{:.0f}'.format(zen))
        ax = plt.gca()
        ax.set_yticks(cosZenBins)
        ax.set_yticklabels(ax_labels)
        plt.xlabel('Energy (log10eV)')
        plt.ylabel('CR Zenith (deg)')
        print(f'rate is {np.sum(rate)}')
        plt.colorbar(label=f'{np.sum(rate):.5f} Evts/Stn/Yr')
        plt.title(title)
#        plt.legend()
        plt.show()

        rate = np.sum(eRate, axis=2)
        rate = np.sum(rate, axis=1)
        plt.plot( (engBins[1:] + engBins[:-1])/2, rate, label=f'{np.sum(rate):.5f} Evts/Stn/Yr')
        plt.title(title)
        plt.xlabel('Energy (log10eV)')
        plt.ylabel('Evts/Stn/Yr')
#        plt.yscale('log')
        plt.show()

        arate = np.sum(Aeff, axis=2)
        arate = np.sum(arate, axis=1)
        plt.plot( (engBins[1:] + engBins[:-1])/2, arate, label=f'{np.sum(arate):.5f} Evts/Stn/Yr')
#        plt.title(title + ' for 35 stations')
        plt.title(title)
        plt.xlabel('Energy (log10eV)')
        plt.ylabel('Aeff (km^2)')
#        plt.yscale('log')
        plt.show()

    return eRate




parser = argparse.ArgumentParser(description='Run file on footprint refraction file')
parser.add_argument('file', type=str, default='None', help='file path to analyze')
parser.add_argument('area', type=float, default=5, help='Length of area thrown over in km, default 5km square')
parser.add_argument('--comment', type=str, default='', help='Title comment to add')

args = parser.parse_args()
file = args.file
area = args.area
comment = args.comment

#with open(f'data/FootprintData/NonRefractedFootprints_SPCR_Layer300.0m_40.0dB_LPDA35.0muV_dipole20muV_Area5.00km_100cores.pkl', 'rb') as fin:
#with open(f'data/MB_FootprintData/NonRefractedFootprints_MB_Layer576m_0dB_LPDA35.0muV_dipole20muV_Area5.00km_100cores.pkl', 'rb') as fin:
with open(file, 'rb') as fin:
    output = pickle.load(fin)


n = []
n_dip_dir = []
n_dip_refl = []
n_lpda_refl = []
n_lpda_dir = []
energy = []
zenith = []
azimuth = []
x = []
y = []
dip_dir_mask = []
dip_refl_mask = []
lpda_dir_mask = []
lpda_refl_mask = []
ant_zen = []
dip_dir_SNR = []
dip_refl_SNR = []
lpda_dir_SNR = []
lpda_refl_SNR = []

for runid in output:
    print(runid)

#    if np.log10(output[runid]['energy']) < 17.7:
#        continue

    n.append(output[runid]['n'])
    n_dip_dir.append(output[runid]['n_dip_dir'])
    n_dip_refl.append(output[runid]['n_dip_refl'])
    n_lpda_refl.append(output[runid]['n_lpda_refl'])
    n_lpda_dir.append(output[runid]['n_lpda_dir'])
    energy.append(output[runid]['energy'])
    zenith.append(output[runid]['zenith'])
    azimuth.append(output[runid]['azimuth'])
    x.append(output[runid]['x_dir_lpda'])
    y.append(output[runid]['y_dir_lpda'])
    dip_dir_mask.append(output[runid]['dip_dir_mask'])
    dip_refl_mask.append(output[runid]['dip_refl_mask'])
    lpda_dir_mask.append(output[runid]['lpda_dir_mask'])
    lpda_refl_mask.append(output[runid]['lpda_refl_mask'])
    ant_zen.append(output[runid]['ant_zen'])
    dip_dir_SNR.append(output[runid]['dip_dir_SNR'])
    dip_refl_SNR.append(output[runid]['dip_refl_SNR'])
    lpda_dir_SNR.append(output[runid]['lpda_dir_SNR'])
    lpda_refl_SNR.append(output[runid]['lpda_refl_SNR'])

n = np.array(n)
n_dip_dir = np.array(n_dip_dir)
n_dip_refl = np.array(n_dip_refl)
n_lpda_refl = np.array(n_lpda_refl)
n_lpda_dir = np.array(n_lpda_dir)
energy = np.array(energy)
zenith = np.array(zenith)
azimuth = np.array(azimuth)
#print(f'x is {x}')
#x = np.array(x, dtype=object)
#y = np.array(y, dtype=object)
#x = np.asarray(x)
#y = np.asarray(y)
x = np.array(list(itertools.zip_longest(*x, fillvalue=np.NaN))).T
y = np.array(list(itertools.zip_longest(*y, fillvalue=np.NaN))).T

dip_dir_mask = np.array(dip_dir_mask, dtype=object)
dip_refl_mask = np.array(dip_refl_mask, dtype=object)
lpda_refl_mask = np.array(lpda_refl_mask, dtype=object)
#print(f'lpda refl mask {lpda_refl_mask}')
lpda_dir_mask = np.array(lpda_dir_mask, dtype=object)
ant_zen = np.array(ant_zen)
#dip_dir_SNR = np.array(dip_dir_SNR)
#dip_refl_SNR = np.array(dip_refl_SNR)
#print(f'dir dir snr {lpda_dir_SNR}')
#lpda_dir_SNR = np.array(lpda_dir_SNR)
#print(f'dir dir snr {lpda_dir_SNR}')
#lpda_refl_SNR = np.array(lpda_refl_SNR)


"""
length = 1000  #m
maxR = np.sqrt( 2 * (length/2)**2 )
dlAeff, dlTrig, dlEngBins, dlZenBins = binnedAeff(energy, zenith, n_lpda_dir, n, lpda_dir_mask, x_dir_lpda, y_dir_lpda, maxR)
dlAeff = dlTrig

eCenter = (dlEngBins[1:] + dlEngBins[:-1])/2
zCenter = (dlZenBins[1:] + dlZenBins[:-1])/2
print(f'Aeff orig {dlAeff}')
print(f'Transposed {dlAeff}')


for iZ, zen in enumerate(zCenter):
    print(f'shape ecent {np.shape(eCenter)} and dl segment {np.shape(dlAeff[:][iZ])}')
    plt.scatter(eCenter, dlAeff[:][iZ], label=f'Zens {dlZenBins[iZ]}-{dlZenBins[iZ+1]}')
    plt.xlabel('Energy (log10 eV)')
    plt.legend()
    plt.title('Aeff of direct LPDA triggers')
    plt.show()
"""

#Code for plotting trigger eff vs SNR
#lpda_dir_SNR = lpda_dir_SNR / 0.001
maxes = max(lpda_dir_SNR)
SNR_bins = np.linspace(0, 0.02, num=24)
#SNR_bins = np.linspace(0, max(maxes)/40, num=100)
#snrDigs = np.digitize(lpda_dir_SNR, SNR_bins) - 1

throws = np.zeros(len(SNR_bins)-1)
trigs = np.zeros(len(SNR_bins)-1)

for iS, SNRs in enumerate(lpda_dir_SNR):
    if len(SNRs) == 0:
        continue
#    print(f'snrs {SNRs} and bins {SNR_bins}')
    snrDigs = np.digitize(SNRs, SNR_bins)-1
    print(f'len snrs {len(SNRs)} and mask {len(lpda_dir_mask[iS])}')
    for iD, dig in enumerate(snrDigs):
        if dig >= len(SNR_bins)-1:
            print(f'dig too big')
            continue
#        print(f'dig {dig} and iS {iS}')
#        trigs[dig] += np.sum(n_lpda_dir[iS])
#        throws[dig] += np.sum(n[iS])
        throws[dig] += 1
        trigs[dig] += int(lpda_dir_mask[iS][iD])


trigs[throws == 0] = 1
throws[throws == 0] = 0
plotSNR = (SNR_bins[1:] + SNR_bins[:-1])/(2*0.0015)
plt.plot(plotSNR, trigs/throws)
plt.xlabel('SNR')
plt.ylabel('Trigger efficiency')
plt.title(comment)
plt.show()
quit()


#GOOD CODE, for event rate
if True:
    RLeRate = getERate(energy, zenith, x, y, lpda_refl_mask, n_lpda_refl, n, area, title='Reflected LPDA '+comment)
    print(f'total eRate refl LPDA {np.sum(RLeRate)}')
    DDeRate = getERate(energy, zenith, x, y, dip_dir_mask, n_dip_dir, n, area, title='Direct Dipole '+comment)
    print(f'total eRate direct dipole {np.sum(DDeRate)}')
    DLeRate = getERate(energy, zenith, x, y, lpda_dir_mask, n_lpda_dir, n, area, title='Direct LPDA '+comment)
    print(f'total eRate direct LPDA {np.sum(DLeRate)}')
    RDeRate = getERate(energy, zenith, x, y, dip_refl_mask, n_dip_refl, n, area, title='Reflected Dipole '+comment)
    print(f'total eRate refl dipole {np.sum(RDeRate)}')

    engBins = np.arange(np.log10(np.min(energy)), 20.1, .1)
    DDengRate = []
    DLengRate = []
    RDengRate = []
    RLengRate = []
    for iE in range(len(engBins)-1):
        DDengRate.append(np.sum(DDeRate[iE]))
        DLengRate.append(np.sum(DLeRate[iE]))
        RDengRate.append(np.sum(RDeRate[iE]))
        RLengRate.append(np.sum(RLeRate[iE]))
    plt.plot( (engBins[1:]+engBins[:-1])/2, DLengRate, label=f'{np.sum(DLengRate):2f} evts/stn/yr')
    plt.legend()
    plt.yscale('log')
    plt.xlabel('Energy (log10eV)')
    plt.ylabel('Event Rate (evts/stn/yr)')
    plt.title('Event Rate per Energy for Upwardfacing LPDA')
    plt.show()

    plt.plot( (engBins[1:]+engBins[:-1])/2, DDengRate, label=f'Direct Dipole {np.sum(DDengRate):.4f} evts/stn/yr', color=next(color))
    plt.plot( (engBins[1:]+engBins[:-1])/2, RDengRate, label=f'Reflected Dipole {np.sum(RDengRate):.4f} evts/stn/yr', color=next(color))
    plt.plot( (engBins[1:]+engBins[:-1])/2, RLengRate, label=f'Reflected LPDA {np.sum(RLengRate):.4f} evts/stn/yr', color=next(color))
    plt.legend()
    plt.yscale('log')
    plt.xlabel('Energy (log10eV)')
    plt.ylabel('Event Rate (evts/stn/yr)')
    plt.title('Event Rate per Energy, SP')
    plt.show()



mask = np.any(lpda_refl_mask, axis=1)
print(f'mask shape {np.shape(mask)} and energy {np.shape(energy)}')
et = energy[mask]
print(f'et {et}')
print(f'true energies are dip dir {energy[np.any(dip_dir_mask, axis=1)]}, dip refl {energy[np.any(dip_refl_mask, axis=1)]}, lpda refl {energy[np.any(lpda_refl_mask, axis=1)]}')
#quit()

print(f' dip dir {np.sum(dip_dir_mask)} dip relf {np.sum(dip_refl_mask)} lpda refl {np.sum(lpda_refl_mask)} lpda dir {np.sum(lpda_dir_mask)}')

plt.scatter(x[dip_dir_mask], y[dip_dir_mask])
plt.title('Direct Dipole xy positions')
plt.show()

plt.scatter(x[dip_refl_mask], y[dip_refl_mask])
plt.title('Reflected Dipole xy Positions')
plt.show()

plt.scatter(x[lpda_refl_mask], y[lpda_refl_mask])
plt.title('Reflected LPDA xy Positions')
plt.show()

plt.scatter(x[lpda_dir_mask], y[lpda_dir_mask])
plt.title('LPDA Triggered Direct core locations')
plt.show()

#area = 3
area = area**2


dd_mask = np.any(dip_dir_mask, axis=1)
rd_mask = np.any(dip_refl_mask, axis=1)
rl_mask = np.any(lpda_refl_mask, axis=1)
dl_mask = np.any(lpda_dir_mask, axis=1)

print(f'n {np.sum(n)} dd {np.sum(n_dip_dir)} rd {np.sum(n_dip_refl)} rl {np.sum(n_lpda_refl)}')
Aeff_dd = area * np.sum(n_dip_dir) / np.sum(n)
Aeff_rd = area * np.sum(n_dip_refl) / np.sum(n)
Aeff_rl = area * np.sum(n_lpda_refl) / np.sum(n)
#Aeff_rl = area * np.sum(n_lpda_refl[rl_mask]) / np.sum(n[rl_mask])
Aeff_dl = area * np.sum(n_lpda_dir) / np.sum(n)

dd_rate = auger.event_rate(np.log10(np.min(energy[dd_mask])), np.log10(np.max(energy[dd_mask])), np.rad2deg(np.max(zenith[dd_mask])), Aeff_dd)
rd_rate = auger.event_rate(np.log10(np.min(energy[rd_mask])), np.log10(np.max(energy[rd_mask])), np.rad2deg(np.max(zenith[rd_mask])), Aeff_rd)
rl_rate = auger.event_rate(np.log10(np.min(energy[rl_mask])), np.log10(np.max(energy[rl_mask])), np.rad2deg(np.max(zenith[rl_mask])), Aeff_rl)
dl_rate = auger.event_rate(np.log10(np.min(energy[dl_mask])), np.log10(np.max(energy[dl_mask])), np.rad2deg(np.max(zenith[dl_mask])), Aeff_dl)

print(f'Event rate dir dip {dd_rate}')
print(f'Event rate refl dip {rd_rate}')
print(f'Event rate refl lpda {rl_rate}')
print(f'Event rate dir lpda {dl_rate}')

#plt.scatter(np.log10(energy[~dd_mask]), np.rad2deg(zenith[~dd_mask]), color='black', marker='x')
#plt.scatter(np.log10(energy[dd_mask]), np.rad2deg(zenith[dd_mask]), c=n_dip_dir[dd_mask] / n[dd_mask])
plt.scatter(np.log10(energy[~dd_mask]), 1-np.cos(zenith[~dd_mask]), color='black', marker='x')
plt.scatter(np.log10(energy[dd_mask]), 1-np.cos(zenith[dd_mask]), c=n_dip_dir[dd_mask] / n[dd_mask])
plt.ylabel('1-cos(zenith)')
plt.title('Footprint Triggers of Direct Dipole')
plt.show()

plt.scatter(np.log10(energy[~rd_mask]), np.rad2deg(zenith[~rd_mask]), color='black', marker='x')
plt.scatter(np.log10(energy[rd_mask]), np.rad2deg(zenith[rd_mask]), c=n_dip_refl[rd_mask] / n[rd_mask])
plt.title('Footprint Triggers of Reflected Dipole')
plt.show()


plt.scatter(np.log10(energy[~rl_mask]), np.rad2deg(zenith[~rl_mask]), color='black', marker='x')
plt.scatter(np.log10(energy[rl_mask]), np.rad2deg(zenith[rl_mask]), c=n_lpda_refl[rl_mask] / n[rl_mask])
plt.title('Footprint Triggers of Reflected LPDA')
plt.show()

plt.scatter(np.log10(energy[~dl_mask]), np.rad2deg(zenith[~dl_mask]), color='black', marker='x')
plt.scatter(np.log10(energy[dl_mask]), np.rad2deg(zenith[dl_mask]), c=n_lpda_dir[dl_mask] / n[dl_mask])
plt.title('Footprint Triggers of Direct LPDA')
plt.show()



plt.hist(np.rad2deg(ant_zen[dd_mask].flatten()))
plt.xlabel('Dipole Direct Recieve Zenith')
plt.show()

plt.hist(np.rad2deg(ant_zen[rd_mask].flatten()))
plt.xlabel('Dipole Reflected Recieve Zenith')
plt.show()

plt.hist(np.rad2deg(ant_zen[rl_mask].flatten()))
plt.xlabel('LPDA Reflected Recieve Zenith')
plt.show()



