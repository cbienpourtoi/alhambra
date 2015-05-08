#!/usr/bin/env python

"""GV_NUV.py: Makes a graph of the Green Valley
from NUV-r color using Alhambra data.
"""

__author__ = "Loic Le Tiran"
__copyright__ = "Copyright 2015"
__credits__ = "Loic Le Tiran"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Loic Le Tiran"
__email__ = "loic.letiran@gmail.com"
__status__ = "Development"

import sys, os
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
from astropy.cosmology import Planck13 as cosmo
from astropy import units as u
from astropy.cosmology import z_at_value
from astropy.io import ascii
import datetime

def main():
    mastername = '/home/loic/Projects/alhambra/data/catalogs/' \
               'master_catalogs/alhambra.Master.ColorProBPZ.cat.npy'
               # Numpy format of the master catalog
                # It contains already a pre-classification of galaxies
    catalog = read_catalog(mastername)

    plot_GV_zfixed(catalog)

def read_catalog(mastername):
    """Reads the Alhambra "master" catalog
    This catalog is created in the catread.py code
    and makes a large catalog from all the different
    subcatalogs downloaded from the Alhambra site,
    with a number of cuts (stellaricity etc...)
    """

    try:
        catalog = Table(np.load(mastername))
    except:
        print "Cannot read the catalog", sys.exc_info()[0]
        raise

    test_mode = False
    if test_mode:
        catalog = catalog[:1000]

    return catalog

def plot_GV_zfixed(catalog):
    """ plots a Green Valley at z=0.8
        Galex NUV: effective lambda = 2267A
        (http://galexgi.gsfc.nasa.gov/docs/galex/Documents/ERO_data_description_2.htm
        table 1.1)
        for z=0.8 : 2267*1.8 = 4081A, about F396W
        and SDSSr = 6231;
        6231*1.8=11216
        6231*[1.75, 1.85] = [10904, 11527]
        Alhambra filters: F954W, or J=1.24 um
        let's take J for the moment

    """

    filterB = "F644W"
    filterR = "J"

    #####################
    # Slice to analyse: #
    #####################

    slice = slice_z(catalog, 0.8, 0.05)
    slice = clean_filter(slice, filterB)
    slice = clean_filter(slice, filterR)

    BmR = slice[filterB] - slice[filterR]

    slice.add_column(Column(data=BmR, name='BmR'))

    # selection in stellar mass:
    mass = 10.
    dmass = 0.2

    slicemass = slice_mass(slice, mass, dmass)


    #####################
    #   Larger slice :  #
    #####################
    larger_slice = slice_z(catalog, 0.8, 0.1)
    massmin_large = 8.4
    larger_slice = slice_param(catalog, "Stell_Mass_1", massmin_large, 10000.)

    #############
    # Densities #
    #############
    envdir = "environment/"
    if not os.path.exists(envdir):
        os.mkdir(envdir)

    do_measurement = True
    if do_measurement:
        # Measures the densities in environment
        densities = get_densities(slicemass, larger_slice)
        # Saves the densities
        print densities
        np.save(envdir+"densities.npy", densities)
        ascii.write(densities, envdir+"densities.txt")
    else:
        densities = Table(np.load(envdir+"densities.npy"))



    #Red, Green, Blue limits:
    BGlim = 1.6
    GRlim = 1.9
    red = select_color(densities, GRlim)
    green = select_color(densities, BGlim, GRlim)
    blue = select_color(densities, -999, BGlim)

    print "Number of galaxies in the redshift slice: "+str(len(slice))
    print "Number of galaxies in the mass bin: "+str(len(slicemass))
    print "Number of blue galaxies: "+str(len(blue))
    print "Number of green galaxies: "+str(len(green))
    print "Number of red galaxies: "+str(len(red))

    print "Mean blue env objects: "+str(np.mean(blue["densities"]))
    print "Mean red env objects: "+str(np.mean(red["densities"]))
    print "Mean green env objects: "+str(np.mean(green["densities"]))


    plt.hist(blue["densities"], color='b', label="blue", normed=True, alpha=.3, range=(0,10), histtype='step')
    plt.hist(red["densities"], color='r', label="red", normed=True, alpha=.3, range=(0,10), histtype='step')
    plt.hist(green["densities"], color='g', label="green", normed=True, alpha=.3, range=(0,10), histtype='step')
    plt.legend()
    plt.show()

    sys.exit()

    plt.hist2d(slice["Stell_Mass_1"], slice['BmR'], bins=100)
    plt.axvline(mass-dmass, color='w')
    plt.axvline(mass+dmass, color='w')
    plt.axhline(GRlim, color='w')
    plt.axhline(BGlim, color='w')
    plt.axvline(massmin_large, color='r')
    plt.show()
    plt.savefig(envdir+"GV.png")
    plt.close()


def slice_z(catalog, z, dz):
    """ select a slice in catalog from z-dz to z+dz """
    slice = catalog[np.where(np.abs(catalog["zb_1"]-z) < dz)]
    return slice

def slice_mass(catalog, mass, dmass):
    """ select a slice in catalog from stelar mass-dmass to mass+dmass """
    slice = catalog[np.where(np.abs(catalog["Stell_Mass_1"]-mass) < dmass)]
    return slice

def clean_filter(slice, filter):
    """ deletes the values set at 99 for the indicated filter """
    slice = slice[np.where(np.abs(slice[filter]) != 99.0)]
    return slice

def select_color(slice, BmRmin = -999, BmRmax = 999):
    """ Selects only the red, blue, or green galaxies """
    slice = slice[np.where(slice["BmR"]>BmRmin)]
    slice = slice[np.where(slice["BmR"]<BmRmax)]
    return slice

def slice_param(catalog, param, pmin, pmax):
    """ select a slice in catalog from
    whatever parameter from min to max"""
    slice = catalog[np.where(catalog[param] < pmax)]
    slice = slice[np.where(catalog[param] > pmin)]
    return slice


def get_densities(slice, large_slice):
    """Measures the density for each galaxy in the slice
    using a larger slice, because of borders effects from z
    (RA Dec not taken in accounts) and also to include all masses
    and not only the current mass bin"""

    radius = 5.*u.Mpc

    densities = np.array([])
    meanmasses = np.array([])
    medianmasses = np.array([])

    print datetime.datetime.now()

    for galaxy in slice:
        z = galaxy["zb_1"]
        RA = galaxy["RA"]
        Dec = galaxy["Dec"]

        Mpc_per_deg = cosmo.kpc_comoving_per_arcmin(z).to(u.Mpc/u.deg)
        dRA = radius / Mpc_per_deg
        dDec = dRA

        slicetmp = large_slice

        slicetmp = slice_RA(slicetmp, RA, dRA)
        slicetmp = slice_Dec(slicetmp, Dec, dDec)
        slicetmp = slice_z_cosmo(slicetmp, z, radius)

        meanmass_comp = np.mean(slicetmp["Stell_Mass_1"][np.where(slicetmp["ID"] != galaxy["ID"])])
        medianmass_comp = np.median(slicetmp["Stell_Mass_1"][np.where(slicetmp["ID"] != galaxy["ID"])])

        densities = np.append(densities, len(slicetmp)-1)
        meanmasses = np.append(meanmasses, meanmass_comp)
        medianmasses = np.append(medianmasses, medianmass_comp)

        #print len(slicetmp)-1

    slice.add_column(Column(data=densities, name="densities"))
    slice.add_column(Column(data=meanmasses, name="meanmasses"))
    slice.add_column(Column(data=medianmasses, name="medianmasses"))

    print datetime.datetime.now()

    return slice

def slice_z_cosmo(slice, z, radius):
    """Slices in redshift considering a distance to cut"""
    zmax = z_at_value(cosmo.comoving_distance, cosmo.comoving_distance(z) + radius)
    zmin = z_at_value(cosmo.comoving_distance, cosmo.comoving_distance(z) - radius)
    slice = slice[np.where(slice["zb_1"] < zmax)]
    slice = slice[np.where(slice["zb_1"] > zmin)]
    return slice

def slice_RA(slice, RA, dRA):
    """ In a table, gets only RA inside RA-dRA, and RA+dRA """
    slice = slice[np.where(np.abs(slice["RA"]-RA)*u.deg < dRA)]
    return slice

def slice_Dec(slice, Dec, dDec):
    """ In a table, gets only Dec inside Dec-dDec, and Dec+dDec """
    slice = slice[np.where(np.abs(slice["Dec"]-Dec)*u.deg < dDec)]
    return slice


if __name__ == '__main__':
    main()

