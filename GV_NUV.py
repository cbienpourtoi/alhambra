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

import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column

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

    # Test mode? If true, only a small number of values will be
    # red from the master catalog
    test_mode = True
    if test_mode:
        number_lines_to_read = 1000
    else:
        number_lines_to_read = -1

    try:
        catalog = np.load(mastername)
    except:
        print "Cannot read the catalog", sys.exc_info()[0]
        raise

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


    slice = slice_z(catalog, 0.8, 0.05)
    slice = clean_filter(slice, filterB)
    slice = clean_filter(slice, filterR)

    BmR = slice[filterB] - slice[filterR]

    slice = Table(slice)
    slice.add_column(Column(data=BmR, name='BmR'))

    # selection in stellar mass:
    mass = 10.
    dmass = 0.2

    slicemass = slice_mass(slice, mass, dmass)

    #Red, Green, Blue limits:
    BGlim = 1.6
    GRlim = 1.9
    red = select_color(slicemass, GRlim)
    green = select_color(slicemass, BGlim, GRlim)
    blue = select_color(slicemass, -999, BGlim)

    print "Number of galaxies in the redshift slice: "+str(len(slice))
    print "Number of galaxies in the mass bin: "+str(len(slicemass))
    print "Number of blue galaxies: "+str(len(blue))
    print "Number of green galaxies: "+str(len(green))
    print "Number of red galaxies: "+str(len(red))

    print len(slice)

    plt.hist2d(slice["Stell_Mass_1"], slice['BmR'], bins=100)
    plt.axvline(mass-dmass, color='w')
    plt.axvline(mass+dmass, color='w')
    plt.axhline(1.9, color='w')
    plt.axhline(1.6, color='w')
    plt.show()



def slice_z(catalog, z, dz):
    """ select a slice in catalog from z-dz to z+dz """
    slice = catalog[np.where(np.abs(catalog["z_ml"]-z) < dz)]
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




if __name__ == '__main__':
    main()

