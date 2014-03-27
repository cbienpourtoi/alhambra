''' Tool forreading catalogs from the Alhambra Survey '''
''' These catalogues are hosted on: https://cloud.iaa.csic.es/alhambra/catalogues/ '''
''' To use the catalogues, see Molino 2013 (soon 2014): http://arxiv.org/abs/1306.4968 (pdf: http://arxiv.org/pdf/1306.4968v2) '''


import numpy as np
#import scipy
import matplotlib.pyplot as plt
from astropy.table import Table


# Due to a bug in astropy.table (see :
# http://stackoverflow.com/questions/22617428/overflowerror-python-int-too-large-to-convert-to-c-long-with-astropy-table
# and
# https://github.com/astropy/astropy/pull/2234 )
# I truncate the first line of my catalogs until the bug is corrected.
# Truncated catalog name has a "MYVERSION." prefix.
catalog_filename = '/home/loic/data/alhambra/catalogs/MYVERSION.alhambra.F02P01C01.ColorProBPZ.cat'
#catalog_filename = '/home/loic/data/alhambra/catalogs/alhambra.F02P01C01.ColorProBPZ.cat'


t = Table.read(catalog_filename, format='ascii')

print "There are "+ str(len(t))+" elements in this catalogue"


# Selects only the objects that are not duplicated, or, if they are, selects only the best SNR object.
t_noduplicates = t[:][np.where(t['DupliDet_Flag']==0)]
print "There are "+ str(len(t_noduplicates))+" elements when low SNR duplicates are deleted"

# Selects only the objects that are not likely to be stars.
t_noduplicates = t[:][np.where(t['DupliDet_Flag']==0)]
print "There are "+ str(len(t_noduplicates))+" elements when low SNR duplicates are deleted"




plt.hist(t['stell'], bins = 10)
plt.show()



#x = np.genfromtxt(catalog_filename, dtype=['i8','i'])#, comments='#')#, names=True)

#print x[1]


#with open(catalog_filename) as file:
#    for line in file:
#        # whatever you wanted to do with each line
