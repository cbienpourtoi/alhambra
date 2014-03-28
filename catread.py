''' Tool forreading catalogs from the Alhambra Survey '''
''' These catalogues are hosted on: https://cloud.iaa.csic.es/alhambra/catalogues/ '''
''' To use the catalogues, see Molino 2013 (soon 2014): http://arxiv.org/abs/1306.4968 (pdf: http://arxiv.org/pdf/1306.4968v2) '''


import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import scipy
from matplotlib.colors import LogNorm

plot_extension = ".png"
plot_directory = "../plots/"
def savemyplot(name):
	fig.savefig(plot_directory+name+plot_extension)
	return


catalogs_directory = '../data/catalogs/'


# Due to a bug in astropy.table (see :
# http://stackoverflow.com/questions/22617428/overflowerror-python-int-too-large-to-convert-to-c-long-with-astropy-table
# and
# https://github.com/astropy/astropy/pull/2234 )
# I truncate the first line of my catalogs until the bug is corrected.
# Truncated catalog name has a "MYVERSION." prefix.
catalog_filename = catalogs_directory+'MYVERSION.alhambra.F02P01C01.ColorProBPZ.cat'
#catalog_filename = catalogs_directory+'alhambra.F02P01C01.ColorProBPZ.cat'


t = Table.read(catalog_filename, format='ascii')

print "There are "+ str(len(t))+" elements in this catalogue"


# Selects only the objects that are not duplicated, or, if they are, selects only the best SNR object.
t = t[:][np.where(t['DupliDet_Flag']==0)]
print "There are "+ str(len(t))+" elements when low SNR duplicates are deleted"

## Selects only the objects that are not likely to be stars.
#t_noduplicates = t[:][np.where(t['DupliDet_Flag']==0)]
#print "There are "+ str(len(t_noduplicates))+" elements when low SNR duplicates are deleted"


####################################
# A short study about stellaricity #
####################################

# Catalog definitions:
# 'stell' corresponds to SExtractor stellaricity
# 'Stellar_Flag' is a much better indicator that takes in account many parameters: F814W, FWHM, J-Ks, F489W-F814W. 
#    See Molino's paper, section "Star/galaxy separation".


# SExtractor 'stell' parameter
fig = plt.figure()
plt.xlabel("SExtractor stellaricity parameter 'stell' (1=star ; 0=galaxy)")
plt.title("Distribution of SExtractor stellaricity parameter 'stell'")
plt.ylabel("Number of occurences")
plt.hist(t['stell'], bins = 10)
#plt.show()
savemyplot("Distribution_parameter_stell")
plt.close()


# Molino's 'Stellar_Flag' parameter
fig = plt.figure()
plt.xlabel("Molino's stellaricity parameter 'Stellar_Flag' (1=star ; 0=galaxy)")
plt.title("Distribution of Molino's stellaricity parameter 'Stellar_Flag'")
plt.ylabel("Number of occurences")
plt.hist(t['Stellar_Flag'], bins = 10)
plt.show()
savemyplot("Distribution_parameter_stellar_flag")
plt.close()

'''

fig = plt.figure()
plt.hist2d(t['stell'], t['Stellar_Flag'], bins=40, norm=LogNorm())
plt.colorbar()
#plt.show()
plt.close()

fig = plt.figure()
plt.plot( t['F814W'], t['stell'],'.')
plt.plot( t['F814W'], t['Stellar_Flag'], '.')
#plt.show()
plt.close()

'''

#x = np.genfromtxt(catalog_filename, dtype=['i8','i'])#, comments='#')#, names=True)

#print x[1]


#with open(catalog_filename) as file:
#    for line in file:
#        # whatever you wanted to do with each line
