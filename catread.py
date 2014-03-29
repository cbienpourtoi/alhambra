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
plt.show()
savemyplot("Distribution_parameter_stell")



# Molino's 'Stellar_Flag' parameter
fig = plt.figure()
plt.xlabel("Molino's stellaricity parameter 'Stellar_Flag' (1=star ; 0=galaxy)")
plt.title("Distribution of Molino's stellaricity parameter 'Stellar_Flag'")
plt.ylabel("Number of occurences")
plt.hist(t['Stellar_Flag'], bins = 10)
plt.show()
savemyplot("Distribution_parameter_stellar_flag")
plt.close()

# SExtractor 'stell' parameter VS Molino's 'Stellar_Flag' parameter
fig = plt.figure()
plt.title("SExtractor 'stell' parameter VS Molino's 'Stellar_Flag' parameter")
plt.xlabel("SExtractor stellaricity parameter 'stell'")
plt.ylabel("Molino's stellaricity parameter 'Stellar_Flag'")
plt.hist2d(t['stell'], t['Stellar_Flag'], bins=40, norm=LogNorm())
plt.colorbar()
#plt.show()
savemyplot("Molino_VS_SExtractor_Stellaricity")
plt.close()
## Conclusions of this graph:
## - At first order, if we forget about the line at y=0.5, the 2 parameters are coherent.
## - The line at y=0.5 is expected: SExtractor is bad at distinguishing for these faint (F814W>22.5) objects. Molino classifies them at .5, but they are probably galaxies (99% chances).   
## - The line at y=0.5 is stronger for low values of x, which is expected: it means SExtractor has a tendancy to classify these points as galaxies more than as stars. 



fig = plt.figure()
plt.title("SExtractor 'stell' and Molino's 'Stellar_Flag' VS F814W")
plt.xlabel("F814W magnitude")
plt.ylabel("stellaricity parameters")
plt.plot( t['F814W'], t['stell'],'.', label='SEx stell')
plt.plot( t['F814W'], t['Stellar_Flag'], '.', label='Molino Stellar_Flag')
plt.legend(loc=2)
#plt.show()
savemyplot("Molino_and_SExtractor_Stellaricity_VS_F814W")
plt.close()
## Conclusions of this graph:
## Molino considers everything with F814W>22.5 as 0.5, but strongly suggests it is a galaxy (99% chance).


#x = np.genfromtxt(catalog_filename, dtype=['i8','i'])#, comments='#')#, names=True)

#print x[1]


#with open(catalog_filename) as file:
#    for line in file:
#        # whatever you wanted to do with each line
