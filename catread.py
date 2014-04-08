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

# If True, the code will study some parameters and make some plots. Useless usually.
study_stellaricity = False
study_Flags = False
study_Odds = False
study_z = False

catalogs_directory = '../data/catalogs/'

# Values of the cuts:
Stellar_Flag_Cut = 0.51
Odds_1_Cut = 0.2
irms_OPT_Flag_Cut = 1
irms_NIR_Flag_Cut = 1
photoflag_Cut = 1


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



####################################
# Selection : no duplicates        #
####################################

# Selects only the objects that are not duplicated, or, if they are, selects only the best SNR object.
t = t[:][np.where(t['DupliDet_Flag']==0)]
print "There are "+ str(len(t))+" unique elements when low SNR duplicates are deleted"


####################################
# A short study about stellaricity #
####################################

if study_stellaricity:

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
	plt.clf()


	# Molino's 'Stellar_Flag' parameter
	fig = plt.figure()
	plt.xlabel("Molino's stellaricity parameter 'Stellar_Flag' (1=star ; 0=galaxy)")
	plt.title("Distribution of Molino's stellaricity parameter 'Stellar_Flag'")
	plt.ylabel("Number of occurences")
	plt.hist(t['Stellar_Flag'], bins = 10)
	#plt.show()
	savemyplot("Distribution_parameter_stellar_flag")
	plt.clf()

	# SExtractor 'stell' parameter VS Molino's 'Stellar_Flag' parameter
	fig = plt.figure()
	plt.title("SExtractor 'stell' parameter VS Molino's 'Stellar_Flag' parameter")
	plt.xlabel("SExtractor stellaricity parameter 'stell'")
	plt.ylabel("Molino's stellaricity parameter 'Stellar_Flag'")
	plt.hist2d(t['stell'], t['Stellar_Flag'], bins=40, norm=LogNorm())
	plt.colorbar()
	#plt.show()
	savemyplot("Molino_VS_SExtractor_Stellaricity")
	plt.clf()
	## Conclusions of this graph:
	## - At first order, if we forget about the line at y=0.5, the 2 parameters are coherent.
	## - The line at y=0.5 is expected: SExtractor is bad at distinguishing for these faint (F814W>22.5) objects. Molino classifies them at .5, but they are probably galaxies (99% chances).   
	## - The line at y=0.5 is stronger for low values of x, which is expected: it means SExtractor has a tendancy to classify these points as galaxies more than as stars. 


	# SExtractor 'stell' parameter and Molino's 'Stellar_Flag' parameter VS F814W magnitude
	fig = plt.figure()
	plt.title("SExtractor 'stell' and Molino's 'Stellar_Flag' VS F814W")
	plt.xlabel("F814W magnitude")
	plt.ylabel("stellaricity parameters")
	plt.plot( t['F814W'], t['stell'],'.', label='SEx stell')
	plt.plot( t['F814W'], t['Stellar_Flag'], '.', label='Molino Stellar_Flag')
	plt.legend(loc=2)
	#plt.show()
	savemyplot("Molino_and_SExtractor_Stellaricity_VS_F814W")
	plt.clf()
	## Conclusions of this graph:
	## Molino considers everything with F814W>22.5 as 0.5, but strongly suggests it is a galaxy (99% chance).



####################################
# Selection : only galaxies        #
####################################

# Selects only the galaxies, i.e. objects who have Stellar_Flag<0.7 (according to Molino's mail 31/3/14)
t = t[:][np.where(t['Stellar_Flag']<Stellar_Flag_Cut)]
print "There are "+ str(len(t))+" unique elements classified as galaxies (Stellar_Flag<0.7)"


####################################
# Study : photoflag,               #
# irms_OPT_Flag, irms_NIR_Flag 	   #
####################################

if study_Flags:

	# photoflag
	fig = plt.figure()
	plt.title("photoflag distribution (galaxies only)")
	plt.xlabel("photoflag")
	plt.ylabel("Number")
	plt.hist(t['photoflag'], bins=20)
	#plt.show()
	savemyplot("photoflag_distribution")
	plt.clf()

	# irms_OPT_Flag
	fig = plt.figure()
	plt.title("irms_OPT_Flag distribution (galaxies only)")
	plt.xlabel("irms_OPT_Flag")
	plt.ylabel("Number")
	plt.hist(t['irms_OPT_Flag'], bins=100)
	#plt.show()
	savemyplot("irms_OPT_Flag_distribution")
	plt.clf()

	# irms_NIR_Flag
	fig = plt.figure()
	plt.title("irms_NIR_Flag distribution (galaxies only)")
	plt.xlabel("irms_NIR_Flag")
	plt.ylabel("Number")
	plt.hist(t['irms_NIR_Flag'], bins=10)
	#plt.show()
	savemyplot("irms_NIR_Flag_distribution")
	plt.clf()



####################################
# Study : Odds                     #
####################################


if study_Odds:

	# Odds vs mag
	fig = plt.figure()
	plt.title("Odds VS F814W")
	plt.xlabel("F814W magnitude")
	plt.ylabel("Odds")
	plt.hist2d(t['F814W'], t['Odds_1'], bins=40, norm=LogNorm())
	plt.plot([0,100],[Odds_1_Cut,Odds_1_Cut], '-', label="Odds Cut")
	plt.colorbar()
	plt.legend()
	#plt.show()
	#savemyplot("Odds_vs_F814W")
	plt.clf()



####################################
# Selection : Odds, PhotoFlag,     #
# irms_OPT_Flag, irms_NIR_Flag     #
####################################

# "Selecting Odds>0.2 will remove several (unreliable) faint galaxies." (according to Molino's mail 31/3/14)
t = t[:][np.where(t['Odds_1']>Odds_1_Cut)]
print "There are "+ str(len(t))+" unique elements classified as galaxies with Odds better than "+str(Odds_1_Cut)


# "irms_OPT_Flag" & "irms_NIR_Flag"

# "irms_OPT_Flag" & "irms_NIR_Flag" indicate the number of
# passbands (out of 23) in the Optican and in the NIR (respect.) a galaxy
# was poorly observed due to an insufficient exposure time. In other words,
# detections with values larger that irms_OPT_Flag>3 and irms_NIR_Flag>1
# could have a suspicious photo-z since the photometry was not as clean as
# the other galaxies." (according to Molino's mail 31/3/14)

t = t[:][np.where(t['irms_OPT_Flag']<irms_OPT_Flag_Cut)]
print "Remains   "+ str(len(t))+" elements with less than "+str(irms_OPT_Flag_Cut)+ " optical bands in which there is no detection"

t = t[:][np.where(t['irms_NIR_Flag']<irms_NIR_Flag_Cut)]
print "Remains   "+ str(len(t))+" elements with less than "+str(irms_NIR_Flag_Cut)+ " NIR bands in which there is no detection"


# photoflag

#"If you are interested in few galaxies you can look at the "photo flag"
#which is a photometric Flag based on SExtractor. It warns you whether the
#galaxy was photometrically compromised when estimating its fluxes (due to
#a neighbor galaxy, has a saturated pixel,...) or the photometry was clean."
#(according to Molino's mail 31/3/14)

#SExtractor manual. Flag parameter at sections 9 and 9.1
#https://www.astromatic.net/pubsvn/software/sextractor/trunk/doc/sextractor.pdf
#Flags:
#1      The object has neighbours, bright and close enough to significantly bias the MAGAUTO photometry, or bad pixels (more than 10% of the integrated area affected)
#2      The object was originally blended with another one,
#4      At least one pixel of the object is saturated (or very close to),
#8      The object is truncated (too close to an image boundary),
#16     Object's aperture data are incomplete or corrupted,
#32     Object's isophotal data are incomplete or corrupted
#64     A memory overflow occurred during deblending,
#128    A memory overflow occurred during extraction.

t = t[:][np.where(t['photoflag']<photoflag_Cut)]
print "Remains   "+ str(len(t))+" elements with photoflag<"+str(photoflag_Cut)+ "."





####################################
# Study : Redshifts                #
####################################


# Bayesian approach :
# zb is the Bayesian value of theredshift (and tb_1 for spectra-type)

# No Prior (to avoid) :
# z_ml (and t_ml) represent again the photometric redshift 
# and its spectral-type but estimated using a maximum 
# likelihood without priors (instead of a Bayesian approach).



if study_z:

	# Redshifts distribution
	fig = plt.figure()
	plt.title("Redshifts distribution")
	plt.xlabel("Redshifts (zb - Bayesian)")
	plt.ylabel("#")
	plt.hist(t['zb_1'], bins=50)
	#plt.hist(t['z_ml'], bins=50)
	#plt.hist2d(t['z_ml'], t['zb_1'])
	plt.show()
	savemyplot("z_dist")
	plt.clf()

	fig = plt.figure()
	plt.title("Bayesian vs no-prior approach Redshifts")
	plt.ylabel("Redshifts (zb - Bayesian)")
	plt.xlabel("Redshifts (z_ml - no prior)")
	plt.hist2d(t['z_ml'], t['zb_1'], bins = 30, norm=LogNorm())
	plt.show()
	savemyplot("z_Bayes_vs_noprior")
	plt.clf()


t_z0 = t[:][np.where(t['zb_1']<0.2)]
print "Remains   "+ str(len(t_z0))+" elements with z<0.2"



# Odds vs mag
fig = plt.figure()
#plt.title("Odds VS F814W")
#plt.xlabel("F814W magnitude")
#plt.ylabel("Odds")
#plt.hist(t_z0['F675W'], bins = 50)
plt.hist(t_z0['F458W'] - t_z0['F675W'], bins = 1000)
#plt.plot(t_z0['F458W'], t_z0['F675W'], '.')
plt.show()
#savemyplot("Odds_vs_F814W")
plt.clf()


print t_z0['F458W'] - t_z0['F675W']




