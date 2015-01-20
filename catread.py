''' File: catread.py '''
''' Copyright: Loic Le Tiran, 2014 '''
''' Contact: loic.le-tiran@obspm.fr '''
''' Licence: GNU GPL v3 '''

''' Description: '''
''' Tool for reading catalogs from the Alhambra Survey '''
''' These catalogues are hosted on: https://cloud.iaa.csic.es/alhambra/catalogues/ '''
''' To use the catalogues, see Molino 2013 (soon 2014): http://arxiv.org/abs/1306.4968 (pdf: http://arxiv.org/pdf/1306.4968v2) '''


import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, join
from astropy.io import ascii
#from astropy import cosmology # Astropy version 0.3
from astropy.cosmology import Planck13 as cosmo # Astropy version 0.4
from astropy import units
import scipy
from matplotlib.colors import LogNorm
import sys
import glob
import os


plot_extension = ".png"
plot_directory = "plots/"
def savemyplot(name):
	fig.savefig(plot_directory+name+plot_extension)
	return



# Make a table of the spectral axis of the Alhambra catalog
filter_names = ['F365W', 'F396W', 'F427W', 'F458W', 'F489W', 'F520W', 'F551W', 'F582W', 'F613W', 'F644W', 'F675W', 'F706W', 'F737W', 'F768W', 'F799W', 'F830W', 'F861W', 'F892W', 'F923W', 'F954W', 'J', 'H', 'KS']
filter_leff = np.array([365., 396., 427., 458., 489., 520., 551., 582., 613., 644., 675., 706., 737., 768., 799., 830., 861., 892., 923., 954., 1216., 1655., 2146.])
filters = Table([filter_names,filter_leff], names=('Filter', 'Lambda eff (A)'), meta={'name': 'table of the filters'})



# Test mode ? If true, only a small number of values will be red from the master catalog
test_mode = False

if test_mode: 
	number_lines_to_read = 1000
else:
	number_lines_to_read = -1


# Will create a Master catalog from the alhambra catalogs if True
# If False, will use the master catalog file created previously.
create_master_catalog = True


# Name of Master catalog (containing all sorted infos from all alhambra catalogs):
master_catalog_name = "alhambra.Master.ColorProBPZ.cat"
master_catalog_directory = "data/catalogs/master_catalogs/"

# Creates the Master catalog here:
if create_master_catalog:

	if not os.path.exists(master_catalog_directory):
		os.makedirs(master_catalog_directory)

	# If True, the code will study some parameters and make some plots. Useless usually.
	study_stellaricity = False
	study_Flags = False
	study_Odds = False
	study_z = False


	# Values of the cuts: # from mail Alberto 24/7/14
	Stellar_Flag_Cut = 0.7
	Odds_1_Cut = 0.2
	irms_OPT_Flag_Cut = 10
	irms_NIR_Flag_Cut = 2
	photoflag_Cut = 1
	Chi2_Cut = 500

	# OLD Values of the cuts:
	#Stellar_Flag_Cut = 0.7
	#Odds_1_Cut = 0.2
	#irms_OPT_Flag_Cut = 3
	#irms_NIR_Flag_Cut = 1
	#photoflag_Cut = 1


	# Expected number of Alhambra catalogues
	exp_numb_cats = 48


	# NB the 6/8/14: This bug has been corrected since astropy 0.4
	# Due to a bug in astropy.table (see :
	# http://stackoverflow.com/questions/22617428/overflowerror-python-int-too-large-to-convert-to-c-long-with-astropy-table
	# and
	# https://github.com/astropy/astropy/pull/2234 )
	# I truncate the first line of my catalogs until the bug is corrected.
	# Truncated catalog name has a "first_col_truncated." prefix, and is in 
	# a different directory "/first_col_truncated/".
	#catalogs_directory = 'data/catalogs/first_col_truncated/first_col_truncated.' #this line contains the repertory AND the prefix.
	# When the bug is solved,use this line:
	catalogs_directory = 'data/catalogs/original_catalogs/'

	catalog_names = glob.glob(catalogs_directory+'alhambra.*.ColorProBPZ.cat')
	print "I have found "+str(len(catalog_names))+" catalogues."
	if exp_numb_cats != len(catalog_names) :
		print "WARNING: The expected number of catalogues is "+str(exp_numb_cats)+". Why do we read a different number?" 

	'''
	t = Table.read(catalog_names[0], format='ascii')
	t_buffer = Table.read(catalog_names[1], format='ascii')
	t2 = join(t,t_buffer)

	print t2
	print len(t2)
	'''


	i=0
	for catalog_filename in catalog_names:
		if catalog_filename == catalog_names[0]:
			t = Table.read(catalog_filename, format='ascii')
			print "first catalog"
		else:
			t_buffer = t
			t = Table.read(catalog_filename, format='ascii')
			#t = join(t_buffer, t, join_type='outer')
		i=i+1
		print "I have red "+ catalog_filename +" , # "+ str(i)

	#ascii.write(t, 'test.dat')

	#sys.exit()


	#= catalogs_directory+'alhambra.F02P01C01.ColorProBPZ.cat'


	#t = Table.read(catalog_filename, format='ascii')



	#print "There are "+ str(len(t))+" elements in this catalogue"



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

	
		
		# Chi2 from mail Alberto 24/7/14
		t = t[:][np.where(t['Chi2']<Chi2_Cut)]
		print "Remains   "+ str(len(t))+" elements with Chi2 being less than "+str(Chi2_Cut)

		


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


		# Joins the new catalog to the master catalog being created
		if catalog_filename != catalog_names[0]:
			print t['Stell_Mass_1'].mean()
			t = join(t_buffer, t, join_type='outer')
			print t['Stell_Mass_1'].mean()

	# Writes the Master catalog:
	ascii.write(t, master_catalog_directory + master_catalog_name)
	np.save(master_catalog_directory + master_catalog_name + ".npy", t)
		
# Opens the Master catalog if we are not in master catalog creation mode:
else :
	
	if os.path.exists(master_catalog_directory + master_catalog_name):
		t = Table.read(master_catalog_directory + master_catalog_name, format='ascii', data_end=number_lines_to_read)
		print "Finished reading Master catalog: "+str(master_catalog_directory + master_catalog_name)
		print str(len(t))+ " lines red."
	else:
		print "Master catalog not found in "+str(master_catalog_directory + master_catalog_name)
		sys.exit()



# Redshift distribution whole sample
fig = plt.figure()
plt.title("z distribution")
plt.xlabel("z")
plt.hist(t['zb_1'], bins = 100)
#plt.show()
savemyplot("z_distribution")
plt.close()

# Stellar Mass distribution whole sample
fig = plt.figure()
plt.title("Stellar Mass distribution")
plt.xlabel("Stell_Mass_1")
plt.hist(t['Stell_Mass_1'], bins = 100)
#plt.show()
savemyplot("Stell_Mass_distribution")
plt.close()




t_keep = t;
i_mass = 0
ntot = 0.
nlowmass = 0.

# Mean z
for z_mean in (np.arange(15)+1.)/10.:
	print z_mean
	#z_mean = 0.3
	z_delta = 0.05

	#print t['Stell_Mass_1'].mean()
	t = t_keep[:][np.where(t_keep['zb_1']>z_mean-z_delta)]
	t = t[:][np.where(t['zb_1']<z_mean+z_delta)]
	#print t['Stell_Mass_1'].mean()
	#t = t[:][np.where(t['zb_1']>0.5)]


	#print "Remains   "+ str(len(t))+" elements with z<0.2"

	# Observed wavelengths for the color, in nm
	lem1 = 350. #nm
	lem2 = 700. #nm

	# Emission wavelengths for the color, in nm
	lobs1 = lem1*(1+z_mean)
	lobs2 = lem2*(1+z_mean)

	filter1 = filter_names[np.argmin(abs(filter_leff - lobs1))]
	filter2 = filter_names[np.argmin(abs(filter_leff - lobs2))]

	print filter1, filter2

	#filter_names = ['F365W']
	#filter_leff = [365.]
	#filters = Table()



	#l1 = 'F551W'
	#l2 = 'F954W'
	#l2 = 'J'

	t = t[:][np.where(t[filter1]<99)]
	t = t[:][np.where(t[filter2]<99)]

	'''
	# Redshifts distribution
	fig = plt.figure()
	#plt.title("Redshifts distribution")
	#plt.xlabel("Redshifts (zb - Bayesian)")
	#plt.ylabel("#")
	plt.hist(t['zb_1'], bins=10)
	plt.show()
	#savemyplot("")
	plt.clf()
	'''


	# Color vs Stellar Mass
	fig = plt.figure()
	plt.title("Color vs Stellar Mass @ z~"+str(z_mean))
	plt.xlabel("Stellar Mass")
	plt.ylabel(filter1+" - "+filter2)
	plt.xlim([6,12])
	plt.ylim([-1,4])
	#plt.hist(t['F675W'], bins = 50)
	plt.hist2d(t['Stell_Mass_1'], t[filter1] - t[filter2], bins = 200, range=np.array([(6, 12), (-1, 4)]))
	#plt.hist(t[l1] - t[l2], bins = 200)
	#plt.hist(t['Stell_Mass_1'], bins = 50)
	#plt.plot(t['F458W'], t['F675W'], '.')
	#plt.show()
	savemyplot("Color_StellarMass_"+str(z_mean))
	plt.close()

	# Color vs Apparent Magnitude
	fig = plt.figure()
	plt.title("Color vs Apparent Magnitude @ z~"+str(z_mean))
	plt.xlabel("m "+filter2)
	plt.ylabel(filter1+" - "+filter2)
	#plt.xlim([27,19])
	#plt.ylim([-1,4])
	#plt.hist(t['F675W'], bins = 50)
	plt.hist2d(t[filter2], t[filter1] - t[filter2], bins = 50, range=np.array([(19, 27), (-1, 4)]), norm=LogNorm())
	plt.axis([27, 19, -1, 4])
	#plt.show()
	savemyplot("Color_ApparentMagnitude_"+str(z_mean))
	plt.close()

	
	#DL = cosmology.luminosity_distance(t['zb_1'], cosmology.Planck13).to(units.pc)/units.pc
	DL = cosmo.luminosity_distance(t['zb_1']).to(units.pc)/units.pc
	absmag = t[filter2] - 5.*np.log10(DL) + 5.

	# Color vs Absolute Magnitude
	fig = plt.figure()
	plt.title("Color vs Absolute Magnitude @ z~"+str(z_mean))
	plt.xlabel("M "+filter2)
	plt.ylabel(filter1+" - "+filter2)
	plt.hist2d(absmag, t[filter1] - t[filter2], bins = 50, range=np.array([(-26, -12), (-1, 4)]), norm=LogNorm())
	plt.axis([-12, -26, -1, 4])
	#plt.show()
	savemyplot("Color_AbsoluteMagnitude_"+str(z_mean))
	plt.close()
	

	# Stellar Mass distribution per sub-sample
	fig = plt.figure()
	plt.title("Stellar Mass distribution for sub-sample z_mean="+str(z_mean))
	plt.xlabel("Stell_Mass_1")
	plt.hist(t['Stell_Mass_1'], bins = 100)
	#plt.show()
	savemyplot("Stell_Mass_distribution_"+str(z_mean))
	plt.close()


	# Stellar_Flag vs Stellar Mass
	fig = plt.figure()
	plt.title("Stellar_Flag vs Stellar Mass @ z~"+str(z_mean))
	plt.xlabel("Stellar Mass")
	plt.ylabel("Stellar_Flag")
	#plt.xlim([6,12])
	#plt.ylim([-1,4])
	#plt.hist(t['F675W'], bins = 50)
	plt.hist2d(t['Stell_Mass_1'], t['Stellar_Flag'], bins = 50, norm=LogNorm(), range=np.array([(6, 12), (-0.05, 0.8)]))
	#plt.hist(t[l1] - t[l2], bins = 200)
	#plt.hist(t['Stell_Mass_1'], bins = 50)
	#plt.plot(t['F458W'], t['F675W'], '.')
	#plt.show()
	savemyplot("Stellar_Flag_StellarMass_"+str(z_mean))
	plt.close()


	smass_limits = [6.2, 7.2, 7.5, 7.7, 8., 8.2, 8.4, 8.4,8.6, 8.7, 8.8, 8.9, 9., 9., 9.2]
	thislimit = smass_limits[i_mass]
	
	t_temp=t
	t = t[:][np.where(t['Stellar_Flag']==0.5)]
	# Stellar Mass distribution per sub-sample for Stellar_Flag = 0.5
	fig = plt.figure()
	plt.title("Stellar Mass distribution for Stellar_Flag=0.5 sub-sample z_mean="+str(z_mean))
	plt.xlabel("Stell_Mass_1")
	plt.hist(t['Stell_Mass_1'], bins = 100)
	plt.plot([thislimit,thislimit], [0,10], '-')
	#plt.show()
	savemyplot("Stell_Mass_distribution_Stellar_Flag_05_"+str(z_mean))
	plt.close()
	
	print "There are "+str( len(t[:][np.where(t['Stell_Mass_1']<thislimit)]) ) +" objects with Stell_Mass_1 <"+str(thislimit)+" at z="+str(z_mean)
	print "There are "+str( len(t[:][np.where(t['Stell_Mass_1']>thislimit)]) ) +" objects with Stell_Mass_1 >"+str(thislimit)+" at z="+str(z_mean)

	ntot = ntot + len(t[:][np.where(t['Stell_Mass_1']>thislimit)])
	nlowmass = nlowmass + len(t[:][np.where(t['Stell_Mass_1']<thislimit)])
	
	print "ratio for all z<"+str(z_mean)+" : "+str(nlowmass/ntot)


	# Odds vs Stellar Mass
	fig = plt.figure()
	plt.title("Odds vs Stellar Mass @ z~"+str(z_mean))
	plt.xlabel("Stellar Mass")
	plt.ylabel("Odds")
	#plt.xlim([6,12])
	#plt.ylim([-1,4])
	#plt.hist(t['F675W'], bins = 50)
	plt.hist2d(t['Stell_Mass_1'], t['Odds_1'], bins = 50, norm=LogNorm(), range=np.array([(6, 12), (-0.05, 0.8)]))
	#plt.hist(t[l1] - t[l2], bins = 200)
	#plt.hist(t['Stell_Mass_1'], bins = 50)
	#plt.plot(t['F458W'], t['F675W'], '.')
	#plt.show()
	savemyplot("Odds_StellarMass_"+str(z_mean))
	plt.close()

	
	i_mass = i_mass + 1


	t=t_temp
	t = t[:][np.where(t['Stellar_Flag']!=0.5)]
	# Stellar Mass distribution per sub-sample for Stellar_Flag != 0.5
	fig = plt.figure()
	plt.title("Stellar Mass distribution for Stellar_Flag!=0.5 sub-sample z_mean="+str(z_mean))
	plt.xlabel("Stell_Mass_1")
	plt.hist(t['Stell_Mass_1'], bins = 100)
	#plt.show()
	savemyplot("Stell_Mass_distribution_Stellar_Flag_not05_"+str(z_mean))
	plt.close()



	'''
	# z distribution sub-sample
	fig = plt.figure()
	plt.title("z distribution for sub-sample z_mean="+str(z_mean))
	plt.xlabel("z")
	plt.hist(t['zb_1'], bins = 20)
	#plt.show()
	savemyplot("z_distribution_"+str(z_mean))
	plt.close()
	'''


	
	
	

sys.exit()

nspec = 10
magspec = np.array(t[filters['Filter'][:].tolist()][nspec]).tolist()
spec = 10.**(-2.5 * np.array(magspec))

'''
# Plot spectra
fig = plt.figure()
#plt.title("Odds VS F814W")
#plt.xlabel("F814W magnitude")
#plt.ylabel("Odds")
#plt.hist(t['F675W'], bins = 50)
#plt.plot(filters['Lambda eff (A)'][:], t[filters['Filter'][:].tolist()][0])
for nspec in range(0,30):
	magspec = np.array(t[filters['Filter'][:].tolist()][nspec]).tolist()
	spec = 10.**(-2.5 * np.array(magspec))
	spec = spec/(spec[9])
	plt.plot(filters['Lambda eff (A)'][:], spec, '-')
plt.show()
#savemyplot("Odds_vs_F814W")
plt.close()
'''





