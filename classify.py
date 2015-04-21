#! /usr/local/bin python
# -*- coding: iso-8859-1 -*-

# Importing packages
import numpy as np
import matplotlib.pylab as plt
from astropy.table import Table, join
#import astropy.cosmology
from sklearn.linear_model import LogisticRegression
import sys
import os
import glob
from astropy.io import ascii

# Catalog definitions:
# 'stell' corresponds to SExtractor stellaricity
# 'Stellar_Flag' is a much better indicator that takes in account many parameters: F814W, FWHM, J-Ks, F489W-F814W. 
#    See Molino's paper, section "Star/galaxy separation".



create_master_catalog = False
master_catalog_name = "alhambra.Master.ALLDATA.cat"
master_catalog_directory = "data/catalogs/master_catalogs/"
# Creates the Master catalog here:
if create_master_catalog:

	if not os.path.exists(master_catalog_directory):
		os.makedirs(master_catalog_directory)

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


		# Selects only the objects that are not duplicated, or, if they are, selects only the best SNR object.
		t = t[:][np.where(t['DupliDet_Flag']==0)]
		print "There are "+ str(len(t))+" unique elements when low SNR duplicates are deleted"


		# Joins the new catalog to the master catalog being created
		if catalog_filename != catalog_names[0]:
			print t['Stell_Mass_1'].mean()
			t = join(t_buffer, t, join_type='outer')
			print t['Stel			l_Mass_1'].mean()

	# Writes the Master catalog:
	ascii.write(t, master_catalog_directory + master_catalog_name)
	np.save(master_catalog_directory + master_catalog_name + ".npy", t)


create_arrays_long = True
if create_arrays_long is True: 

	number_lines_to_read = 100000
	t = Table.read("data/catalogs/master_catalogs/alhambra.Master.ALLDATA.cat", format='ascii', data_end=number_lines_to_read)

	y_stell_sex = t["stell"]
	y_stell_mol = t["Stellar_Flag"]

	t.remove_column('stell')
	t.remove_column('Stellar_Flag')
	t.remove_column('ID')

	X = []
	for i in t:
		X.append(list(i))
	X = np.array(X)

	np.save("ml/master_np_X.npy", X)
	np.save("ml/master_np_y_stell_sex.npy", y_stell_sex)
	np.save("ml/master_np_y_stell_mol.npy", y_stell_mol)


X = np.load("ml/master_np_X.npy") # the matrix
ymol = np.load("ml/master_np_y_stell_mol.npy") # The stellaricity from Molino
ysex = np.load("ml/master_np_y_stell_sex.npy") # The stellaricity from Bertin's sextractor

limit = 0.1
center = 0.5
ymol_unsure_mask = np.where(np.abs(ymol-center)>center-limit)

ymol_sure = ymol[ymol_unsure_mask]
X_sure = X[ymol_unsure_mask,:]
X_sure = X_sure[0]

ymol_bool = []
for i in ymol_sure:
	if i<limit: ymol_bool.append(-1)
	if i>1-limit: ymol_bool.append(1)

#plt.figure(figsize=(10,10), facecolor='w', edgecolor='k')
#plt.hist(ymol[ymol_unsure_mask])
#plt.show()
#plt.ion()
#plt.close()

# Will learn from the even positions:
X_sure_learn = X_sure[0::2]
ymol_learn = ymol_bool[0::2]

# Sill predict the odd positions
X_sure_target = X_sure[1::2]
y_target = ymol_bool[1::2]

# ML
model = LogisticRegression(C=1.)
model.fit(X_sure_learn, ymol_learn)
y_prediction = model.predict(X_sure_target)

# Results:
print "total number of objects to predict: "+ str(len(y_target))
N_correct = len(np.where(y_prediction - y_target == 0)[0])
print "I have predicted correctly "+str(N_correct)+" objects"
N_error = len(np.where(y_prediction - y_target != 0)[0])
print "there are " + str(N_error) + " errors (" +str(float(N_error)/float(N_correct)*100.)+"%)"


