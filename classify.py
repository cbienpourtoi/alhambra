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
from astroML.plotting import hist
from sklearn.decomposition import PCA
from sklearn import preprocessing
import machine
from sklearn.metrics import confusion_matrix

# Catalog definitions:
# 'stell' corresponds to SExtractor stellaricity
# 'Stellar_Flag' is a much better indicator that takes in account many parameters: F814W, FWHM, J-Ks, F489W-F814W. 
#    See Molino's paper, section "Star/galaxy separation".


def table_cleaner(table):

    table.remove_column('ID')
    table.remove_column('Field')
    table.remove_column('Pointing')
    table.remove_column('CCD')
    table.remove_column('RA')
    table.remove_column('Dec')
    table.remove_column('x')
    table.remove_column('y')
    table.remove_column('photoflag')
    table.remove_column('nfobs')
    table.remove_column('Satur_Flag')
    table.remove_column('DupliDet_Flag')

    return table



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
            print t['Stel            l_Mass_1'].mean()

    # Writes the Master catalog:
    ascii.write(t, master_catalog_directory + master_catalog_name)
    np.save(master_catalog_directory + master_catalog_name + ".npy", t)


create_arrays_long = True
if create_arrays_long is True: 

    number_lines_to_read = 10000
    t = Table.read("data/catalogs/master_catalogs/alhambra.Master.ALLDATA.cat", format='ascii', data_end=number_lines_to_read)

    ########
    # Let's begin by looking a bit at the data
    
    F814W_limit = 22.5 # The limit over which the Stellar_Flag is always 0.5
    t_faint = t[np.where(t["F814W"]>=F814W_limit)]
    t_bright = t[np.where(t["F814W"]<F814W_limit)]
    if np.where(t_faint["Stellar_Flag"] != 0.5)[0].size == 0: print "Sanity check OK: all faint sources have Stellar_Flag = 0.5"
    if np.where(t_bright["Stellar_Flag"] == 0.5)[0].size != 0: print "Some bright sources also have Stellar_Flag = 0.5"
    # Shows that some random values of the bright sample can also have Stellar_Flag = 0.5 (which is OK)
    # print t_bright[np.where(t_bright["Stellar_Flag"] == 0.5)[0]]["F814W"] # Shows teh values of F814W for these bright objects with Stellar_Flag = 0.5

    plt.figure(figsize=(10,10), facecolor='w', edgecolor='k')
    plt.hist(t_bright["Stellar_Flag"], bins=100)
    plt.title("Bright objects Stellar_Flag distribution")
    #plt.show()
    plt.close()

    # Gets only the objects with Stellar_Flag<0.1 or Stellar_Flag>0.9: the "sure" catalog
    limit = 0.1
    center = 0.5
    t_bright_sure = t_bright[np.where(np.abs(t_bright["Stellar_Flag"]-center)>center-limit)]

    # Sanity check:
    plt.figure(figsize=(10,10), facecolor='w', edgecolor='k')
    plt.hist(t_bright_sure["Stellar_Flag"], bins=100)
    plt.title("Bright objects Stellar_Flag distribution, only for very certain objects")
    #plt.show()
    plt.close()

    # Selects stars or galaxies from the "sure" catalog
    t_bright_sure_star = t_bright_sure[np.where(t_bright_sure["Stellar_Flag"]>1.-limit)]
    t_bright_sure_gal = t_bright_sure[np.where(t_bright_sure["Stellar_Flag"]<limit)]

    # Sanity check:
    plt.figure(figsize=(10,10), facecolor='w', edgecolor='k')
    plt.hist(t_bright_sure_star["Stellar_Flag"], label = "star")
    plt.hist(t_bright_sure_gal["Stellar_Flag"], label = "galaxy")
    plt.title("Bright objects Stellar_Flag distribution for very certain stars vs galaxies")
    plt.legend()
    #plt.show()
    plt.close()

    # Plots the hists of the different variables:
    make_all_plots = False
    if make_all_plots:
        for c in t.columns:
            plt.figure(figsize=(10,10), facecolor='w', edgecolor='k')
            plt.title("Distribution variable: "+c)
            hist(t_bright_sure_star[c], label = "star", alpha = 0.5, bins='blocks', normed = True)
            hist(t_bright_sure_gal[c], label = "gal", alpha = 0.5, bins='blocks', normed = True)
            plt.legend()
            plt.show()
            plt.close()
        







    t_bright_sure = table_cleaner(t_bright_sure)

    X = []
    for t2 in t_bright_sure:
        X.append(list(t2))
    X = np.array(X)

    """

    #Normalization:
    scaler = preprocessing.StandardScaler().fit(X)
    Xnorm = scaler.transform(X)


    #
    # PCA
    #

    print "begins PCA analysis"
    pca = PCA()
    pca.fit(Xnorm)
    X_new = pca.transform(Xnorm)

    plt.figure()
    plt.plot(X_new[:, 0], X_new[:, 1], 'x')
    plt.xlabel('column 0')
    plt.ylabel('column 1')
    plt.title('column 0 x 1')

    plt.figure()
    plt.plot(X_new[:, 1], X_new[:, 2], 'x')
    plt.xlabel('column 1')
    plt.ylabel('column 2')
    plt.title('column 1 x 2')


    plt.figure()
    plt.plot(pca.components_[0,:], 'r')
    plt.plot(pca.components_[1,:], 'g')
    plt.plot(pca.components_[2,:], 'b')
    plt.title('loadings 0..2')

    plt.show()

    print(pca.explained_variance_ratio_)

    
    sys.exit()







    y_stell_sex = t["stell"]
    y_stell_mol = t["Stellar_Flag"]

    t.remove_column('stell')
    t.remove_column('Stellar_Flag')
    t.remove_column('ID')

    X = []
    for i in t:
        X.append(list(i))
    X = np.array(X)

    zb_1 = t["zb_1"]
    Stell_Mass_1 = t["Stell_Mass_1"]

    np.save("ml/master_np_X.npy", X)
    np.save("ml/master_np_y_stell_sex.npy", y_stell_sex)
    np.save("ml/master_np_y_stell_mol.npy", y_stell_mol)
    np.save("ml/master_zb_1.npy", zb_1)
    np.save("ml/master_Stell_Mass_1.npy", Stell_Mass_1)

"""
X = np.load("ml/master_np_X.npy") # the matrix
ymol = np.load("ml/master_np_y_stell_mol.npy") # The stellaricity from Molino
ysex = np.load("ml/master_np_y_stell_sex.npy") # The stellaricity from Bertin's sextractor
zb_1 = np.load("ml/master_zb_1.npy") # z
Stell_Mass_1 = np.load("ml/master_Stell_Mass_1.npy") # stellar mass



########
# Test on the known part of the sample:
# We try only on "good" objects: Stellar_Flag>0.9 and Stellar_Flag<0.1
# half of this sample is the learning sample
# half is the target
# Target is predicted with 1% error (for a sample of 1000 "good" objects)


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


model = machine.logistic()
model.train(X_sure_learn, ymol_learn)
y_prediction = model.use(X_sure_target)


confm = confusion_matrix(np.array(y_target), y_prediction, labels=[1, -1])
print confm[0,0]/float(np.sum(confm[0,:]))
print confm[1,1]/float(np.sum(confm[1,:]))



sys.exit()


"""
# ML
model = LogisticRegression(C=1.)
model.fit(X_sure_learn, ymol_learn)
y_prediction = model.predict(X_sure_target)
"""

# Results:
print "total number of objects to predict: "+ str(len(y_target))
N_correct = len(np.where(y_prediction - y_target == 0)[0])
print "I have predicted correctly "+str(N_correct)+" objects"
error_places = np.where(y_prediction - y_target != 0)[0]
print error_places
N_error = len(error_places)
print "there are " + str(N_error) + " errors (" +str(float(N_error)/float(N_correct)*100.)+"%)"


smsure = Stell_Mass_1[ymol_unsure_mask]
smsure_target = smsure[1::2]
print smsure_target[error_places]




Stell_Mass_1_z08 = Stell_Mass_1[np.where(np.abs(zb_1-0.8)<0.05)[0]]



tz = t[np.where(np.abs(t["zb_1"]-0.8)<0.05)[0]]
plt.hist(tz['Stell_Mass_1'], bins = 100)
plt.hist(Stell_Mass_1_z08_sure, bins = 100, alpha = 0.6)
plt.show()

sys.exit()

########
# Now we use "good" objects to predict Stellar_Flag = 0.5 objects:

model = LogisticRegression(C=1.)
model.fit(X_sure, ymol_bool) # Fit on all objects where we know the stellaricity very well (Stellar_Flag>0.9 and Stellar_Flag<0.1)

faint = np.where(ymol==0.5)[0]
print "Total number of objects: " + str(len(ymol))
print "Number of faint objects: " + str(len(faint))
X_faint = X[faint, :]

y_faint_target = model.predict(X_faint)
faint_star = np.where(y_faint_target == 1)[0]
faint_gal = np.where(y_faint_target == -1)[0]

print len(faint_star)
print len(faint_gal)

#plt.hist()
