# Copyright (c) 2016 Joseph D. Eschweiler
#
#This file is part of CIUSUITE.

#CIUSUITE is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#CIUSUITE is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with CIUSUITE.  If not, see <http://www.gnu.org/licenses/>.

# Authors : Joseph D. Eschweiler; joeesch@umich.edu and Jessica Rabuck-Gibbons

###############################################################################################################
###############################################################################################################
#PLEASE CITE:

#CIUSuite: A Quantitative Analysis Package for Collision Induced Unfolding Measurements of Gas-Phase Protein Ions
#Joseph D. Eschweiler, Jessica N. Rabuck-Gibbons, Yuwei Tian, and Brandon T. Ruotolo*
#Anal. Chem., 2015, 87 (22), pp 11516–11522
#DOI: 10.1021/acs.analchem.5b03292
###############################################################################################################
###############################################################################################################

import os
import numpy as np
import csv
import scipy.stats as stats

##### Functions #########################
### Get data from raw_CSV files ###
def getCIUData(fname):
	fid=open(fname)
	fileName= str(fname)
	rawData=np.genfromtxt(fid, missing_values=[""], filling_values=[0], delimiter=",")
	return rawData, fileName        ###### Returns Numpy Array of the Raw Data(with axes), followed by the filename

### Generate lists of trap collision energies and drift times used for the plots ###
def getAxisParameters(rawData):
	yaxis = rawData[1:,0]
    	xaxis = rawData[0,1:]
	return xaxis, yaxis

### Normalize the files ###
def CIUnormalize(rawData, fileName):
	strippedMat=rawData[1:,1:] ### strip the axes off
	maxInt = np.amax(strippedMat, axis = 0) ### find the maximum value for each energy
	norm = strippedMat/maxInt  ### Normalize to max value
        return  norm  ### Output an array with and without axes
        
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin() ###get the index of the value nearest to the input value
    return idx

################################################################################



		
### Batch data into respective typeData arrays. Also store axis parameters for reconstructing graphs later ###
files=os.listdir(".")
files=[x for x in files if x.endswith('_raw.csv')]

#allData = []
typeIData = []
typeIIData = []
typeITitles = []
typeIITitles = []
axisParametersI = []
axisParametersII = []
UKData = []
UKTitles = []
axisParametersUK = []


######## now trim all the data to the given range in CV space ############
print "This program will give scores useful for classifying CIU Fingerprints."
answer = raw_input("If you would like to scan a specific range of collision energies, enter y" + "\n" +"If you would like to scan the whole spectra, press n:  ")
lower = answer.lower()
assert lower == 'y' or lower =='n', "you must enter y or n"
if lower == 'y':
    print "User Defined Energy Range"
    L = int(raw_input("Please enter the lowest collision energy you would like to scan: "))
    H = int(raw_input("Please enter the highest collision energy you would like to scan: "))
    for item in files:
	OD = getCIUData(item)
	Mat = OD[0]
	AP = getAxisParameters(Mat)
	assert np.any(AP[0] == L) == True, "The lowest collision energy you entered does not match any energies in your data"
	WL = np.where(AP[0] == L)
	WL = WL[0]
	WL = int(WL)
	assert np.any(AP[0] == H) == True, "The highest collision energy you entered does not match any energies in your data"
	WH = np.where(AP[0] == H)
	WH = int(WH[0])
	FN = OD[1]
	ND = CIUnormalize(Mat, FN)
	ND1 = np.swapaxes(ND, 0, 1)
	ND1 = ND1[WL:WH+1]
	ND1 = np.swapaxes(ND1, 0, 1)
	if FN.endswith('_typeI_raw.csv'):
		FN = FN[:-8]
		typeIData.append(ND1)
		typeITitles.append(FN)
		axisParametersI.append(AP)
	elif FN.endswith('_typeII_raw.csv'):
		FN = FN[:-8]
		typeIIData.append(ND1)
		typeIITitles.append(FN)
           	axisParametersII.append(AP)
        elif FN.endswith('_UK_raw.csv'):
		FN = FN[:-8]
		UKData.append(ND1)
		UKTitles.append(FN)
           	axisParametersUK.append(AP)
        
elif lower == 'n':
    print "Using the entire spectrum"
    for item in files:
	OD = getCIUData(item)
	Mat = OD[0]
	AP = getAxisParameters(Mat)
	FN = OD[1]
	ND = CIUnormalize(Mat, FN)
	if FN.endswith('_typeI_raw.csv'):
		FN = FN[:-8]
		typeIData.append(ND)
		typeITitles.append(FN)
		axisParametersI.append(AP)
	elif FN.endswith('_typeII_raw.csv'):
		FN = FN[:-8]
		typeIIData.append(ND)
		typeIITitles.append(FN)
		axisParametersII.append(AP)
	elif FN.endswith('_UK_raw.csv'):
		FN = FN[:-8]
	        UKData.append(ND)
		UKTitles.append(FN)
		axisParametersUK.append(AP)
		

    
############### NOW TRIM THE DATA IN DRIFT TIME SPACE ###########
typeIData = np.array(typeIData)
t1crop = []
t2crop = []
UKcrop = []
answer = raw_input("If you would like to scan specific drift time regions, press y" + "\n" +"If you would like to scan the entire drift time spectrum, press n: ")
lower = answer.lower()
assert lower == 'y' or lower == 'n', "you must enter y or n"
if lower == "y":
        lowdt = round(float(raw_input("Please enter the lowest value of drift time you would like to scan(use one decimal place e.g. 5.0): ")),1)
        assert type(lowdt) == float, "you must enter a number to one decimal place"
        highdt = round(float(raw_input("please enter the highest value of drift time you would like to scan(use one decimal place e.g. 40.6): ")),1)
        assert type(highdt) == float, "you must enter a number to one decimal place"
        assert np.any(np.around(axisParametersI[0][1],decimals = 1) == lowdt) == True, "The lowest drift time you entered does not correlate to any drift times in your data"
        assert np.any(np.around(axisParametersI[0][1],decimals = 1) == highdt) == True, "The highest drift time you entered does not correlate to any drift times in your data"
        L = np.where(np.around(axisParametersI[0][1],decimals = 1) == lowdt)
        H = np.where(np.around(axisParametersI[0][1],decimals = 1) == highdt)
        WL = L[0]
        WH = H[0]
        i = 0
        for item in typeIData:
            item1 = item[WL:WH+1]
            t1crop.append(item1)
        for item in typeIIData:
            item1 = item[WL:WH+1]
            t2crop.append(item1)
        for item in UKData:
            item1 = item[WL:WH+1]
            UKcrop.append(item1)
        typeIData = t1crop
        typeIIData = t2crop
        UKData = UKcrop
######################################################################################

### More functions that rely on the data above ####
### Calculate the average and standard dev. ###
def groupwise_statistics(group):
        mean = np.mean(group, axis=0)
        std = np.std(group, axis=0)
        return mean,std
        
### Calculate the Chi2 as a function of CV
def chisq2d(data, group, groupname):
	AV = groupwise_statistics(group)
	difmat = data - AV[0]
	difmat2 = difmat**2
	v = difmat2/AV[0]
	findnan = np.isnan(v)
	v[findnan]=0
	findinf = np.isinf(v)
	v[findinf]= 0
	v=np.array(v)
	#x = axisParametersI[0][0]
	colsum = np.sum(v, axis = 0)
	colsum = np.array(colsum)
	return colsum

#### Scoring METHOD 1 ####


print "METHOD 1" + "\n"
i =0
dsItoI = []
dsItoII = []
#### Calculate the Chi^2 matrices

with open("Scores.csv", "wb") as csvfile:
    writer = csv.writer(csvfile, delimiter = ",")
    writer.writerow(["File Name","TypeI Score", "Type II Score"]) 
    writer.writerow([""])
    print "Type I Scores" +"\n"
    ### Calculate scores for type Is by Method I (SEE README)###
    
    t12sl = []
    t11sl = []
    for data in typeIData:
           I = chisq2d(data, typeIData, typeITitles[i]+" vs typeI Average")
	   II = chisq2d(data, typeIIData, typeITitles[i]+" vs typeII Average")
	   print typeITitles[i]
	   print "Type  I Score = " , np.sum(I)
	   print "Type II Score = ", np.sum(II)
	   print I
	   avg1, std1 = groupwise_statistics(typeIData)
	   avg2, std2 = groupwise_statistics(typeIIData)
	   t1zs = (data - avg1)/std1*avg1
	   t1zs[np.isnan(t1zs)] = 0
	   t2zs = (data - avg2)/std2*avg2
	   t2zs[np.isnan(t2zs)] = 0
	   t1s = np.average(t1zs, axis = 0)
	   t2s = np.average(t2zs, axis = 0)
	   t1s = np.sum(t1s)
	   t2s = np.sum(t2s)
	   #writer.writerow([typeITitles[i], np.sum(I), np.sum(II), t1s, t2s])
	   #dsItoI.append(I)
	   #dsItoII.append(II)
	   t12sl.append(t2s)
	   t11sl.append(t1s)
	   i = i+1
    print "\n" + "Type II Scores" + "\n"
    writer.writerow([""])
    AItoI = np.average(dsItoI)
    AItoII = np.average(dsItoII)
    dsIItoI = []
    dsIItoII = []	


    ### Calculate scores for type II's by Method I (SEE README)###
    i = 0
    t21sl = []
    t22sl = []
    for data in typeIIData:
       	   I = chisq2d(data, typeIData, typeIITitles[i]+" vs typeI Average")
	   II = chisq2d(data, typeIIData, typeIITitles[i]+" vs typeII Average")
	   print typeIITitles[i]
	   print "Type  I Score = " , np.sum(I)
	   print "Type II Score = ", np.sum(II)
	   avg1, std1 = groupwise_statistics(typeIData)
	   avg2, std2 = groupwise_statistics(typeIIData)
	   t1zs = (data - avg1)/std1*avg1
	   t1zs[np.isnan(t1zs)] = 0
	   t2zs = (data - avg2)/std2*avg2
	   t2zs[np.isnan(t2zs)] = 0
	   t1s = np.average(t1zs, axis = 0)
	   t1s[np.isnan(t1s)] = 0
	   t2s = np.average(t2zs, axis = 0)
	   t1s = np.sum(t1s)
	   t2s = np.sum(t2s)
	   t21sl.append(t1s)
	   t22sl.append(t2s)
	   #writer.writerow([typeIITitles[i], np.sum(I), np.sum(II),t1s, t2s])
	   dsIItoI.append(I)
	   dsIItoII.append(II)
	   i = i+1
	   
    writer.writerow([""])	   
	   
    UKI = []
    UKII = []	

    print "\n" + "Uknown Scores" + "\n"
    ### Calculate scores for type Uknowns's by Method I (SEE README)###
    i = 0
    UK1 = []
    UK2 = []
    for data in UKData:
       	   I = chisq2d(data, typeIData, UKTitles[i]+" vs typeI Average")
	   II = chisq2d(data, typeIIData, UKTitles[i]+" vs typeII Average")
	   print UKTitles[i]
	   print i
	   print "Type  I Score = " , np.sum(I)
	   print "Type II Score = ", np.sum(II)
	   avg1, std1 = groupwise_statistics(typeIData)
	   avg2, std2 = groupwise_statistics(typeIIData)
	   t1zs = (data - avg1)/std1*avg1
	   t1zs[np.isnan(t1zs)] = 0
	   t2zs = (data - avg2)/std2*avg2
	   t2zs[np.isnan(t2zs)] = 0
	   t1s = np.average(t1zs, axis = 0)
	   t1s[np.isnan(t1s)] = 0
	   t2s = np.average(t2zs, axis = 0)
	   t1s = np.sum(t1s)
	   t2s = np.sum(t2s)
	   UK1.append(t1s)
	   UK2.append(t2s)
	   #writer.writerow([UKTitles[i], np.sum(I), np.sum(II),t1s, t2s])
	   UKI.append(I)
	   UKII.append(II)
	   i = i+1
	   
    t11av, t11std = groupwise_statistics(t11sl)
    t12av, t12std = groupwise_statistics(t12sl)
    t21av, t21std = groupwise_statistics(t21sl)
    t22av, t22std = groupwise_statistics(t22sl)
    
    i = 0

    writer.writerow(["","Zscore1", "Zscore2", "Pvalue1", "Pvalue2"])
    writer.writerow([""])
    writer.writerow(["Type I Data"])
    writer.writerow([""])
    for item in t12sl:
            zscore1 = (item - t12av)/t12std
            zscore2 = (item - t22av)/t22std
            pdf1 = stats.norm.sf(abs(zscore1))*2
            pdf2 = stats.norm.sf(abs(zscore2))*2
            
            s1 = pdf1
            s2 = pdf2
            
            
            
            writer.writerow([typeITitles[i], zscore1, zscore2,s1, s2])
            i = i+1    
    
    
    i = 0
    writer.writerow([""])
    writer.writerow(["Type II Data"])
    writer.writerow([""])

    for item in t22sl:
            zscore1 = (item - t12av)/t12std
            zscore2 = (item - t22av)/t22std
            pdf1 = stats.norm.sf(abs(zscore1))
            pdf2 = stats.norm.sf(abs(zscore2))
            
            s1 = pdf1
            s2 = pdf2
            
            
            
            writer.writerow([typeIITitles[i], zscore1, zscore2,s1, s2])
            i = i+1

    i = 0 
    writer.writerow([""])
    writer.writerow(["UKNOWN Data"])
    writer.writerow([""])
    for item in UK2:
            zscore1 = (item - t12av)/t12std
            zscore2 = (item - t22av)/t22std
            pdf1 = stats.norm.sf(abs(zscore1))
            pdf2 = stats.norm.sf(abs(zscore2))
            
            s1 = pdf1
            s2 = pdf2
            
            
            
            writer.writerow([UKTitles[i], zscore1, zscore2,s1, s2])
            i = i+1

csvfile.close()
            
    
#    writer.writerow([""])
#    writer.writerow(["Scoring Method 2"])
#    writer.writerow([""])



########## METHOD 2 ###################
#
#
# ##### SEE README ####
# #### First Acverage the Data for each type ###
#    print "\n" + "METHOD 2" + "\n"
#    adatalist1 = []
#    for data in typeIData:
#        adata1 = np.sum(data, axis = 1)*100
#        amax = np.max(adata1)
#        adata1 = adata1/amax
#        findnan = np.isnan(adata1)
#	adata1[findnan]=0
#	findinf = np.isinf(adata1)
#	adata1[findinf]= 0
#        adatalist1.append(adata1)
#    
#    adatalist2 = []
#    for data in typeIIData:
#        adata2 = np.sum(data, axis = 1)*100
#        amax = np.max(adata2)
#        adata2 = adata2/amax
#        findnan = np.isnan(adata2)
#	adata2[findnan]=0
#	findinf = np.isinf(adata2)
#	adata2[findinf]= 0
#        adatalist2.append(adata2)
#        
#    adatalistUK = []
#    for data in UKData:
#        adata2 = np.sum(data, axis = 1)*100
#        amax = np.max(adata2)
#        adata2 = adata2/amax
#        findnan = np.isnan(adata2)
#	adata2[findnan]=0
#	findinf = np.isinf(adata2)
#	adata2[findinf]= 0
#        adatalistUK.append(adata2)
#    
#    average1 = np.sum(adatalist1, axis = 0)
#    average2 = np.sum(adatalist2, axis = 0)
#    averageUK = np.sum(adatalistUK, axis = 0)
#    amax1 = np.max(average1)
#    amax2 = np.max(average2)
#    amaxUK = np.max(averageUK)
#    average2 = average2/amax2  ### normalize the averaged data
#    average1 = average1/amax1
#    averageUK = averageUK/amaxUK
#
#    i = 0
#    t1t1 = []
#    t1t2 = []
#    t2t1 = []
#    t2t2 = []
#    UK1 = []
#    UK2 = []
# ### Calculate the scores for type Is ###
#    print "Type I Scores"+ "\n"
#    for data in adatalist1:
#        difsq = (data - average1)**2
#        chisq11 = difsq/average1
#        chisq11[np.isnan(chisq11)] = 0
#        t11score = np.sum(chisq11)
#        difsq = (data - average2)**2
#        chisq12 = difsq/average2
#        chisq12[np.isnan(chisq12)] = 0
#        t12score = np.sum(chisq12)
#        t1t1.append(t11score)
#        t1t2.append(t12score)
#        print typeITitles[i]
#        print "Type  I Score = ", t11score
#        print "Type II Score = ", t12score
#        writer.writerow([typeITitles[i], t11score, t12score])
#        i = i + 1
#    writer.writerow([""])
#    
#    #### Calculate the Scores for type IIs ####
#    print "\n"+"Type II Scores" + "\n"
#    i = 0
#    for data in adatalist2:
#        difsq = (data - average1)**2
#        chisq21 = difsq/average1
#        chisq21[np.isnan(chisq21)] = 0
#        chisq21[np.isinf(chisq21)] = 0
#        t21score = np.sum(chisq21)
#        difsq = (data - average2)**2
#        chisq22 = difsq/average2
#        chisq22[np.isnan(chisq22)] = 0
#        chisq22[np.isinf(chisq22)] = 0
#        t22score = np.sum(chisq22)
#        print typeIITitles[i]
#        print "Type  I Score = ",t21score
#        print "Type II Score = ",t22score
#        writer.writerow([typeIITitles[i], t21score, t22score])
#        t2t1.append(t21score)
#        t2t2.append(t22score)
#        i = i + 1
#    writer.writerow([""])
#    print "\n"+"Uknown Scores" "\n"
#    
#    i = 0
#    for data in adatalistUK:
#        difsq = (data - average1)**2
#        chisq21 = difsq/average1
#        chisq21[np.isnan(chisq21)] = 0
#        chisq21[np.isinf(chisq21)] = 0
#        t21score = np.sum(chisq21)
#        difsq = (data - average2)**2
#        chisq22 = difsq/average2
#        chisq22[np.isnan(chisq22)] = 0
#        chisq22[np.isinf(chisq22)] = 0
#        t22score = np.sum(chisq22)
#        print UKTitles[i]
#        print "Type  I Score = ",t21score
#        print "Type II Score = ",t22score
#        writer.writerow([UKTitles[i], t21score, t22score])
#        UK1.append(t21score)
#        UK2.append(t22score)
#        i = i + 1
        
        
