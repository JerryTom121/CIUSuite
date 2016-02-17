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
import sys
import numpy as np
import matplotlib.pyplot as plt
import csv
import scipy.stats as stats
### Functions ####
def getCIUData(fname):
	fname = str(fname)
	rawData=np.genfromtxt(fname, missing_values=[""], filling_values=[0], delimiter=",")
	return rawData, fname

# get axes for a matrix
def getAxisParameters(rawData):
	yaxis = rawData[1:,0]
	xaxis = rawData[0,1:]
	return xaxis, yaxis

# Normalize a matrix across drift bins
def CIUnormalize(rawData, fileName):
	strippedMat=rawData[1:,1:]
	maxInt = np.amax(strippedMat, axis = 0)
	norm = strippedMat/maxInt
	return norm

# crop a matrix
def crop(EL, EH, DL, DH, AP, file1n, file2n):
	EL = np.where(AP[0] == EL)
	EL = int(EL[0])
	EH = np.where(AP[0] == EH)
	EH = int(EH[0])
	DL = find_nearest(AP[1], DL)
	DH = find_nearest(AP[1], DH)
	file1n = file1n[DL:DH+1]
	f1 = file1n[:,EL:EH+1]
	file2n = file2n[DL:DH+1]
	f2 = file2n[:,EL:EH+1]
	EN = AP[0][EL:EH+1]
	DT = AP[1][DL:DH+1]
	return f1, f2, EN, DT

# Find the matrix value closest to a user input
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin() ###get the index of the value nearest to the input value
    return idx

# subtract matrices and comput the scaled RMSD
def difference(f1, f2):
	f1[f1 < 0.1] = 0
	f1v = np.count_nonzero(f1)
	f2[f2 < 0.1]= 0
	f2v = np.count_nonzero(f2)
	dif = f1 - f2
	RMSD = ((np.sum(dif**2)/(f1v + f2v))**0.5)*100
	return dif, RMSD

# Make a plot of the difference matrix
def plot(title, mat,EN, DT, v1, v2):
	plt.clf()
	plt.title(title)
	plt.contourf(EN,DT,mat, v1, cmap = "bwr", ticks = "none")
	plt.tick_params(axis = 'x', which = 'both', bottom = 'off', top = 'off', left = 'off', right = 'off')
	plt.tick_params(axis = 'y', which = 'both', bottom = 'off', top = 'off', left = 'off', right = 'off')
	plt.annotate(rtext, xy=(260,10), xycoords='axes points')
	plt.xlabel(titlex)    
	plt.ylabel(titley)
	plt.colorbar(ticks = v2)
	plt.savefig(title +'.png')
	plt.close()
    
###############################

####### Introduction and options ################
print("CIU Fingerprint(Comparison")
print("This tool will generate difference plots between CIU fingerprints. There are three modes: Basic, Batch, and Cluster.")
print("Basic mode requires two files as inputs, and will generate a plot of File1 - File2" + '\n' + "Batch mode will generate a plot of Ref - i, for all files i in the directory."+"\n"+"Cluster mode will calculate the RMSD for each pair of files, and use k-medoids clustering to sort them")
answer = input("Please enter basic or batch ")
lower = answer.lower()
assert lower == 'basic' or lower == 'batch' , "You must enter basic or batch"

def optionsmenu():
	print("Settings:")
	print("1. Cropping")
	print("2. Intensity Scaling")
	print("3. Axis Titles")
	print("4. Exit Menu")
	print("\n")
	option = input("Please choose an option by entering its numerical value : ")
	return option
    
def cropinputs():
	EL = int(input("Please enter the lowest collision energy you would like to scan: "))
	EH = int(input("Please enter the highest collision energy you would like to scan: "))
	DL = round(float(input("Please enter the lowest value of drift time you would like to scan(use one decimal place e.g. 5.0): ")),1)
	DH = round(float(input("please enter the highest value of drift time you would like to scan(use one decimal place e.g. 40.6): ")),1)
	print("\n")
	return EL, EH, DL, DH 


#### Default Params ####
titleyn = 'n'
cropyn = 'n'
titlex = "Trap Collision Energy (V)"
titley = "Drift Time (ms)" 
s = 1

#### Menu ####       

print("Default settings are no cropping, no scaling of intensities, and axes titles of Trap Collision Voltage(V) and Drift Time(ms)")
optionsyn = input("Please press y to continue or n to change these settings  ")
optionsyn = optionsyn.lower()
if optionsyn == 'n':
	while True:
		option = optionsmenu()
		if option == "1":
			cropyn = 'y'
			EL, EH, DL, DH = cropinputs() 
			print("Crop Inputs Saved, Please choose another option, or exit the menu")
			continue
          
		if option == "2":
			print("This is the intensity scale for the difference plot. Lower values will exaggerate differences, higher ones will make them appear more subtle") 
			s = input("Please enter a value from 0 to 1: ")
			s = float(s)
			continue
		if option == "3":
			title = input("Plot title: ")
			titleyn = "y"
			titlex = input("X axis title: ")
			titley = input("Y axis title: ")
			continue
		if option == "4":
			break
		else:
			continue
            

################################################
################################################

if lower == "Basic" or lower == "basic":
	print("Basic Mode:")
    
	print("Please Select File1 and File2")
	files=os.listdir(".")
	files=[x for x in files if x.endswith('_raw.csv')]
	i = 1
	for item in files:
		print(i, item)
		i = i +1
        
	fileid1 = input("File1 Number : ")
	file1 = files[int(fileid1) -1]
	fileid2 = input("File2 Number : ")
	file2 = files[int(fileid2) -1]

    ### Do the work ###    
	file1 = getCIUData(file1)
	file2 = getCIUData(file2) ## load the data
	file1n = CIUnormalize(file1[0],file1[1])
	file2n = CIUnormalize(file2[0],file2[1])
	AP = getAxisParameters(file1[0]) ## look at the axes
	if cropyn == 'n': ### If no cropping, use the whole matrix
		f1, f2, EN, DT = crop(AP[0][0],AP[0][-1],AP[1][0], AP[1][-1],AP,file1n, file2n)
	if cropyn == 'y':  ### If cropping, find the indices that correspond to the extrema that were input
		f1, f2, EN, DT = crop(EL, EH, DL, DH, AP, file1n, file2n)

	dif, RMSD = difference(f1, f2)

	rtext = "RMSD = " +  '%2.2f' % RMSD

	FN1 = str(file1[1])
	FN2 = str(file2[1])
	if titleyn != 'y':
		title = FN1[:-8] + "-" + FN2[:-8]

	v1 = np.linspace(-s,s, 100, endpoint=True)
	v2 = np.linspace(-s,s,6,endpoint=True)
	plot(title,dif,EN, DT, v1, v2)
	print("\n" + "Processing ", file1[1] +" - " + file2[1])
	print("%RMSD = " , '%2.2f' % RMSD )
    

######################################################################
######################################################################

if lower == "Batch" or lower == "batch":
	#subroutine 2
	print("Batch Mode:" + '\n' + "all files must end with '_raw.csv' to be read with Batch Mode")
	files=os.listdir(".")
	files=[x for x in files if x.endswith('_raw.csv')]
	print("Please select a reference file")
	for item in enumerate(files):
		print(item[0], item[1])
        
	fileid = input("File Number : ")
	ref = files[int(fileid) -1]
	del files[int(fileid) - 1]

	ref = getCIUData(ref)
	refn =CIUnormalize(ref[0],ref[1])  ##normalized reference file
	AP = getAxisParameters(ref[0])
	with open(ref[1].rstrip("_raw.csv") + " RMSDs.csv", "w") as csvfile:
		writer = csv.writer(csvfile, delimiter = ",")
		writer.writerow(["File Name","RMSD"]) 
		writer.writerow([""]) 
		for item in files:
			i = getCIUData(item)
			t = i[1] ## the name of the file for later
			inorm = CIUnormalize(i[0],i[1])  #normalized i file
			if cropyn == 'n':
				f1, f2, EN, DT = crop(AP[0][0],AP[0][-1],AP[1][0], AP[1][-1],AP,refn, inorm)
			if cropyn == 'y':
				f1, f2, EN, DT = crop(EL, EH, DL, DH, AP, refn, inorm)
                    
			dif, RMSD = difference(f1, f2)
			title1 = ref[1][:-8]
			title2 = t[:-8]
			title = title1 + " - " + title2
			v1 = np.linspace(-s,s, 100, endpoint=True)
			rtext = "RMSD = " +  '%2.2f' % RMSD
			v2 = np.linspace(-s,s,6,endpoint=True)
			plot(title,dif,EN, DT, v1, v2)
			writer.writerow([i[1],RMSD]) 
			print("\n" + "Processing ", ref[1].rstrip("_raw") +" - " + i[1])
			print("%RMSD = " ,  "%.2f" % RMSD)
		print("\n" + "A list of these values has been output to a CSV file for further analysis")
 
##CLUSTER HAS BEEN DEPRECATED##     
#############################################################
#if lower == "Cluster" or lower == "cluster":
#    from Pycluster import  kmedoids
#    files=os.listdir(".")
#    files=[x for x in files if x.endswith('_raw.csv')]
#    c = int(input("Please enter the number of clusters to search for: "))
#    rmsdmat = []
#    for item in files: 
#        item = getCIUData(item) 
#        inorm = CIUnormalize(item[0],item[1])
#        rmsdmat1 = []
#        for item1 in files:
#            item1 = getCIUData(item1)
#            i1norm = CIUnormalize(item1[0],item1[1])
#            AP = getAxisParameters(item1[0])
#            
#            if cropyn == 'n':
#                ELi = 0
#                EHi = len(AP[0])
#                DLi = 0
#                DHi = len(AP[1])
#            if cropyn == 'y':
#                ELi, EHi, DLi, DHi = getcroplimits(EL, EH, DL, DH, AP)
#                   
#            i1norm1 = i1norm[DLi:DHi+1]    
#            i1norm1 = np.swapaxes(i1norm1, 0, 1)
#            i1norm1 = i1norm1[ELi:EHi+1]
#            i1norm1 = np.swapaxes(i1norm1, 0, 1)
#            
#            inorm1 = inorm[DLi:DHi+1]   
#            inorm1 = np.swapaxes(inorm1, 0, 1)
#            inorm1 = inorm1[ELi:EHi+1]
#            inorm1 = np.swapaxes(inorm1, 0, 1)
#            file1n[file1n < 0.1] = 0
#    	    file1nvalues = np.count_nonzero(file1n)

#    	    file2n[file2n < 0.1]= 0
#              file2nvalues = np.count_nonzero(file2n)
#            difmat = inorm1[1] - i1norm1[1]     
#            difmat2 = difmat**2
#            RMSD = (np.average(difmat2))**0.5
#            pRMSD = RMSD*100
#            rmsdmat1.append(pRMSD)
#        rmsdmat.append(rmsdmat1)
#    rmsdmat = np.array(rmsdmat)    
#    rmsdmat = 1 - rmsdmat
#    np.savetxt("DifferenceMatrix.csv", rmsdmat, delimiter = ",")
#    clusterid, error, nfound = kmedoids (rmsdmat, nclusters=c,npass = 10)
#    print("Error = ", error)
#    print("Found this configuration " , nfound, " out of 10 times")
#    with open("Clusters.csv", "wb") as csvfile:
#        writer = csv.writer(csvfile, delimiter = ",")
#        writer.writerow(["item #","File Name","Cluster #"]) 
#        writer.writerow([""]) 
#        for i in range(0, len(clusterid)):
#            print(i, files[i], " cluster = ", clusterid[i])
#            writer.writerow([i,files[i], clusterid[i]])
#        writer.writerow([""])
#        writer.writerow(["Note : The cluster number is defined as the item number of the centroid of the cluster."])    
#            
#    print("\n" + "The cluster number is defined as the item number of the centroid of the cluster.")
#    print("\n" + "A list of the cluster assignments has been output to Clusters.csv.")
#    print("\n" + "If you would like to use a different clustering algorithm, a copy of the distance matrix has been output to DifferenceMatrix.csv")
#        

#else:
#    sys.exit()
