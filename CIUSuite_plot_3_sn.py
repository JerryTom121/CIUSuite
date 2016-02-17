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
import matplotlib.pyplot as plt
import scipy.signal as signal
import scipy.stats as stats
np.seterr(divide='ignore', invalid='ignore')

##### Introduction and Options ##########
print("Collision Induced Unfolding Fingerprint Generator")
print("all files must end with '_raw.csv'  to be read")
answer = input("Axis titles are in Trap Collision Voltage (V) and Drift Time (ms) by default." +"\n" +"Data is smoothed using a Savitsky-Golay Filter with a window length of 3 by default." + "\n" + "The data is not cropped by default" + "\n" + "Please type y to continue or n to change these settings: ")
lower = answer.lower() 
assert lower == 'y' or lower =='n', "Answer must be y or n" 

def optionsmenu():
    print("Settings:")
    print("1. Cropping")
    print("2. Savitsky Golay Smoothing")
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

titlex = "Trap Collision Voltage (V)"
titley = "Drift Time (ms)"
savyn = "y"
window = 3
cropyn = "n"

if lower == 'n':
    while True:
        a = optionsmenu()
        if a == "1":
            cropyn = 'y'
            EL, EH, DL, DH = cropinputs() 
            print("Crop Inputs Saved, Please choose another option, or exit the menu")
            continue
        if a == "2":
            savyn =input("Would you like to use the Savitsky Golay algorithm to smooth your data ? y or n : ") 
            if savyn.lower() == "y":
                window = int(input("Please enter the window width, must be an odd integer greater than or equal to 3: "))
                assert window % 2 != 0, "You must enter an odd integer"
                continue
        if a == "3":
                titlex = input("Please title your X axis:") #For our lab normally Trap Collision Voltage
                titley = input("Please title your Y axis:") #Normally Drift time or CCS
                continue
        if a == "4":
                break
            	    


#########################################    
    
    
##### Functions #########################
### Get data from raw_CSV files ###
def getCIUData(fname):
	fileName= str(fname)
	print(fileName)
	rawData=np.genfromtxt(fileName, missing_values=[""], filling_values=[0], delimiter=",")
	return rawData, fileName        ###### Returns Numpy Array of the Raw Data(with axes), followed by the filename

### Generate lists of trap collision energies and drift times used for the plots ###
def getAxisParameters(rawData):
	yaxis = rawData[1:,0]
	xaxis = rawData[0,1:]
	return xaxis, yaxis

### Normalize the files ###
def CIUnormalize(rawData):
	strippedMat=rawData[1:,1:] ### strip the axes off
	allmax = np.amax(strippedMat)
	maxInt = np.amax(strippedMat, axis = 0) ### find the maximum value for each energy
	norm = strippedMat/maxInt  ### Normalize to max value
	#plt.clf()
	#plt.scatter(np.linspace(0,len(maxInt), len(maxInt)),np.log10(maxInt/allmax))
	#plt.show()
	if savyn == 'y':	
		norm = signal.savgol_filter(norm, window,2, axis = 1)
	#np.savetxt(fileName.rstrip('_raw.csv') +"_norm.csv", normOut, delimiter=",")
	return norm  ### Output an array with and without axes

### Plot the CIU Fingerprint
def generateCIUPlot(data, AP, fname, crops, noise):
	plt.clf()  ### make a plot
	title = fname[:-8]  ### Generate a title based on the filename
	plt.title(title) ### plot the title
	#plt.axis([crops[0],crops[1],crops[2],crops[3]])

	plt.contourf(AP[0][crops[0]:crops[1]+1],AP[1][crops[2]:crops[3]+1],data,100,cmap='jet')  ### plot the data
	plt.ylabel(titley)
	plt.xlabel(titlex)    ### label the axes
	s = plt.twinx()
	#print(np.log10(1/np.array(noise)))
	plt.scatter(AP[0][crops[0]:crops[1]+1],np.log10(1/np.array(noise)),zorder = 1, c = "white", s = 100)
	axes = plt.gca()
	axes.set_xlim([AP[0][crops[0]],AP[0][crops[1]-1]])
	axes.set_xlim([AP[0][crops[0]],AP[0][crops[1]-1]])
	s.set_ylim([0,8.2])
	s.set_ylabel("log signal/noise")
	plt.savefig(fname[:-8]+'_sn.png') ### save the figure
	plt.close()  ### close the file
        
def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin() ###get the index of the value nearest to the input value
        return idx
        
def getcroplimits(EL, EH, DL, DH, AP):
	assert np.any(AP[0] == EL) == True, "The lowest collision energy you entered does not match any energies in your data"
	assert np.any(AP[0] == EH) == True, "The highest collision energy you entered does not match any energies in your data"
	EL = np.where(AP[0] == EL)
	EL = int(EL[0])
	EH = np.where(AP[0] == EH)
	EH = int(EH[0])
	DL = find_nearest(AP[1], DL)
	DH = find_nearest(AP[1], DH)
	return EL, EH, DL, DH

##########################################

#### Now do the actual work
files=os.listdir(".")
files=[x for x in files if x.endswith('_raw.csv')]  ### Make a list of files in the directory ending with _raw.csv
allData =[] ### Make an empty list to dump the data
golay = []
for item in files:
	OD = getCIUData(item)  ### Read the file
	Mat = OD[0]  ### Pull out the data
	FN = OD[1]  
	ND = CIUnormalize(Mat) ### Normalize
	AP = getAxisParameters(Mat) ### The axes
	if cropyn == 'n': ### If no cropping, use the whole matrix
            ELi = 0
            EHi = len(AP[0])
            DLi = 0
            DHi = len(AP[1])
	if cropyn == 'y':  ### If cropping, find the indices that correspond to the extrema that were input
	    ELi, EHi, DLi, DHi = getcroplimits(EL, EH, DL, DH, AP)
	ND1 = ND[DLi:DHi+1]    ### now crop the rows
	ND1 = ND1[:,ELi:EHi+1] ###crop the columns
	crops = ELi, EHi, DLi, DHi

	noises = []
	for i in ND1.T:
		split = np.array_split(i, 10)
		means = []
		for s in split:
			mean = np.mean(s)
			means.append(mean)
			#print(mean)
		
		noise = np.mean(sorted(means)[:3])
		print(noise)
		noises.append(noise)

	noises = np.array(noises)
	noises[noises == 0] = 1e-8
	noises[noises < 0] = 1e-8
	
	#sn.append(noises)
	#inoises = 1/noises
	#lnoises = np.log10(inoises)
	#print(noises)
	#print(inoises)
	#print(lnoises)
	allData.append((ND1, AP, FN, crops, noises)) ### dump the data without axes and the filename to the list above

for item in allData:
	generateCIUPlot(item[0], item[1], item[2], item[3], item[4])  #make the plot
	print("Processing", item[2])
	
	

