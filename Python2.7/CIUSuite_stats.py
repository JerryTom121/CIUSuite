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


#### Functions #####
def optionsmenu():
    print "Settings:"
    print "1. Cropping"
    print "2. Savitsky Golay Smoothing"
    print "3. Axis Titles"
    print "4. Standard Deviation Intensity Scaling"
    print "5. Exit Menu"
    print "\n"
    option = raw_input("Please choose an option by entering its numerical value : ")
    return option
    
def cropinputs():
    EL = int(raw_input("Please enter the lowest collision energy you would like to scan: "))
    EH = int(raw_input("Please enter the highest collision energy you would like to scan: "))
    DL = round(float(raw_input("Please enter the lowest value of drift time you would like to scan(use one decimal place e.g. 5.0): ")),1)
    DH = round(float(raw_input("please enter the highest value of drift time you would like to scan(use one decimal place e.g. 40.6): ")),1)
    print "\n"
    return EL, EH, DL, DH 
    
    
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

def getCIUData(fname):
	fid=open(fname)
	fileName= str(fname)
	rawData=np.genfromtxt(fid, missing_values=[""], filling_values=[0], delimiter=",")
	return rawData, fileName

### Get the high and low values for the axes ###
def getAxisParameters(rawData):
	yaxis = rawData[1:,0]
    	xaxis = rawData[0,1:]
	return xaxis, yaxis

### Normalize the files ###
def CIUnormalize(rawData, fileName):
	strippedMat=rawData[1:,1:]
	maxInt = np.amax(strippedMat, axis = 0)
	norm = strippedMat/maxInt
        return  norm

### Plot the CIU Fingerprint
def generateCIUPlot(fname):
	OD = getCIUData(fname)
	AP = getAxisParameters(OD[0])
        plt.clf()
        plt.title(item[1])
	plt.contourf(AP[0],AP[1],item[0],100,cmap='jet')
        plt.xlabel(titlex)    
        plt.ylabel(titley)
        plt.colorbar()
        plt.savefig(item[1].rstrip('_raw.csv')+'.png')
        plt.close()





print '\n'
print "Collision Induced Unfolding Statistical Fingerprint Generator"
print "all files must end with '_raw.csv'  to be read"
answer = raw_input("The default titles of your average and standard deviation plots will be 'Average CIU Fingerprint' and 'Fingerprint Standard Deviation.'"+ '\n' + "Axis titles are in Trap Collision Energy (V) and Drift Time (ms) by default." + '\n' + "The standard deviation plots will be scaled to 0.6 relative to the raw data" + "\n" + "Please type y to continue or n to change these settings:")
lower = answer.lower() 
assert lower == "y" or lower == 'n', "You must enter y or n" 

savyn = "y"
window = 3
cropyn = "n"
s = 0.6
titleav = "Average CIU Fingerprint"
titlestd = "Fingerprint Standard Deviation"
titlex = "Trap Collision Energy (V)"
titley = "Drift Time (ms)"


if lower == 'n':
    while True:
        a = optionsmenu()
        if a == "1":
            cropyn = 'y'
            EL, EH, DL, DH = cropinputs() 
            print "Crop Inputs Saved, Please choose another option, or exit the menu"
            continue
        if a == "2":
            savyn =raw_input("Would you like to use the Savitsky Golay algorithm to smooth your data ? y or n : ") 
            if savyn.lower() == "y":
                window = raw_input("Please enter the window width, must be an odd integer greater than or equal to 3: ")
                assert window % 2 != 0, "You must enter an odd integer"
                continue
        if a == "3":
                titleav = raw_input("Please enter title for the average plot: ")
                titlestd = raw_input("Please enter title for the standard deviation plot: ")
                titlex = raw_input("Please title your X axis: ") #For our lab normally Trap Collision Energy
                titley = raw_input("Please title your Y axis: ") #Normally Drift time or CCS
                continue
        if a == "4":
                print "This is the intensity scale for the standard deviation plot. Lower values will exaggerate deviations, higher ones will make them appear more subtle" 
                s = raw_input("Please enter a value from 0 to 1: ")
                s = float(s)
                continue
        if a == "5":
                break
            	    


### Get data from raw_CSV files ###


files=os.listdir(".")
files=[x for x in files if x.endswith('_raw.csv')]
allData =[]
data = []
for item in files:
	OD = getCIUData(item)
	Mat = OD[0]
	FN = OD[1]
	ND = CIUnormalize(Mat, FN)
	AP = getAxisParameters(Mat)
	if cropyn == 'n': ### If no cropping, use the whole matrix
            ELi = 0
            EHi = len(AP[0])
            DLi = 0
            DHi = len(AP[1])
        if cropyn == 'y':  ### If cropping, find the indices that correspond to the extrema that were input
            ELi, EHi, DLi, DHi = getcroplimits(EL, EH, DL, DH, AP)
	ND1 = ND[DLi:DHi+1]    ### now crop the rows
        ND1 = np.swapaxes(ND1, 0, 1)
        ND1 = ND1[ELi:EHi+1] ###crop the columns
        ND1 = np.swapaxes(ND1, 0, 1)
        crops = ELi, EHi, DLi, DHi
	allData.append((ND1, AP, FN, crops))
	data.append(ND1)



mean = np.mean(data, axis = 0)
std = np.std(data, axis = 0)
#stdmean = np.mean(std)
#print stdsum

plt.clf()
plt.title(titleav)
plt.contourf(AP[0][ELi:EHi+1],AP[1][DLi:DHi+1],mean,100,cmap="jet")
plt.xlabel(titlex)
plt.ylabel(titley)
plt.colorbar()
plt.savefig("Average CIU Fingerprint.png")
plt.close()

plt.clf()

##### Change v in order to properly scale the standard deviation color scale ######
v = np.linspace(0,s, 100, endpoint=True)
plt.title(titlestd)
plt.contourf(AP[0][ELi:EHi+1],AP[1][DLi:DHi+1],std,v,cmap="jet",vmin = 0, vmax = 1)
plt.xlabel(titlex)
plt.ylabel(titley)
v = np.linspace(0,s,6,endpoint=True)
plt.colorbar(ticks=v)
plt.savefig("Fingerprint Standard Deviation.png")
plt.close()

a = np.vstack([AP[0][ELi:EHi+1],mean])
b = np.insert(AP[1][DLi:DHi+1],0,np.nan)
b = np.hsplit(b,len(b))
c = np.hstack((b,a))
np.savetxt("Average CIU Fingerprint.csv", c, delimiter=",")

a = np.vstack((AP[0][ELi:EHi+1],std))
c = np.hstack((b,a))
np.savetxt("Fingerprint Standard Deviation.csv", c, delimiter=",")
