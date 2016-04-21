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

import numpy as np
import scipy.ndimage as ndi
import scipy.signal as signal
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as color
from matplotlib.backends.backend_pdf import PdfPages
import csv

##### Introduction and Options ##########
print("Collision Induced Unfolding Fingerprint Feature Detection and Characterization")
print("all files must end with '_raw.csv'  to be read")
answer = input("Axis titles are in Trap Collision Voltage (V) and Drift Time (ms) by default. Please type y to continue or n to change these settings: ")
lower = answer.lower()
assert lower == 'y' or lower == 'n', "Your input must be y or n"
if lower == 'n':
    titlex = input("Please title your X axis:") #For our lab normally Trap Collision Voltage
    titley = input("Please title your Y axis:") #Normally Drift time or CCS
elif lower == 'y':
    titlex = "Trap Collision Voltage (V)"
    titley = "Drift Time (ms)"

window = input("Please choose a Savitzky-Golay Smoothing Window(must be odd integer greater than 2, we recommend 5 to start) ")
try:
    window = int(window)
except (TypeError, ValueError):
    print("you must enter an odd integer")
    sys.exit()
assert window % 2 != 0, "number must be odd"
threshold = float(input("Please enter the signal intensity threshold(e.g. '80'):  "))
assert 0 < threshold < 100, "You must enter a value between 0 and 100"
threshold = float(threshold)/100
scalef = input("Please enter a scaling factor for your data(e.g. '0.5', '1', or '2'):  ")
try:
    scalef = float(scalef)
except (TypeError, ValueError):
    print("Your input must be a number")
    sys.exit()
    

###### If you want to stop the annoying messages at the the front of this program, comment out everthing under Introduction and Options ######
###### You will then want to uncomment the parameters below and enter your values #######
###### Examples(optimized for APT) #####
#titlex = "Trap Collision Voltage (V)"
#titley = "Drift Time (ms)"
#window = 5
#threshold = 80
#scalef = 2
#####################



#############Definitions############################  



def getCIUData(fname):
	rawData=np.genfromtxt(fname, missing_values=[""], filling_values=[0], delimiter=",")
	return rawData
	
def getAxisParameters(rawData):
        yaxis = rawData[1:,0]
        xaxis = rawData[0,1:]
        return xaxis, yaxis
		
def CIUnormalize(rawData):
	strippedMat=rawData[1:,1:] ### strip the axes off
	maxInt = np.amax(strippedMat, axis = 0) ### find the maximum value for each energy
	norm = strippedMat/maxInt  ### Normalize to max value
	return norm
	
def generateCIUPlot(fname, features):
        OD = getCIUData(fname[1]) ### Read the file
        AP = getAxisParameters(OD)  ### get the axes
        title = fname[1][:-8]### Generate a title based on the filename
        pdf = PdfPages(title + ".pdf")  ### create a pdf
        plt.clf()  ### make a plot
        plt.title(title) ### plot the title
        plt.contourf(AP[0],AP[1],fname[0],100,cmap='jet')  ### plot the data
        plt.xlabel(titlex)    ### label the axes
        plt.ylabel(titley)
        plt.colorbar()  ### plot a colorbar
        pdf.savefig() ### save the figure
        
### Here is the template for making another plot ####
        #plt.clf()
        #plt.title(title+": Detected Features")
        #plt.contourf(AP[0],AP[1],mk,100, cmap = 'jet')
        #plt.colorbar()
        #pdf.savefig() ### save the figure
        #pdf.close()  ### close the file
#######################################################
        
        plt.clf()
        plt.title(title+":1 Detected Features")
        #plt.contourf(AP[0],AP[1],features,3, cmap = color.ListedColormap(['white', 'red', 'orange', 'green', 'blue', 'indigo','violet']))
        features = np.flipud(features)
        plt.imshow(features, aspect = 'auto',cmap = color.ListedColormap(['white', 'red','orange','green','blue', 'indigo', 'violet', 'black'], N = np.max(features)+1), extent = [AP[0][0], AP[0][-1],AP[1][0],AP[1][-1]], interpolation = "spline36")
        plt.tick_params(axis = 'x', which = 'both', bottom = 'off', top = 'off', left = 'off', right = 'off')
        plt.tick_params(axis = 'y', which = 'both', bottom = 'off', top = 'off', left = 'off', right = 'off')
        plt.colorbar()
        #plt.imshow(features, aspect = 'auto',cmap = color.ListedColormap(['white', 'blue','blue','blue','blue','blue','blue']), extent = [AP[0][0], AP[0][-1],AP[1][0],AP[1][-1]], interpolation = "spline36")
        pdf.savefig() ### save the figure
        pdf.close()  ### close the file
        
     
            
### Now do the work ###

### generate a list of files ###    
files=os.listdir(".")
fnames=[x for x in files if x.endswith('_raw.csv')]  ### Make a list of files in the directory ending with _raw.csv


### Extract all the data and store it in an array array ###

allData = [] ## an array to store the data ##
print("Files to be processed")
for item in fnames:
        print(item)
        OD = getCIUData(item)  ### Read the file
        ND = CIUnormalize(OD) ### Normalize
        AP = getAxisParameters(OD) ### The axes
        allData.append((ND, item, AP))  ### put the normalized data, the filename, and the axis parameters into the above array
print("\n")
print("\n")

### Find the features, write them to a csv, also create individual PDFs ###

with open('Features.csv', 'w') as csvfile:  ### Create and open a csv file
    writer = csv.writer(csvfile, delimiter = ",")  ### define the function that can write in the csv file
    writer.writerow(["File Name","# Features","Appearance(CV)","Termination(CV)","Stability(V)",titley])    #### Format the CSV File
    for data in allData:
        oimg = data[0]  ###Data [0] is the "original image" or original matrix
        AP = data[2]    ###Axis parameters
        simg = signal.savgol_filter(oimg, window,2, axis = 1) ### smooth the original image with the window parameter defined above
        #simg = CIUnormalize(simg)
        gimg = np.gradient(simg)  ### calculate the gradient of the smoothed image
        gimgy = gimg[0]  ###we only want the first part of gimag, which is dy(when you view the matrix as a CIU plot) 
        gimgy1 = gimgy < .001 ### now pick out only the areas with near zero gradient (these are the "features")
        gimgy1 = gimgy1 > -.001 
        label_im, nb_labels = ndi.label(gimgy1.astype(int))  ### now label the features in the array
        #label_im = np.array(label_imy) 
        f = ndi.find_objects(label_im)### extract the footprints of each feature
        
        
        #### Go back and refine the features  ###
        
        nf = []  ### an array for refined features
        for i in f: ### for every feature in f
            oslice = simg[i] ### find the detected feature back in the original smoothed image
            nslice = oslice**scalef ### apply a scaling factor to separate or combine features as you see fit
            nslice = oslice > np.max(oslice)*threshold  ### Take only the signal above a certain threshold
            newimage = np.zeros_like(oimg) ### now make a blank array the same shape as the original data
            newimage[i] = nslice ### now past the refined feature onto the blank image 
            ni = newimage.astype(int)  ### Change the feature to 1s or 0s for now
            nf.append(ni) ###add the new feature to the array above
            
        newimage = np.sum(nf, axis = 0) ### the new image is the sum of all the features we added to the array above
        label_im_new, nb_labels_new = ndi.label(newimage) ### label the refined features
        find = ndi.find_objects(label_im_new)  ### find their new footprints
        
        
#### Comment this line if you want to turn the plots off ##############################
        generateCIUPlot(data, label_im_new) ### make the plots for each file
#######################################################################################
        
        writer.writerow([data[1][:-8], nb_labels_new]) ### write a new line for your filename
        i = 1
        print(data[1])
        for s in find: ### For every feature in the plot
                    dts = s[0].indices(len(AP[1])) ### Get the starting and stopping points in drift time space
                    cvs = s[1].indices(len(AP[0])) ### get the start and stop in CV space
                    print("feature ", i)
                    print("Stability = ", AP[0][cvs[1]-1] - AP[0][cvs[0]] , "V from ")#, AP[0][cvs[0]], " to " #, AP[0][cvs[1]-1]##Stability is the total Voltage the feature survives through
                    print(titley, " = ", (AP[1][dts[1]-1] + AP[1][dts[0]])/2)  
                    writer.writerow(["feature " + str(i)," ",AP[0][cvs[0]],AP[0][cvs[1]-1], AP[0][cvs[1]-1] - AP[0][cvs[0]],(AP[1][dts[1]-1] + AP[1][dts[0]])/2])  ### DT is reported as an average
                    i = i+ 1
        print("\n")
        print("\n")
                    
csvfile.close()
