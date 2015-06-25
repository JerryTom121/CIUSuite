import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import csv
import scipy.stats as stats
### Functions ####
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
        return norm

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

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin() ###get the index of the value nearest to the input value
    return idx
    
###############################

####### Introduction and options ################
print "CIU Fingerprint Comparison"
print "This tool will generate difference plots between CIU fingerprints. There are three modes: Basic, Batch, and Cluster."
print "Basic mode requires two files as inputs, and will generate a plot of File1 - File2" + '\n' + "Batch mode will generate a plot of Ref - i, for all files i in the directory."+"\n"+"Cluster mode will calculate the RMSD for each pair of files, and use k-medoids clustering to sort them"
answer = raw_input("Please enter basic, batch, or cluster: ")
lower = answer.lower()
assert lower == 'basic' or lower == 'batch' or lower == 'cluster', "You must enter basic, batch, or cluster"

def optionsmenu():
    print "Settings:"
    print "1. Cropping"
    print "2. Intensity Scaling"
    print "3. Axis Titles"
    print "4. Exit Menu"
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


#### Default Params ####
titleyn = 'n'
cropyn = 'n'
titlex = "Trap Collision Energy (V)"
titley = "Drift Time (ms)" 
s = 1

#### Menu ####       

print "Default settings are no cropping, no scaling of intensities, and axes titles of Trap Collision Voltage(V) and Drift Time(ms)"
optionsyn = raw_input("Please press y to continue or n to change these settings  ")
optionsyn = optionsyn.lower()
if optionsyn == 'n':
    while True:
        option = optionsmenu()
        if option == "1":
            cropyn = 'y'
            EL, EH, DL, DH = cropinputs() 
            print "Crop Inputs Saved, Please choose another option, or exit the menu"
            continue
          
        if option == "2":
            print "This is the intensity scale for the difference plot. Lower values will exaggerate differences, higher ones will make them appear more subtle" 
            s = raw_input("Please enter a value from 0 to 1: ")
            s = float(s)
            continue
        if option == "3":
            title = raw_input("Plot title: ")
            titleyn = "y"
            titlex = raw_input("X axis title: ")
            titley = raw_input("Y axis title: ")
            continue
        if option == "4":
            break
        else:
            continue
            

################################################
################################################

if lower == "Basic" or lower == "basic":
    print "Basic Mode:"
    
    print "Please Select File1 and File2"
    files=os.listdir(".")
    files=[x for x in files if x.endswith('_raw.csv')]
    i = 1
    for item in files:
        print i, item
        i = i +1
        
    fileid1 = raw_input("File1 Number : ")
    file1 = files[int(fileid1) -1]
    fileid2 = raw_input("File2 Number : ")
    file2 = files[int(fileid2) -1]

    ### Do the work ###    
    file1 = getCIUData(file1)
    file2 = getCIUData(file2) ## load the data
    file1n = CIUnormalize(file1[0],file1[1])
    file2n = CIUnormalize(file2[0],file2[1]) ## normalize the data
    AP = getAxisParameters(file1[0]) ## look at the axes
    if cropyn == 'n': ### If no cropping, use the whole matrix
        ELi = 0
        EHi = len(AP[0])
        DLi = 0
        DHi = len(AP[1])
    if cropyn == 'y':  ### If cropping, find the indices that correspond to the extrema that were input
        ELi, EHi, DLi, DHi = getcroplimits(EL, EH, DL, DH, AP)
    file1n = file1n[DLi:DHi+1]    ### now crop the rows
    file1n = np.swapaxes(file1n, 0, 1)
    file1n = file1n[ELi:EHi+1] ###crop the columns
    file1n = np.swapaxes(file1n, 0, 1)
    
    file2n = file2n[DLi:DHi+1]    
    file2n = np.swapaxes(file2n, 0, 1)
    file2n = file2n[ELi:EHi+1]
    file2n = np.swapaxes(file2n, 0, 1)
    difmat = file1n - file2n  ### the difference matrix
    difmat2 = difmat**2
    RMSD = (np.average(difmat2))**0.5
    pRMSD = RMSD*100
    rtext = "RMSD = " +  '%2.2f' % pRMSD
    AP = getAxisParameters(file1[0])
    FN1 = str(file1[1])
    FN2 = str(file2[1])
    if titleyn != 'y':
        title = FN1[:-8] + "-" + FN2[:-8]
    plt.clf()  ### Make the plot
    plt.title(title)
    v = np.linspace(-s,s, 100, endpoint=True)
    plt.contourf(AP[0][ELi:EHi+1],AP[1][DLi:DHi+1],difmat,v,cmap='jet', ticks = "none")
    plt.tick_params(axis = 'x', which = 'both', bottom = 'off', top = 'off', left = 'off', right = 'off')
    plt.tick_params(axis = 'y', which = 'both', bottom = 'off', top = 'off', left = 'off', right = 'off')
    plt.annotate(rtext, xy=(260,10), xycoords='axes points')
    plt.xlabel(titlex)    
    plt.ylabel(titley)

    v = np.linspace(-s,s,6,endpoint=True)
    plt.colorbar(ticks=v)

    plt.savefig(title +'.png')
    plt.close()
    difmat2 = difmat**2
    RMSD = (np.average(difmat2))**0.5
    pRMSD = RMSD*100
    print "\n" + "Processing ", file1[1] +" - " + file2[1]
    print "%RMSD = " , '%2.2f' % pRMSD
    

######################################################################
######################################################################

if lower == "Batch" or lower == "batch":
    #subroutine 2
    print "Batch Mode:" + '\n' + "all files must end with '_raw.csv' to be read with Batch Mode"
    files=os.listdir(".")
    files=[x for x in files if x.endswith('_raw.csv')]
    print "Please select a reference file"
    i = 1
    for item in files:
        print i, item
        i = i +1
        
    fileid = raw_input("File Number : ")
    ref = files[int(fileid) -1]
    del files[int(fileid) - 1]

    ref = getCIUData(ref)
    refn1 =CIUnormalize(ref[0],ref[1])  ##normalized reference file
    AP = getAxisParameters(ref[0])
    with open(ref[1].rstrip("_raw.csv") + " RMSDs.csv", "wb") as csvfile:
        writer = csv.writer(csvfile, delimiter = ",")
        writer.writerow(["File Name","RMSD"]) 
        writer.writerow([""]) 
        for item in files:
            i = getCIUData(item)
            t = i[1] ## the name of the file for later
            inorm = CIUnormalize(i[0],i[1])  #normalized i file
            AP = getAxisParameters(ref[0])
            if cropyn == 'n':
                ELi = 0
                EHi = len(AP[0])
                DLi = 0
                DHi = len(AP[1])
            if cropyn == 'y':
                ELi, EHi, DLi, DHi = getcroplimits(EL, EH, DL, DH, AP)
                    
            refn = refn1[DLi:DHi+1]  ### do the cropping here
            refn = np.swapaxes(refn, 0, 1)
            refn = refn[ELi:EHi+1]
            refn = np.swapaxes(refn, 0, 1)
            inorm = inorm[DLi:DHi+1]  
            inorm = np.swapaxes(inorm, 0, 1)
            inorm = inorm[ELi:EHi+1]
            inorm = np.swapaxes(inorm, 0, 1)
            
            difmat = refn - inorm
            difmat2 = difmat**2
            RMSD = (np.average(difmat2))**0.5
            pRMSD = RMSD*100
            plt.clf()
            title1 = ref[1][:-8]
            title2 = t[:-8]
            title = title1 + " - " + title2
            plt.title(title)
            v = np.linspace(-s,s, 100, endpoint=True)
            plt.contourf(AP[0][ELi:EHi+1],AP[1][DLi:DHi+1],difmat,v,cmap='jet', ticks = "none")
            rtext = "RMSD = " +  '%2.2f' % pRMSD
            plt.annotate(rtext, xy=(260,10), xycoords='axes points')
            plt.xlabel(titlex)    
            plt.ylabel(titley)
            v = np.linspace(-s,s,6,endpoint=True)
            plt.colorbar(ticks=v)
            plt.savefig(title +'.png')
            plt.close()

            writer.writerow([i[1],pRMSD]) 
            
            print "\n" + "Processing ", ref[1].rstrip("_raw") +" - " + i[1]
            print "%RMSD = " ,  "%.2f" % pRMSD
        print "\n" + "A list of these values has been output to a CSV file for further analysis"
 
        
#############################################################
if lower == "Cluster" or lower == "cluster":
    from Pycluster import  kmedoids
    files=os.listdir(".")
    files=[x for x in files if x.endswith('_raw.csv')]
    c = int(raw_input("Please enter the number of clusters to search for: "))
    rmsdmat = []
    for item in files: 
        item = getCIUData(item) 
        inorm = CIUnormalize(item[0],item[1])
        rmsdmat1 = []
        for item1 in files:
            item1 = getCIUData(item1)
            i1norm = CIUnormalize(item1[0],item1[1])
            AP = getAxisParameters(item1[0])
            
            if cropyn == 'n':
                ELi = 0
                EHi = len(AP[0])
                DLi = 0
                DHi = len(AP[1])
            if cropyn == 'y':
                ELi, EHi, DLi, DHi = getcroplimits(EL, EH, DL, DH, AP)
                   
            i1norm1 = i1norm[DLi:DHi+1]    
            i1norm1 = np.swapaxes(i1norm1, 0, 1)
            i1norm1 = i1norm1[ELi:EHi+1]
            i1norm1 = np.swapaxes(i1norm1, 0, 1)
            
            inorm1 = inorm[DLi:DHi+1]   
            inorm1 = np.swapaxes(inorm1, 0, 1)
            inorm1 = inorm1[ELi:EHi+1]
            inorm1 = np.swapaxes(inorm1, 0, 1)
            
            difmat = inorm1[1] - i1norm1[1]     
            difmat2 = difmat**2
            RMSD = (np.average(difmat2))**0.5
            pRMSD = RMSD*100
            rmsdmat1.append(pRMSD)
        rmsdmat.append(rmsdmat1)
    rmsdmat = np.array(rmsdmat)    
    rmsdmat = 1 - rmsdmat
    np.savetxt("DifferenceMatrix.csv", rmsdmat, delimiter = ",")
    clusterid, error, nfound = kmedoids (rmsdmat, nclusters=c,npass = 10)
    print "Error = ", error
    print "Found this configuration " , nfound, " out of 10 times"
    with open("Clusters.csv", "wb") as csvfile:
        writer = csv.writer(csvfile, delimiter = ",")
        writer.writerow(["item #","File Name","Cluster #"]) 
        writer.writerow([""]) 
        for i in range(0, len(clusterid)):
            print i, files[i], " cluster = ", clusterid[i]
            writer.writerow([i,files[i], clusterid[i]])
        writer.writerow([""])
        writer.writerow(["Note : The cluster number is defined as the item number of the centroid of the cluster."])    
            
    print "\n" + "The cluster number is defined as the item number of the centroid of the cluster."
    print "\n" + "A list of the cluster assignments has been output to Clusters.csv."
    print "\n" + "If you would like to use a different clustering algorithm, a copy of the distance matrix has been output to DifferenceMatrix.csv"
        

else:
    sys.exit()
