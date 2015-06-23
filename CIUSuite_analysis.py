import os
import numpy as np
import matplotlib.pyplot as plt

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
	x = rawData[0][1:]
	y = rawData[:,0]
	normOut = np.vstack((x,norm))
	normOut = np.column_stack((y,normOut)) ### Add the axes back to the array
	normOut[0][0] = np.nan
        return normOut, norm  ### Output an array with and without axes

################################################################################
	
### Batch data into allData, and respective typeData. Also store axis parameters for reconstructing graphs later ###
files=os.listdir(".")
files=[x for x in files if x.endswith('_raw.csv')]  #get the files

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
for item in files:
	OD = getCIUData(item)
	Mat = OD[0]
	FN = OD[1]
	ND = CIUnormalize(Mat, FN)
	ND1= np.array(ND[1])
	AP = getAxisParameters(Mat)
	if FN.endswith('_typeI_raw.csv'):### Put type I into the type I bins
		FN = FN[:-8]
		typeIData.append(ND1)
		typeITitles.append(FN)
		axisParametersI.append(AP)
	elif FN.endswith('_typeII_raw.csv'): ### Put type II into the type II bins
		FN = FN[:-8]
		typeIIData.append(ND1)
		typeIITitles.append(FN)
		axisParametersII.append(AP)
	elif FN.endswith('_UK_raw.csv'): ### Put type II into the type II bins
	        print "UK FOUND"
		FN = FN[:-8]
		UKData.append(ND1)
		UKTitles.append(FN)
		axisParametersUK.append(AP)



### More functions that rely on the data above ####
### Calculate the average and standard dev. ###
def groupwise_statistics(group):
        mean = np.mean(group, axis=0)
        std = np.std(group, axis=0)
        return mean,std
	
	
### Plot the Chi^2 of a drug compared to some average spectrum	
def chisq2d(data, group, groupname, title):
	AV = groupwise_statistics(group)
	AV = AV[0]
	AV[AV==0] = np.nan
	difmat = data - AV
	difmat2 = difmat**2
	v = difmat2/AV
	findnan = np.isnan(v)
	v[findnan]=0
	findinf = np.isinf(v)
	v[findinf]= 0
	v=np.array(v)
	np.savetxt(title + " v.csv", v, delimiter = ",")
	x = axisParametersI[0][0]
	colsum = np.sum(v, axis = 0)
	colsum = np.array(colsum)
	plt.clf()
	plt.title(str(groupname)+" Groupwise X2")
        plt.xlabel("Collision Energy(V)")    
        plt.ylabel("X^2")
	plt.plot(x, colsum)
	plt.savefig(groupname + " Groupwise X2.png")
	plt.close()
	return colsum
i =0
dsItoI = []
dsItoII = []## we output the Chi2 scores to these lists
#### Calculate the Chi^2 matrices for both bins

typeIav, typeIstd = groupwise_statistics(typeIData)
typeIIav, typeIIstd = groupwise_statistics(typeIIData)
z1to1 = []
z1to2 = []

for data in typeIData:
        print typeITitles[i]
	#I = chisq2d(data, typeIData, typeITitles[i]+" vs typeI Average", typeITitles[i])
	#II = chisq2d(data, typeIIData, typeITitles[i]+" vs typeII Average", typeITitles[i])
	np.seterr(divide='ignore', invalid='ignore')
	zscore1 = (data - typeIav)/typeIstd*typeIav
        zscore1[np.isnan(zscore1)] = 0
	zscore2 = (data - typeIIav)/typeIIstd*typeIIav
	zscore2 = np.abs(zscore2)
	zscore2[np.isnan(zscore2)] = 0
	zscore1 = np.average(zscore1, axis = 0)
	zscore1 = abs(zscore1)
	zscore2 = np.average(zscore2, axis = 0)
	zscore2 = abs(zscore2)
	z1to1.append(zscore1)
	z1to2.append(zscore2)
	#dsItoI.append(I)
	#dsItoII.append(II)
	i = i+1
	
print np.shape(z1to2)
dsIItoI = [] ## output scores to these lists
dsIItoII = []

z2to1 = []
z2to2 = []
i = 0
for data in typeIIData:
        print typeIITitles[i]
    	#I = chisq2d(data, typeIData, typeIITitles[i]+" vs typeI Average", typeIITitles[i])
	#II = chisq2d(data, typeIIData, typeIITitles[i]+" vs typeII Average", typeIITitles[i])
	zscore1 = (data - typeIav)/typeIstd*typeIav
	zscore1[np.isnan(zscore1)] = 0
	zscore2 = (data - typeIIav)/typeIIstd*typeIIav
	zscore2 = np.abs(zscore2)
	zscore2[np.isnan(zscore2)] = 0
	zscore1 = np.average(zscore1, axis = 0)
	zscore1 = abs(zscore1)
	zscore2 = np.average(zscore2, axis = 0)
	zscore2 = abs(zscore2)
	z2to1.append(zscore1)
	z2to2.append(zscore2)
	#dsIItoI.append(I)
	#dsIItoII.append(II)
	i = i+1
zUK = []
i = 0	
for data in UKData:
        print UKTitles[i]
    	#I = chisq2d(data, typeIData, typeIITitles[i]+" vs typeI Average", typeIITitles[i])
	#II = chisq2d(data, typeIIData, typeIITitles[i]+" vs typeII Average", typeIITitles[i])
	#zscore1 = (data - typeIav)/typeIstd*typeIav
	#zscore1[np.isnan(zscore1)] = 0
	zscore2 = (data - typeIIav)/typeIIstd*typeIIav
	zscore2[np.isnan(zscore2)] = 0
	#zscore1 = np.average(zscore1, axis = 0)
	zscore2 = np.average(zscore2, axis = 0)
	zscore2 = abs(zscore2)
	zUK.append(zscore2)
	#z2to1.append(zscore1)
	#z2to2.append(zscore2)
	#dsIItoI.append(I)
	#dsIItoII.append(II)
	i = i+1
	
	
	
print np.shape(z2to1)
print np.shape(z2to2)	
### Plot X^2 against type I Average for each CE 
avg1 = np.average(z1to1, axis=0)
avg2 = np.average(z2to1, axis = 0)
error1 = np.std(z1to1, axis = 0)*2
error2 = np.std(z2to1, axis = 0)*2
x = axisParametersI[0][0]

plt.clf()
plt.title("Zscore Vs CE against type I")
for item in z1to1:
    a = plt.scatter(x, item, c="blue", s= 50)
for item in z2to1:
    b = plt.scatter(x, item, c="red", s= 50)
plt.scatter(x, avg1, c = "blue", s = 100)
plt.scatter(x, avg2, c = "red", s = 100)
colors = ["orange", "black", "purple"]
i = 0
for item in zUK:
    plt.scatter(x, item, c= colors[i], s = 30)
    i = i+1
plt.errorbar(x, avg1, error1, c="blue",linestyle='None')
plt.errorbar(x, avg2, error2, c="red", linestyle='None')
plt.xlabel("Collision Voltage(V)")    
plt.ylabel("Scaled Z-Score")
plt.legend([a, b],("TypeI Value", "TypeII Zscore"),scatterpoints=1,loc='best',ncol=1,fontsize=15)
plt.savefig("Zscore Vs CE against type I.pdf")
plt.close()


### Plot X^2 Against type II Average for each CE
avg1 = np.average(z1to2, axis = 0)
avg2 = np.average(z2to2, axis = 0)
#avg1 = np.average(dsItoII, axis = 0)
#avg2 = np.average(dsIItoII, axis = 0)
print z2to2
print avg2
error1 = np.std(z1to2, axis = 0)*2
error2 = np.std(z2to2, axis = 0)*2
plt.clf()
plt.title("Scaled Type II Deviation Score Vs Collision Voltage")
for item in z1to2:
    plt.scatter(x, item, c="blue",  s= 25)
for item in z2to2:
    plt.scatter(x, item, c="red",  s= 25)
i = 0
#colors = ["orange", "black", "purple"]
#for item in zUK:
#    plt.scatter(x, item, c= colors[i], s = 30)
#    i = i+1
plt.scatter(x, avg1, c = "blue", s = 100)
plt.scatter(x, avg2, c = "red", s = 100)
plt.errorbar(x,avg1, error1, c="blue",linestyle='None')
plt.errorbar(x, avg2, error2, c="red", linestyle='None')
plt.xlabel("Collision Energy(V)")    
plt.ylabel("Scaled Deviation Score")
plt.legend([a, b],("TypeI Average", "TypeII Average"),scatterpoints=1,loc='best',ncol=1,fontsize=15)
plt.savefig("Zscore Vs CE against type II.pdf")
plt.close()