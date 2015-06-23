CIUSuite is a suite of tools for the analysis of CIU data. These tools are in the form of python scripts with unique functionality.

CIUSuite was authored by Joseph Eschweiler and Jessica Rabuck and is awaiting publication, Please do not distribute these programs without prior consent

Correspondence to Joe at joeesch@umich.edu

*******

All scripts written in python2.7, and some have different dependencies. If you're savy with python, you may have all of the dependent libraries installed already. If you are not savy, We recommend using The Enthough Canopy python environment to deal with dependencies. https://store.enthought.com/downloads/#default  The free version should do, and you can also upgrade with an academic license. 

********


###CIUSuite_plot.py###

Dependencies : numpy, matplotlib, scipy

The foundation of the CIUSuite is CIUSuite_plot.py. This script will read all files in the working directory
ending in "_raw.csv" and generate CIU fingerprints for them. When you run this script, you'll be asked
if you want to use default or custom options for the axis titles, the titles will always be the filename
stripped of its "_raw.csv" component. The outpout of CIUSuite_plot is a file with the same name as the input file,
stripped of "_raw.csv" and replaced with ".png"

Additionally, CIUSuite_plot normalizes and smoothes data. By default, data is smoothed using a Savitsky-Golay filter with a polynomial order of 2 and window length of three(very very minor smoothing). You can increase the smoothing by increasing the window length, using odd integers.

The format for the input _raw.csv file must be as follows:

(blank)   trap1    trap2    trap3   etc

Dt1       Data     Data     Data

Dt2       Data     Data     Data

Dt3       Data     Data     Data

etc


where trap values are the trap collision energies, and Dt values are the drift times. 

The input file format for all CIUgen tools will be the same.

***note*** 
The dimensions of the input files do not need to be the same for CIUgen to work. In all subsequent scripts, 
uniform dimensionality is required.
**********




###CIUSuite_stats.py###

Dependencies : numpy, matplotlib

CIUSuite_stats is a simple program to generate average and standard deviation plots from CIU data.
CIUSuite_stats will calculate the average and standard deviation matrices for all _raw.csv files in
the working directory. 

To run CIUSuite_stats, first make sure all the files are of the same dimensionality. Upon executing the
script, you'll be asked to keep or change the default title and axis title parameters. Next, you'll be
asked set the scale of the standard deviation colorbar. To accurately reflect the deviation relative to the
data, the best option is 1, however we have found that to produce adequate contrast a value of 0.6 is desired.

The output of CIUSuite_stats are two files, an average and standard deviation plot in .png format, and normalized
.csv files for both the average and standard deviation. The normalized .csv files will be useful for the following
scripts.  



###CIUSuite_compare.py###

Dependencies : numpy, matplotlib, scipy, csv, PyCluster(for cluster mode only)

This script is useful for comparison of CIU fingerprints. It operates in three modes, Basic, Batch, and Cluster.
Upon executing the script, you'll be asked to choose which mode to run in. 

BASIC MODE:
basic mode requires input of two files, file1 and file2. the script will than calculate the difference
matrix file1 - file2 and present it in a CIU-style contour plot. The contour plot has similar options to those
above, allowing one to change the axis labels and the plot title. 

The difference plots also require a colorbar scaling factor, set as 1 by default. Increase this value if you find that
your data is cut off, and decrease it if you do not see enough contrast in your data.

Additionally, basic mode will calculate the root mean square deviation (RMSD) between the two datasets. The RMSD
is reported as RMSD*100% for ease of use. The RMSD gives a the relative deviation of file2 from file1. The usefulness of
RMSD may be low as an absolute measure, but can be useful when comparing many datasets. 



BATCH MODE:
Batch mode performs the same basic functionality as basic mode, however applies it on a larger scale. Batch mode requires
one input file, termed the Reference File. This file is what all other files in the directory will be compared to. Formally,
batch mode calculates the distance matrix and RMSD for REF - i for all i_raw.csv in the directory. 

The output of batch mode will be a difference plot for each file i, and the RMSD between REF and I.
In Addition to the command prompt outputs, A list of files and their RMSDs will be output to a CSV.


CLUSTER MODE:
Calculates RMSDs for each pair of fingerprints in the directory. Then uses K-Medioids clustering to find optimal cluster
configurations. Takes as input the expected number of clusters. If you don't know the number of clusters, see the "ERROR"
output on the command line and try to minimize it by guessing different number of clusters within a reasonable range. 
The output of this mode is a csv file listing the clusters each file belongs to, and another csv file containing the distance
matrix if you want to do some other form of clustering(for example, hierarchical). 




###CIUSuite_Analysis.py###

Dependencies : Numpy, matplotlib

CIUSuite_analysis provides statistical analysis useful in narrowing the search window for high throughput screens and allowing the user to choose high information content regions for classifying spectra. CIU_analysis takes as input raw or normalized files
ending in _typeX_raw.csv. The type tag is necessary for grouping the "training" data in the beginning of the program. 

Once the training data is grouped and averaged, each piece of data will then be compared back to its parent group and to the opposing group via scaled deviation scores across each energy.

(see upcoming paper for further explaination of the score) 

By calculating an average SDS for each group compared to its own average and compared to the opposing group average
a more useful score is produced. These  values are plotted in two figures that are the final output of the program,  X2 vs CE
against type I and against type II. These figure show, as a function of Collision Voltage, where the greatest differences in the
type I and type II SDS occur. The errorbars on these figures represent the 2 standard deviations. The error bars are EXTREMELY
important, especially when working with small datasets. 




###CIUSuite_Score.py ###

Dependencies : numpy, matplotlib, csv

This program shares much functionality with CIUSuite_Analysis, however in addition to training data, it also accepts unclassified data,
tagged with _UK_Raw.csv. The ultimate function of the program is to provide scores corresponding to the TYPEI or TYPEII character of 
a spectrum. 

We have found that rather than using the entire CIU spectrum for SDS scoring, focusing on specific collision energies that show the
greatest difference in SDS2 score (AS SHOWN BY CIUSuite_ANALYSIS) yeilds more robust scores. Thus, CIUSuite_score provies options to crop
the data in both drift time and CV space. The drift time cropping function will search for the drift time input, rounded to 1 decimal place

example: If the drif time value in the RAW data is 10.122, the command line input should be 10.1. 

The collision voltage input is an integer. 


After the data has been cropped, Each spectrum is grouped according to its tag, and the normalized average matrix is calculated for each of the types.
Using these matrices as expected values, a type I and type II score will be calculated for each spectrum, including the training data and uknowns. 

The typeI and typeII scores are calculated via two methods, both output to the command window and to Score.csv.

METHOD 1

SDS is calculated for each collision energy and the final score is the sum of SDS accross all energies for typeI and typeII. 



###CIUSuite_detect###

Dependencies : Numpy, scipy, matplotlib, csv

CIUSuite_detect is an elementary feature detection algorithm for CIU fingerprints. We have found CIUSuite works well for detecting and characterizing features
that are well resolved in DT or CV, however the accuracy of the feature detection is poor when dealing with poorly resolved features. CIU_detect works by using the first derivative test to identify local maxima in the data. Then, using user-defined scaling factors(exponents, ie 2 for square the data, 0.5 square root the data) the features can be refined to separate or combine features(not an exact science,see Zhong et al, DOI: 10.1002/anie.201403784 for more info)

Once features are refined and validated they are characterized by centroid drift time, centroid energy, and energy range. 

Inputs:
	Savgol window : Smoothing the data is very important for accurate detection of features. Higher windows give better results, a good start is 11 or 13. 

	Scaling Factor: 2 to square the data, helps to separate some features in conjunction with threshold.

	Intensity Threshold: %signal intensity threshold to take as a feature. ~80 usually gives accurate representation and separation of features. 

The script is handy for batch processing, it takes all files _raw.csv.

Improvements to this script are planned








