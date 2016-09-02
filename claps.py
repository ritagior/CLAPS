#! /usr/bin/env python


__author__ = "Rita Giordano and Ricardo Leal"
__copyright__ = "Copyright 2013"
__credits__ = ["Rita Giordano", "Ricardo Leal"]
__license__ = ""
__version__ = "0.1"
__maintainer__ = "Rita Giordano"
__email__ = "rgiordano@gmx.com"
__status__ = "Beta"



# Insert parse  to change the file path from command line

from optparse import OptionParser

parser = OptionParser(usage="%prog --XSCALEfile=<LP filename> --outname=<output dendogram>")
parser.add_option("-i","--XSCALEfile", dest="XSCALEfile", help="XSCALE log file")
parser.add_option("-o","--outname", dest="outname", help="output dendogram file name")
(options, args) = parser.parse_args()

if options.XSCALEfile is None:
    parser.print_help()
    parser.error('All options are mandatory.')



import sys
import os
import numpy as np

import rpy2
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as np2r


def parseXscaleLpFile(filename=options.XSCALEfile):
    """
        
        Will parse this bit of XSCALE:
        
        ******************************************************************************
        OVERALL SCALING AND CRYSTAL DISORDER CORRECTION
        ******************************************************************************
        
        CORRELATIONS BETWEEN INPUT DATA SETS AFTER CORRECTIONS
        
        DATA SETS  NUMBER OF COMMON  CORRELATION   RATIO OF COMMON   B-FACTOR
        #i   #j     REFLECTIONS     BETWEEN i,j  INTENSITIES (i/j)  BETWEEN i,j
        
        1    2        3614           0.983            0.9560         0.4878
        
        """
    
    startTableTxt = "CORRELATIONS BETWEEN INPUT DATA SETS AFTER CORRECTIONS"
    dataStartingAfterLines = 5
    dataArr = None
    
    if os.path.exists(options.XSCALEfile) is False:
        print options.XSCALEfile + ' does not exist...'
        sys.exit()
    else :
        linesLeftToParseTableContent = dataStartingAfterLines
        data = []
        f = open(options.XSCALEfile,"r")
        for line in f:
            if line.find(startTableTxt) >= 0 :
                linesLeftToParseTableContent -= 1
            elif  linesLeftToParseTableContent != dataStartingAfterLines and linesLeftToParseTableContent is not 0 :
                linesLeftToParseTableContent -= 1
            elif linesLeftToParseTableContent == 0 :
                dataLine = line.rstrip().split();
                if dataLine == [] :
                    break
                else :
                    data.append(dataLine)
        
        
        f.close()
        dataArr = np.array(data,dtype=(float))
    return dataArr

ccTable = parseXscaleLpFile(options.XSCALEfile)


# This part of the code will read the elements from XSCALE.LP and put in a matrix form.
# This matrix will be read by R to produce the cluster analysis

data = ccTable

Matrix = {}

for line in data:
	if line is not None and len(line) is not 0:
		Matrix[(line[0],line[1])] = line[3]


keys_x = [line[0] for line in Matrix.keys()]
keys_y = [line[1] for line in Matrix.keys()]


dataMatrix = []

for index_row,row in enumerate(range(int(min(keys_x)),  int(max(keys_x)+1))):
    line = []
    for index_col,col in enumerate(range(int(min(keys_y)),  int(max(keys_y)+1))):
	
        # First columns of the matrix must be 1 followed by 0s
        if index_row == 0 and index_col == 0 :
            line.append(1)
        elif index_row > 0 and index_col == 0 :
            line.append(0)
		
		# Print file contents
        if row == col :
            line.append(1)
        elif Matrix.has_key((row,col)):     # check for key before fetch
            line.append(Matrix[(row,col)])
        else:
            line.append(0)
                
    dataMatrix.append(line)



line = []

# Add a new line at the end of the matrix with 0s and 1 at the end of the line
for col in range(int(min(keys_y)),  int(max(keys_y)+1)):
	if col < int(max(keys_y)+1) :
            line.append(0)
line.append(1)
dataMatrix.append(line)
dataMatrixArr = np.array(dataMatrix, dtype=(float))

print len(dataMatrix)


r = robjects.r

    
# Open the data
data = np2r.numpy2ri(dataMatrixArr)


data_m = robjects.r['as.matrix'](data)


    
# transpose the data's matrix
cc=robjects.r['t'](data)

    
cm=robjects.r['as.matrix'](cc)


## Convert thr robjects in numpy objects
cm_np=np.array(cm)
    
# define distance function
d = np.sqrt(1-cm_np**2)
    
## convert from numpy to rpy2
dist = np2r.numpy2ri(d)
    
dm=robjects.r['as.matrix'](dist)
    
## Define that dit is a distance matrix
dist_m=robjects.r['as.dist'](dm)
    # Perform the cluster analysis and plot the dendrogram

r.X11()
p=r.plot(robjects.r['hclust'](dist_m, "average"), sub="", xlab="", ylab="dist(i,j)")
raw_input()

r.pdf(options.outname)
p=r.plot(robjects.r['hclust'](dist_m, "average"), sub="", xlab="", ylab="dist(i,j)")
r.dev_off()


"""
		print "%2.3linef  " % 0 ,
print "%2.3f  " % 1 ,
print
print
"""

