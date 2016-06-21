#!/usr/bin/env python

# made for comparing unfiltered and filtered scorefiles for Rosetta enzdes post analysis

import argparse
import collections
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages       


def getDataFromSc(axes, f, uf):
    """initializes two dictionaries and poplulates them based on -f and -u options"""
    fComboDict = collections.defaultdict(list)
    ufComboDict = collections.defaultdict(list)
    maxX = -10000
    maxY = -10000
    minX = 10000
    minY = 10000

    for fileType in [uf, f]:
        for i, item in enumerate(fileType):
            with open(item) as f:
                header = f.readline().split()
                indices = [header.index(a) for a in axes]
                
                for line in f:
                    lineList = line.split()
                    if (not lineList) or (lineList[0].startswith("#")) or (lineList[0][0].isalpha()):\
                       continue
                    try:
                        descStr = lineList[indices[-1]]
                        foundDesc = re.search('A([0-9]+)_P([0-9]+)', descStr).group()
                    except AttributeError:
                        continue
                    
                    pointList = [lineList[i] for i in indices[:-1]]
                    pointTuple = tuple(map(float, pointList))
                    if pointTuple[0] > maxX:
                        maxX = pointTuple[0]
                    if pointTuple[0] < minX:
                        minX = pointTuple[0]
                    if pointTuple[1] > maxY:
                        maxY = pointTuple[1]
                    if pointTuple[1] < minY:
                        minY = pointTuple[1]
                    
                    if fileType == uf:
                        ufComboDict[foundDesc].append(pointTuple)
                    else:
                        fComboDict[foundDesc].append(pointTuple)
    return ufComboDict, fComboDict, minX, maxX, minY, maxY


def genPlots(ufDict, fDict, minX, maxX, minY, maxY, axes, name):
    """makes pdf of plots - one plot for each A[0-9]_P[0-9]"""
    with PdfPages(name) as pdf:
        for entry in ufDict:
            print 'Making plot for ' + entry
            xuf, yuf = zip(*ufDict[entry])
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.scatter(xuf,yuf,c='b',marker='o',label='Unfiltered structures')

            try:
                xf, yf = zip(*fDict[entry])
                ax1.scatter(xf,yf,c='orange',marker='x',label='Filtered structures')
            except:
                pass
            plt.legend(loc='upper right')
            plt.title(entry, fontsize = 30)
            plt.xlim(minX, maxX)
            plt.ylim(minY, maxY)
            plt.xlabel(axes[0], fontsize = 20)
            plt.ylabel(axes[1], fontsize = 20)
            pdf.savefig(fig)
            plt.close()       


def main(xaxis, yaxis, filteredFiles, unFilteredFiles, name):
    """create axes variable and calls previous functions"""
    axes = [xaxis,yaxis,'description']
    ufDict, fDict, minX, maxX, minY, maxY = getDataFromSc(axes,filteredFiles,unFilteredFiles)
    genPlots(ufDict, fDict, minX, maxX, minY, maxY, axes, name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generates scatter plot of data from rosetta score files")
    parser.add_argument("-x", "--xaxis", \
                        help="criterion to be plotted on x-axis (default: total_score)", \
                        default='total_score')
    parser.add_argument("-y", "--yaxis", \
                        help="criterion to be plotted on y-axis (default: SR_1_total_score)", \
                        default='SR_1_total_score')
    parser.add_argument("-n", "--name", default='postProcessPlot.pdf',
                        help='name of output pdf (default: postProcessPlot.pdf')
    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-f", "--filtered", nargs='*', required=True,\
                        help="one or more filtered score files from which data is pulled")
    requiredO.add_argument("-u", "--unfiltered", nargs='*', required=True,\
                        help="one or more unfiltered score files from which data is pulled")
    args = parser.parse_args()

    main(args.xaxis,args.yaxis,args.filtered,args.unfiltered, args.name)

    


