#!/usr/bin/env python


### source: https://github.com/nf-core/chipseq/blob/master/bin/macs2_merged_expand.py
### own changes and adjustments for snakemake-workflow chipseq are marked with "# AVI: " in comment


# MIT License
#
# Copyright (c) Philip Ewels
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


#######################################################################
#######################################################################
## Created on June 29th 2018 to annotate merged peaks
#######################################################################
#######################################################################

import os
import errno
import argparse
import sys  # AVI: added to create log files
import re

sys.stderr = open(snakemake.log[0], "w")  # AVI


############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

# Description = 'Add sample boolean files and aggregate columns from merged MACS narrow or broad peak file.'
# Epilog = """Example usage: python macs2_merged_expand.py <MERGED_INTERVAL_FILE> <SAMPLE_NAME_LIST> <OUTFILE> --is_narrow_peak --min_replicates 1"""
#
# argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)
#
# ## REQUIRED PARAMETERS
# argParser.add_argument('MERGED_INTERVAL_FILE', help="Merged MACS2 interval file created using linux sort and mergeBed.")
# argParser.add_argument('SAMPLE_NAME_LIST', help="Comma-separated list of sample names as named in individual MACS2 broadPeak/narrowPeak output file e.g. SAMPLE_R1 for SAMPLE_R1_peak_1.")
# argParser.add_argument('OUTFILE', help="Full path to output directory.")
#
# ## OPTIONAL PARAMETERS
# argParser.add_argument('-in', '--is_narrow_peak', dest="IS_NARROW_PEAK", help="Whether merged interval file was generated from narrow or broad peak files (default: False).",action='store_true')
# argParser.add_argument('-mr', '--min_replicates', type=int, dest="MIN_REPLICATES", default=1, help="Minumum number of replicates per sample required to contribute to merged peak (default: 1).")
# args = argParser.parse_args()

############################################
############################################
## HELPER FUNCTIONS
############################################
############################################

def makedir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


############################################
############################################
## MAIN FUNCTION
############################################
############################################

## MergedIntervalTxtFile is file created using commands below:
## 1) broadPeak
## sort -k1,1 -k2,2n <MACS_BROADPEAK_FILES_LIST> | mergeBed -c 2,3,4,5,6,7,8,9 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse > merged_peaks.txt
## 2) narrowPeak
## sort -k1,1 -k2,2n <MACS_NARROWPEAK_FILE_LIST> | mergeBed -c 2,3,4,5,6,7,8,9,10 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse > merged_peaks.txt

def macs2_merged_expand(MergedIntervalTxtFile, SampleNameList, OutFile, isNarrow=False, minReplicates=1):
    makedir(os.path.dirname(OutFile))

    combFreqDict = {}
    totalOutIntervals = 0
    SampleNameList = sorted(SampleNameList)
    fin = open(MergedIntervalTxtFile, 'r')
    fout = open(OutFile, 'w')
    oFields = ['chr', 'start', 'end', 'interval_id', 'num_peaks', 'num_samples'] + [x + '.bool' for x in
                                                                                    SampleNameList] + [x + '.fc' for x
                                                                                                       in
                                                                                                       SampleNameList] + [
                  x + '.qval' for x in SampleNameList] + [x + '.pval' for x in SampleNameList] + [x + '.start' for x in
                                                                                                  SampleNameList] + [
                  x + '.end' for x in SampleNameList]
    if isNarrow:
        oFields += [x + '.summit' for x in SampleNameList]
    fout.write('\t'.join(oFields) + '\n')
    while True:
        line = fin.readline()
        if line:
            lspl = line.strip().split('\t')

            chromID = lspl[0];
            mstart = int(lspl[1]);
            mend = int(lspl[2]);
            starts = [int(x) for x in lspl[3].split(',')];
            ends = [int(x) for x in lspl[4].split(',')]
            names = lspl[5].split(',');
            fcs = [float(x) for x in lspl[8].split(',')]
            pvals = [float(x) for x in lspl[9].split(',')];
            qvals = [float(x) for x in lspl[10].split(',')]
            summits = []
            if isNarrow:
                summits = [int(x) for x in lspl[11].split(',')]

            ## GROUP SAMPLES BY REMOVING TRAILING *_R*
            groupDict = {}
            for sID in [re.split(".narrow_peak|.broad_peak",x)[0] for x in names]:
                gID = re.split("_R.",sID)[0]
                if gID not in groupDict:
                    groupDict[gID] = []
                if sID not in groupDict[gID]:
                    groupDict[gID].append(sID)

            ## GET SAMPLES THAT PASS REPLICATE THRESHOLD
            passRepThreshList = []
            for gID, sIDs in groupDict.items():
                if len(sIDs) >= minReplicates:
                    passRepThreshList += sIDs

            ## GET VALUES FROM INDIVIDUAL PEAK SETS
            fcDict = {};
            qvalDict = {};
            pvalDict = {};
            startDict = {};
            endDict = {};
            summitDict = {}
            for idx in range(len(names)):
                sample = re.split(".narrow_peak|.broad_peak",names[idx])[0]
                if sample in passRepThreshList:
                    if sample not in fcDict:
                        fcDict[sample] = []
                    fcDict[sample].append(str(fcs[idx]))
                    if sample not in qvalDict:
                        qvalDict[sample] = []
                    qvalDict[sample].append(str(qvals[idx]))
                    if sample not in pvalDict:
                        pvalDict[sample] = []
                    pvalDict[sample].append(str(pvals[idx]))
                    if sample not in startDict:
                        startDict[sample] = []
                    startDict[sample].append(str(starts[idx]))
                    if sample not in endDict:
                        endDict[sample] = []
                    endDict[sample].append(str(ends[idx]))
                    if isNarrow:
                        if sample not in summitDict:
                            summitDict[sample] = []
                        summitDict[sample].append(str(summits[idx]))

            samples = sorted(fcDict.keys())
            if samples != []:
                numSamples = len(samples)
                boolList = ['TRUE' if x in samples else 'FALSE' for x in SampleNameList]
                fcList = [';'.join(fcDict[x]) if x in samples else 'NA' for x in SampleNameList]
                qvalList = [';'.join(qvalDict[x]) if x in samples else 'NA' for x in SampleNameList]
                pvalList = [';'.join(pvalDict[x]) if x in samples else 'NA' for x in SampleNameList]
                startList = [';'.join(startDict[x]) if x in samples else 'NA' for x in SampleNameList]
                endList = [';'.join(endDict[x]) if x in samples else 'NA' for x in SampleNameList]
                oList = [str(x) for x in [chromID, mstart, mend, 'Interval_' + str(totalOutIntervals + 1), len(names),
                                          numSamples] + boolList + fcList + qvalList + pvalList + startList + endList]
                if isNarrow:
                    oList += [';'.join(summitDict[x]) if x in samples else 'NA' for x in SampleNameList]
                fout.write('\t'.join(oList) + '\n')

                tsamples = tuple(sorted(samples))
                if tsamples not in combFreqDict:
                    combFreqDict[tsamples] = 0
                combFreqDict[tsamples] += 1
                totalOutIntervals += 1

        else:
            fin.close()
            fout.close()
            break

    ## WRITE FILE FOR INTERVAL INTERSECT ACROSS SAMPLES.
    ## COMPATIBLE WITH UPSETR PACKAGE.
    fout = open(OutFile[:-4] + '.intersect.txt', 'w')
    combFreqItems = sorted([(combFreqDict[x], x) for x in combFreqDict.keys()], reverse=True)
    for k, v in combFreqItems:
        fout.write('%s\t%s\n' % ('&'.join(v), k))
    fout.close()


############################################
############################################
## RUN FUNCTION
############################################
############################################

# AVI: arguments adapted to snakemake
macs2_merged_expand(MergedIntervalTxtFile=snakemake.input[0],
                    SampleNameList=list(snakemake.params.get("sample_control_peak")), OutFile=snakemake.output.get("bool_txt"),
                    isNarrow=snakemake.params.get("narrow_param"),
                    minReplicates=int(snakemake.params.get("min_reps_consensus")))
