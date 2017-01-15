##Editor: Diyue Bu####
###last edition: 09/29/2015###

from copy import deepcopy
from math import log10
from sys import exit
from os import walk
from multiprocessing import Pool, Process
from random import seed, shuffle, random
import numpy as np
# from numpy import random
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def personlist():
    ##remove family numbers
    frelate = '20140625_related_individuals.txt'
    f = open(frelate, 'r')
    lremove = []
    for l in f:
        if 'Sample' in l:
            continue
        lremove.append(l.strip().split()[0])
    f.close()

    fpinfo = 'integrated_call_samples.20130502.ALL.ped.txt'
    lpersonid = []
    lrelate = set([])
    f = open(fpinfo, 'r')
    for l in f:
        if 'Population' in l:
            continue
        r = l.strip().split()
        pid = r[1]
        if pid in lremove or pid in lrelate:
            continue
        relation = r[7]
        if relation == 'child':
            pass
        elif relation == 'mother' or relation == 'father':
            if r[8] == '0':
                lpersonid.append(pid)
            else:
                pass
        else:
            lpersonid.append(pid)
    del lrelate
    f.close()

    return lpersonid


def readcol(sfile):
    popSize = 500
    # int(sfile.split('\\')[-2].split('Size')[1])

    loutput = []

    with open(sfile, 'r') as f:
        for l in f:
            l = l.strip()
            if l.startswith('##'):
                continue

            if l.startswith('#CHROM'):
                continue

            r = l.split()

            linfo = r[7].split(';')

            for i, elem in enumerate(linfo):
                if 'AC=' in elem:
                    index = i

            saf = linfo[index + 1].split('=')[1]
            if ',' in saf:
                laf = saf.split(',')
                faf = float(laf[0])
            else:
                faf = float(saf)
            # print faf

            # lsnpaf.append(faf)


            c = 0
            # check total # of snps in test set
            m = 0
            for i in range(popSize):

                if '1' in r[9 + i]:
                    m += 1
            if m < t:
                x = 0
                faf = 0.0
            elif m == t:
                freq = random()
                if freq < delta:
                    x = 0
                    faf = 0.0
                else:
                    x = 1
            else:
                x = 1

            loutput.append([x, faf])

    return loutput


# def main():
###walk file names and pair every two files to multi-process list####
# lfName = []
# for root, dirs, files in walk('vcfSize500\\vcfSize500\\'):
#     for fname in files:
#         lfName.append(fname)
#
# return lfName


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-t', '--threshold')
    parser.add_argument('-d', '--delta')
    parser.add_argument('-f', '--filename')
    parser.add_argument('-e', '--epsilon')
    parser.add_argument('-n', '--N')
    args = parser.parse_args()

    t = int(args.threshold)
    delta = float(args.delta)
    fname = args.filename
    epsilon = float(args.epsilon)
    N = int(args.N)

    # t = 1
    # delta = float(1e-02)
    cutoff = 1000000
    tmp = 0
    # epsilon = 1
    # N = 500
    # l10 = readcol('/data/diybu/challengeData16/c1beacon/wholeRecords/vcfSize500/chr10_selectedInd.vcf')
    # fname = 'vcfSize500\\vcfSize500\\chr10_selectedInd.vcf'
    l10 = readcol(fname)
    schr = fname.split('\\')[-1].split('chr')[1].split('_')[0]
    fout = open('chr' + schr + 'returnValue_laplace_delta_' + str(delta) + '.txt', 'w')
    fout.write('returnValue\tminorAlleleFreq\n')
    false_pos = 0
    false_neg = 0
    total_pos = 0
    for lpos in l10:
        if lpos[0] == 1:
            lpos[1] = lpos[1] + np.random.laplace(0, (1/epsilon))
            if lpos[1] < 0 or lpos[1] > 1:
                lpos[1] = 0.0
                lpos[0] = 0
                false_neg += 1

        fout.write(str(lpos[0]) + '\t' + str(lpos[1]) + '\n')

    fout.close()
    print('total ' + str(total_pos) + " false pos " + str(false_pos) + ' false neg ' + str(false_neg))
    print('chr' + schr + 'returnValue_laplace_delta_' + str(delta) + '.txt in output')