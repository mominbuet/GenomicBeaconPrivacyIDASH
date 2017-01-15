####Editor: Diyue Bu, Indiana University####

from copy import deepcopy
from math import log10
from sys import exit
from os import walk, mkdir, path
from multiprocessing import Pool, Process
from random import seed, shuffle, random
# from numpy import random
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from subprocess import call


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


def writeVCFbyN(sfile):
    # print 'population size:', N
    # print sfile
    lpersonid = personlist()
    lindex = []
    print(sfile)
    schr = sfile.split('.')[0].split('chr')[1]
    fout = open('output\\vcfSize' + str(N) + '\\chr' + schr + '.vcf', 'w')

    crare = 0
    ctotal = 0
    with open(sfile, 'r') as f:
        for l in f:
            l = l.strip()
            if l.startswith('##'):
                lnew = l + '\n'
                fout.write(lnew)
                continue

            if l.startswith('#CHROM'):
                lpid = l.split()[9:]
                seed(99)
                lnofamily = list(set(lpid).intersection(set(lpersonid)))
                shuffle(lnofamily)
                ntotal = len(lnofamily)
                lselectid = lnofamily[:N]
                lnew = '\t'.join(map(str, l.split()[:9]))
                lnew = lnew + '\t' + '\t'.join(map(str, lselectid)) + '\n'
                fout.write(lnew)
                for i in lselectid:
                    lindex.append(lpid.index(i))
                del lselectid
                del lpid
                continue

            # r = l.split()
            # c1 = 0
            # lgt = r[:9]
            # for i in lindex:
            #     gt = r[9 + i]
            #     if '.' in gt:
            #         gt = '0|0'
            #     elif '1|1' in gt:
            #         gt = '0|1'
            #         c1 += 1
            #     elif '2' in gt:
            #         gt = '0|0'
            #     elif '0|1' in gt or '1|0' in gt:
            #         c1 += 1
            #     else:
            #         pass
            #     lgt.append(gt)
            #
            # if c1 == 0:
            #     continue
            # ctotal += 1
            # if c1 == 1:
            #     crare += 1
            # # lgt[7] = 'beaconAF=' + str(float(c1)/len(lindex))
            # lnew = '\t'.join(lgt)
            # lnew = lnew + '\n'
            # fout.write(lnew)

    f.close()

    fout.close()

    print
    'rare snp in chr', schr, ':', crare
    print
    'total snp in chr', schr, ':', ctotal


def main():
    ###walk file names and pair every two files to multi-process list####
    lfName = []
    for root, dirs, files in walk('vcfSize500\\vcfSize500'):
        for fname in files:
            lfName.append('vcfSize500\\vcfSize500\\' + fname)

    return lfName


if __name__ == '__main__':
    # parser = ArgumentParser()
    # parser.add_argument('-n', '--N')
    # args = parser.parse_args()

    N = int(500)

    lfName = main()
    writeVCFbyN(lfName[1])

    # p = Pool(6)
    # p.map(writeVCFbyN, lfName)
