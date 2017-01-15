#!/usr/bin/env python


import bitarray as ba
import pickle as pl
from multiprocessing import Pool, Process
from argparse import ArgumentParser


def main(sIn):
    lid = []
    ######################################################
    # change the number of individuals in input vcf file #
    ######################################################
    lRes = [ba.bitarray() for i in range(N)]
    count = 0
    tmp = 0
    with open(sIn, 'r') as f:
        # lines = f.readlines()
        for l in f:
            # print(l)
            if l.startswith('#'):
                continue
            # print(l)
            count += 1
            r = l.split()
            index = 0
            for gt in r[9:]:
                if '1' in gt:
                    lRes[index].append(True)
                    tmp += 1
                else:
                    lRes[index].append(False)
                # if lRes[index] == True:
                index += 1
    print('positive ' + str(tmp))
    fname = sIn.split('\\')[-2]
    output = open(fname + 'wholeRecords.pkl', 'wb')
    print(fname + 'wholeRecords.pkl output')
    pl.dump(lRes, output, -1)
    output.close()


if __name__ == '__main__':
    ##########################
    #  change input vcf file #
    ##########################
    fname = 'vcfSize500\\vcfSize500\\chr10_selectedInd.vcf'
    # parser = ArgumentParser()
    #    parser.add_argument('-f', '--filename')
    #    parser.add_argument('-n', '--N')
    #    args = parser.parse_args()

    # fname = args.filename
    # N = int(args.N)
    N = 500

    main(fname)
