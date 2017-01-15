#!/usr/bin/env python


import bitarray as ba
import pickle as pl
import random
import numpy as np
from multiprocessing import Pool, Process
from argparse import ArgumentParser


def returnResult(schr):
    lreturn = []

    # sfile = 'chr' + str(schr) + 'returnValue_500wholeR.txt'
    sfile = rare_fname
    print('load beacondb results', sfile)
    f = open(sfile, 'r')
    n1 = 0
    total_data = 0
    for l in f:
        if l.startswith('return'):
            continue
        lreturn.append(l.strip().split('\t'))
        if l.strip().split('\t')[0] == '1':
            n1 += 1
        total_data += 1
    f.close()

    return lreturn, n1, total_data


def main(sIn, frequency, positives, total_records):
    lid = []
    ######################################################
    # change the number of individuals in input vcf file #
    ######################################################
    lRes = [ba.bitarray() for i in range(N)]
    count = 0
    tmp = 0
    print('opening file ', sIn, ' epsilon ', epsilon)
    with open(sIn, 'r') as f:
        for l in f:
            if l.startswith('#'):
                continue

            count += 1
            r = l.split()
            index = 0
            # epsilon = 5
            t = 1
            total_data = int(N * positives / total_records)
            p = np.random.laplace(0, (2 * total_data / epsilon))
            for gt in r[9:]:
                vi = np.random.laplace(0, (4 * total_data / epsilon))

                if '1' in gt and tmp < total_data:
                    # print(str(vi) + ":" + str(p))

                    if frequency[index][0] == '1':
                        t = np.random.laplace(0, (2 * total_data / epsilon))
                    else:
                        t = 1
                    if vi >= p + t:  # q(d)=1
                        lRes[index].append(True)
                        p = np.random.laplace(0, (2 * total_data / epsilon))
                        tmp += 1
                        # if tmp > total_data:
                        #     break
                    else:
                        lRes[index].append(False)
                else:
                    if vi >= p:
                        lRes[index].append(False)
                    else:
                        lRes[index].append(True)
                        p = np.random.laplace(0, (2 * total_data / epsilon))
                        tmp += 1

                index += 1

    print('positive values ' + str(tmp))
    fname = sIn.split('\\')[-2]
    output = open(fname + 'WholeRecords_sparse_vector_limited_' + str(epsilon) + '.pkl', 'wb')
    print(fname + 'WholeRecords_sparse_vector_limited_' + str(epsilon) + '.pkl output')
    pl.dump(lRes, output, -1)
    output.close()


if __name__ == '__main__':
    ##########################
    #  change input vcf file #
    ##########################
    fname = 'vcfSize500\\vcfSize500\\chr10_selectedInd.vcf'
    parser = ArgumentParser()
    # parser.add_argument('-f', '--filename')
    # parser.add_argument('-n', '--N')
    # parser.add_argument('-e', '--Epsilon')
    # parser.add_argument('-r', '--rare_file')
    args = parser.parse_args()

    fname = args.filename
    # N = int(args.N)
    # epsilon = float(args.N)
    rare_fname = args.rare_file
    N = 500
    t = 1
    epsilon = float(1e-06)
    frequency, n1, total_records = returnResult(rare_fname)
    main(fname, frequency, n1, total_records)
