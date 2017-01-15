#!/usr/bin/env python


import bitarray as ba
import pickle as pl
import random
from multiprocessing import Pool, Process
from argparse import ArgumentParser


def returnResult(schr):
    lreturn = []

    sfile = 'chr' + str(schr) + 'returnValue_500wholeR.txt'
    print('load beacondb results', sfile)
    f = open(sfile, 'r')
    n1 = 0
    for l in f:
        if l.startswith('return'):
            continue
        lreturn.append(l.strip().split('\t'))
        if l.strip().split('\t')[0] == '1':
            n1 += 1
    f.close()

    return lreturn


def main(sIn):
    lid = []
    ######################################################
    # change the number of individuals in input vcf file #
    ######################################################
    lRes = [ba.bitarray() for i in range(N)]
    count = 0
    tmp = 0
    print('opening file ', sIn, ' bias ', rand_bias)
    with open(sIn, 'r') as f:
        for l in f:
            if l.startswith('#'):
                continue

            count += 1
            r = l.split()
            index = 0

            for gt in r[9:]:
                rand = random.randint(0, 100)
                if '1' in gt:
                    # rand = random.randint(1, 100)
                    # if
                    #     rand -= random.randint(rand_bias, 100)
                    # if rand > rand_bias:
                    #     lRes[index].append(True)
                    #     tmp += 1
                    # else:
                    # if frequency[index][0] == '1':
                    #     if rand > 80:
                    #         lRes[index].append(True)
                    #         tmp += 1
                    #     else:
                    #         lRes[index].append(False)
                    # else:

                    if rand > rand_bias:
                        lRes[index].append(True)
                        tmp += 1
                    else:
                        rand = random.randint(0, 100)
                        if rand > rand_bias:
                            lRes[index].append(True)
                            tmp += 1
                        else:
                            lRes[index].append(False)
                else:
                    if rand > rand_bias:
                        lRes[index].append(False)
                    else:
                        rand = random.randint(0, 100)
                        if rand > rand_bias:
                            lRes[index].append(False)
                        else:
                            lRes[index].append(True)
                            tmp += 1

                index += 1

    print('positive values' + str(tmp))
    fname = sIn.split('\\')[-2]
    output = open(fname + 'WholeRecords_rand_res_bias_' + str(rand_bias) + '_both.pkl', 'wb')
    print(fname + 'WholeRecords_rand_res_bias_' + str(rand_bias) + '_both.pkl output')
    pl.dump(lRes, output, -1)
    output.close()


if __name__ == '__main__':
    ##########################
    #  change input vcf file #
    ##########################
    # fname = 'vcfSize500\\vcfSize500\\chr10_selectedInd.vcf'
    parser = ArgumentParser()
    parser.add_argument('-f', '--filename')
    parser.add_argument('-n', '--N')
    parser.add_argument('-b', '--rand_bias')
    args = parser.parse_args()

    fname = args.filename
    N = int(args.N)
    # N = 500
    rand_bias = int(args.rand_bias)  # bias < 50 then more true results
    main(fname)
