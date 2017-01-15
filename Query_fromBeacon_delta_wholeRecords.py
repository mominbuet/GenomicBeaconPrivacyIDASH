##########################USAGE############################################################
# this script is to perform allele frequency attack based on log likelihood ratio test(lrt)#
# python Query_fromBeacon_delta.py -t returnThreshold -d delta				###
# input:  case(not in beacon) & test(in beacon) individual variation info		###
#	 beacon return result & allele freq (prepared by user)				###
#        										###
# output: lrt scores for each individual every 100 queries				###
#											###
###########################################################################################

# from pandas import read_csv
# import pandas as pd
from copy import deepcopy
from math import log10, pow
from sys import exit
from os import walk
from multiprocessing import Pool, Process
from random import seed, shuffle
from numpy import random, array, mean
from argparse import ArgumentParser
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
import pickle


def lrtTest(lfName):
    fnameIn = lfName[0]
    fnameOut = lfName[1]
    t = lfName[2]
    delta = lfName[3]
    print(fnameIn)
    ####################################
    #   customized input for 2 dataset #
    #   case: individuals not in beacon#
    #   test: individuals in beacon    #
    ####################################
    # l_caseIndex = range(500, 1000)[:100]
    # l_testIndex = range(0, 500)[:100]

    # N = 500

    print('population size:', N)

    ################
    ##set chr NO. ##
    ################
    schr = '10'

    ##############################
    ##returnResult from beacondb##
    ## 1st col: 0/1; 2nd col: af##
    ##############################
    lreturn = returnResult(t, schr, delta)

    ######################################################################################
    ##load preprocessed variation info per person					    ##
    ##need to be preprocessed by user            					    ##
    ##each element in Lcol is a list of variation distribution (True) of an individual  ##
    ######################################################################################
    with open(fnameIn, 'rb') as pkl_file:
        # pkl_file = open(fnameIn, 'rb')
        LcolIn = pickle.load(pkl_file)
        pkl_file.close()
        # pkl_file = open(fnameOut, 'rb')
        with open(fnameOut, 'rb') as pkl_file:
            LcolOut = pickle.load(pkl_file, encoding='latin1')
            pkl_file.close()

            Lrtcase = []
            lncase = []
            ##calculate lrt score for case dataset(not in beacon)
            for i in range(N):
                lperson = LcolOut[int(i)]
                # print 'per individual lcol length', len(lperson)
                lrt_case = cal_lrt(lreturn, lperson, N, delta)

                Lrtcase.append(lrt_case)
                lrt_case = []

            Lrttest = []
            lntest = []
            ##calculate lrt score for test dataset(in beacon)
            for i in range(N):
                lperson = LcolIn[int(i)]
                # print 'per individual lcol length', len(lperson)
                lrt_test = cal_lrt(lreturn, lperson, N, delta)
                Lrttest.append(lrt_test)
                lrt_test = []

            print('Lrtcase/test length:', len(Lrtcase), len(Lrttest))

            save_result(Lrtcase, Lrttest, t, schr)


def cal_lrt(lreturn, lcol, N, delta):
    # print(str((int)(len(lcol) / 100)))
    ##calculate lrt score every 100 queries
    lrt_step = [0.0 for i in range((int)(len(lcol) / 100))]
    count = 0
    lrt = 0.0
    psuedo = float(pow(10, -290))
    # delta = float(pow(10, -6))

    ntotal = 0
    n1 = 0
    # print(count / 100 == len(lrt_step))
    for i in range((int)(len(lcol))):
        if count / 100 == len(lrt_step):
            break
        if i > len(lreturn):
            print(str(i) + ":" + str(len(lreturn)) + ":" + str(len(lcol)))
            break
        if lcol[i] == True:
            ntotal += 1
            linfo = lreturn[i]
            x = int(linfo[0])
            faf = float(linfo[1])

            if x == 1:
                n1 += 1

            D1 = dprob1(faf, t, N, delta)
            D0 = dprob0(faf, t, N)
            if D1 == 0.0:
                LH1 = x * log10(1 - D1 - psuedo) + (1 - x) * log10(D1 + psuedo)
            elif D1 == 1.0:
                LH1 = x * log10(1 - D1 + psuedo) + (1 - x) * log10(D1 - psuedo)
            else:
                LH1 = x * log10(1 - D1) + (1 - x) * log10(D1)
            if D0 == 0.0:
                LH0 = x * log10(1 - D0 - psuedo) + (1 - x) * log10(D0 + psuedo)
            elif D0 == 1.0:
                LH0 = x * log10(1 - D0 + psuedo) + (1 - x) * log10(D0 - psuedo)
            else:
                LH0 = x * log10(1 - D0) + (1 - x) * log10(D0)
            L = LH0 - LH1

            nstep = (int)(count / 100)
            # print(nstep)
            lrt_step[nstep] = lrt_step[nstep] + L
            count += 1

    lrt_step = list(filter(lambda x: x != 0.0, lrt_step))

    return lrt_step


def dprob1(f, t, N, delta):
    if t == 1:
        D = delta * pow((1 - f), (2 * (N - 1)))
        return D

    D = delta * pow((1 - f), (2 * (N - t)))
    for i in range(1, t):
        D += pow(f, (i - 1)) * pow((1 - f), (2 * (N - i)))

    return D


def dprob0(f, t, N):
    D = 0
    for i in range(t):
        D += pow(f, i) * pow((1 - f), (2 * (N - i)))

    return D


def returnResult(t, schr, delta):
    lreturn = []

    sfile = rare_file
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

    print('# of queries beacon will return True', n1)
    print('size of beacondb', len(lreturn))

    return lreturn


def main(t, delta):
    ###walk file names and pair every two files to multi-process list####
    lfName = []
    for root, dirs, files in walk('query_pickle/'):
        for fname in files:
            lfName.append([fname, t, delta])

    return lfName


def save_result(Lcase, Ltest, t, schr):
    ###################################
    # customize output file dir & name#
    ###################################
    sfile = 'chr' + schr + 'lrtScoreByStep_laplace_delta_' + str(delta) + '.txt'
    print(sfile + ' output')
    fout = open(sfile, 'w')
    fout.write('#query different length of genomes\n')
    fout.write('#Not in beacon lrt score(case)\tIn beacon lrt score(test)\n')
    length = min(len(Lcase), len(Ltest))

    # print len(Lcase), len(Ltest)

    step = len(Lcase[0])
    # print step

    b = 0
    for i in range(length):
        fout.write('###a pair of individuals lrt score by 100 queries per step\n')
        for j in range(len(Lcase[i])):
            try:
                fout.write(str(Lcase[i][j]) + '\t' + str(Ltest[i][j]) + '\n')
            except:
                pass

    fout.close()


if __name__ == '__main__':
    parser = ArgumentParser()
    # parser.add_argument('-n', '--N')
    # parser.add_argument('-t', '--threshold')
    # parser.add_argument('-d', '--delta')
    # parser.add_argument('-f1', '--file1')
    # parser.add_argument('-f2', '--file2')
    # parser.add_argument('-r', '--rare_file')
    args = parser.parse_args()

    # N = int(args.N)
    # t = int(args.threshold)
    # delta = float(args.delta)
    # f1 = args.file1
    # f2 = args.file2
    # rare_file = args.rare_file
    t = int(1)
    delta = float(1e-02)
    N = 500
    f1 = "vcfSize500wholeRecords.pkl"
    f2 = "500notInBeaconwholeRecords.pkl"
    rare_file = "chr10returnValue_laplace_delta_1e-06.txt "
    ###############################
    # perform query attack        #
    ###############################
    # use chr10 as example

    lrtTest([f1, f2, t, delta])
