#!/usr/bin/env python


# import bitarray as ba
import pickle
from argparse import ArgumentParser

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-f', '--main_file')
    parser.add_argument('-n', '--solution_file')
    args = parser.parse_args()

    fname = args.main_file
    solution_file = args.solution_file
    print()
    # fname = 'vcfSize500wholeRecords.pkl'
    # solution_file = 'vcfSize500WholeRecords_rand_res_bias_2_both.pkl'
    with open(fname, 'rb') as pkl_file:
        pkl_file = open(fname, 'rb')
        LcolIn = pickle.load(pkl_file)
        pkl_file.close()
        with open(solution_file, 'rb') as pkl_file2:
            pkl_file2 = solution_file
            print(pkl_file2)
            pkl_file = open(pkl_file2, 'rb')
            LcolMethod = pickle.load(pkl_file)
            pkl_file.close()
            total_record = 0
            ones = 0
            zero_one = 0
            one_zero = 0
            # print(str(len(LcolMethod[0])))
            for i in range(500):
                smaller = len(LcolIn[i]) if len(LcolMethod[i]) > len(LcolIn[i]) else len(LcolMethod[i])
                # print(str(smaller))

                total_record += len(LcolIn[i])
                for j in range(smaller):
                    if LcolIn[i][j] == 1:
                        ones += 1
                    if LcolIn[i][j] > LcolMethod[i][j]:
                        one_zero += 1
                    elif LcolIn[i][j] < LcolMethod[i][j]:
                        zero_one += 1

            # print('zero_one ', zero_one, ' one_zero ', one_zero, ' total ', total_record, ' ones ', ones, 'ones remaining',
            #       str(ones - one_zero))
            true_positive = ones - one_zero
            true_negetive = total_record - ones - zero_one
            false_positive = zero_one
            false_negetive = one_zero
            print("true_positive ", true_positive, " true_negetive ", true_negetive, " false_positive ", false_positive,
                  " false_negetive ", false_negetive, " total ", total_record)
            print("sensitivity ", true_positive / (true_positive + false_negetive), " specificity ",
                  true_negetive / (true_negetive + false_positive))
