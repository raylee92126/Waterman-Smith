#!/usr/bin/python

__author__ = "Raymond Lee"
__email__ = "raymond92126@gmail.com"
__copyright__ = "None"
__license__ = "MIT"
__version__ = "0.0.1"

import argparse
import sys
import numpy as np
import pandas as pd


## Usage: python hw1.py -i <input file> -s <score file> -op <output file>
## Example: python hw1.py -i input.txt -s blosum62.txt -op outputfile.txt
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-op', '--output', help='output file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False, default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False, default=-1)
args = parser.parse_args()

def writedf(df, opened_file):
    colnames_list = df.columns.to_list()
    printed_colnames = "\t".join(map(str, colnames_list))
    opened_file.write(printed_colnames)
    opened_file.write("\n")
    
    for i, row in df.iterrows():
        printed_row = "\t".join(map(str, row.tolist()))
        opened_file.write(printed_row)
        opened_file.write("\n")

def printdf(df, opened_file):
    colnames_list = df.columns.to_list()
    printed_colnames = "\t".join(map(str, colnames_list))
    print(printed_colnames)
    
    for i, row in df.iterrows():
        printed_row = "\t".join(map(str, row.tolist()))
        print(printed_row)


def runSW(inputFile, scoreFile, outputFile, openGap, extGap):
    ### load scoreFile
    score_df = pd.read_csv(scoreFile, header=0, delim_whitespace=True, index_col=0)

    #load inputFile
    with open(inputFile, "r") as f:
        seq1_in, seq2_in = f.read().split("\n")
    seq1_list = list(seq1_in) #on top of matrix
    seq2_list = list(seq2_in) #on side of matrix

    #initialize necessary arrays (main_arr, left_arr, up_arr, dir_arr)
    main_arr = np.zeros([(len(seq2_list)+1),(len(seq1_list)+1)], dtype="int64")
    main_arr[0, 0] = 0
    main_arr[0, 1] = 0
    main_arr[1, 0] = 0
    for j in range(1, len(seq1_list)):
        main_arr[0, j+1] = 0
    for i in range(1, len(seq2_list)):
        main_arr[i+1, 0] = 0
        
    left_arr = np.zeros([(len(seq2_list)+1),(len(seq1_list)+1)], dtype="int64")
    left_arr[0, 0] = 0
    left_arr[0, 1] = 0
    left_arr[1, 0] = 0
    for j in range(1, len(seq1_list)):
        left_arr[0, j+1] = 0
    for i in range(1, len(seq2_list)):
        left_arr[i+1, 0] = 0

    up_arr = np.zeros([(len(seq2_list)+1),(len(seq1_list)+1)], dtype="int64")
    up_arr[0, 0] = 0
    up_arr[0, 1] = 0
    up_arr[1, 0] = 0
    for j in range(1, len(seq1_list)):
        up_arr[0, j+1] = 0
    for i in range(1, len(seq2_list)):
        up_arr[i+1, 0] = 0

    dir_arr = np.empty([(len(seq2_list)+1),(len(seq1_list)+1)], dtype='object')
    dir_arr[0, 0] = "stop"
    for j in range(len(seq1_list)):
        dir_arr[0, j+1] = "stop"
    for i in range(len(seq2_list)):
        dir_arr[i+1, 0] = "stop"

    #run algorithm for filling in arrays
    for i, bp2 in enumerate(seq2_list):
        i += 1
        for j, bp1 in enumerate(seq1_list):
            j += 1
            
            left_arr[i, j] = max(left_arr[i, j-1] + extGap, #extend existing gap
                                main_arr[i, j-1] + openGap, #create new gap
                                0) 
            up_arr[i, j] = max(up_arr[i-1, j] + extGap, #extend existing gap
                            main_arr[i-1, j] + openGap, #create new gap
                            0) 
            match_score = main_arr[i-1, j-1] + score_df.loc[bp2, bp1]
            
            max_val = max(0, left_arr[i, j], up_arr[i, j], match_score)
            main_arr[i, j] = max_val
            
            if max_val==0:
                dir_arr[i, j] = "stop"
            elif max_val==match_score:
                dir_arr[i, j] = "diag"
            elif max_val==left_arr[i, j]:
                dir_arr[i, j] = "left"
            elif max_val==up_arr[i, j]:
                dir_arr[i, j] = "up"

    #backtracing
    start_indicies = np.where(main_arr == np.amax(main_arr))
    alignment_score = np.amax(main_arr)
    curr_row = start_indicies[0][0]
    curr_col = start_indicies[1][0]
    seq1 = []
    seq2 = []

    while curr_row>0 and curr_col>0:
        if dir_arr[curr_row, curr_col] == "stop":
            break
        elif dir_arr[curr_row, curr_col] == "diag":
            seq1.append(seq1_list[curr_col-1])
            seq2.append(seq2_list[curr_row-1])
            curr_row-=1
            curr_col-=1
        elif dir_arr[curr_row, curr_col] == "left":
            seq1.append(seq1_list[curr_col-1])
            seq2.append("-")
            curr_col -= 1
        elif dir_arr[curr_row, curr_col] == "up":
            seq1.append("-")
            seq2.append(seq2_list[curr_row-1])
            curr_row-=1


    ### parsing output sequences
    seq1.reverse()
    seq2.reverse()

    bars_aligned = []
    for seq1_i, seq2_i in zip(seq1, seq2):
        if seq1_i==seq2_i:
            bars_aligned.append("|")
        else:
            bars_aligned.append(" ")

    seq1_prefix = seq1_list[0:(curr_col)]
    seq2_prefix = seq2_list[0:(curr_row)]
    if len(seq1_prefix) >= len(seq2_prefix):
        seq1_pre_ws = []
        bar_pre_ws = [" " for _ in range(len(seq1_prefix))]
        seq2_pre_ws = [" " for _ in range((len(seq1_prefix) - len(seq2_prefix)))]
    else:
        seq1_pre_ws = [" " for _ in range((len(seq2_prefix) - len(seq1_prefix)))]
        bar_pre_ws = [" " for _ in range(len(seq2_prefix))]
        seq2_pre_ws = []

    seq1_suffix = seq1_list[start_indicies[1][0]:]
    seq2_suffix = seq2_list[start_indicies[0][0]:]
    if len(seq1_suffix) >= len(seq2_suffix):
        seq1_post_ws = []
        bar_post_ws = [" " for _ in range(len(seq1_suffix))]
        seq2_post_ws = [" " for _ in range((len(seq1_suffix) - len(seq2_suffix)))]
    else:
        seq1_post_ws = [" " for _ in range((len(seq2_suffix) - len(seq1_suffix)))]
        bar_post_ws = [" " for _ in range(len(seq2_suffix))]
        seq2_post_ws = []


    seq1_out = "".join(seq1_pre_ws + seq1_prefix + ["("] + seq1 + [")"] + seq1_suffix + seq1_post_ws)
    bars_out = "".join(bar_pre_ws + [" "] + bars_aligned + [" "] + bar_post_ws)
    seq2_out = "".join(seq2_pre_ws + seq2_prefix + ["("] + seq2 + [")"] + seq2_suffix + seq2_post_ws)

    #printing output
    output1 = """-----------
|Sequences|
-----------
sequence1
{input_seq1}
sequence2
{input_seq2}
--------------
|Score Matrix|
--------------
""".format(input_seq1=seq1_in, input_seq2=seq2_in)

    #print score matrix as tsv
    main_df = pd.DataFrame(main_arr)
    main_df.columns = [""] + seq1_list
    main_df.insert(0, "", [""]+seq2_list, allow_duplicates=True)



    #print second part of output
    output2 = """----------------------
|Best Local Alignment|
----------------------
Alignment Score:{max_score}
Alignment Results:
{output_seq1}
{bars}
{output_seq2}
""".format(max_score=alignment_score, output_seq1=seq1_out, bars=bars_out, output_seq2=seq2_out)
    
    # print(output1)
    # #printdf(main_df)
    # print(output2)
    with open(outputFile, "w") as outfile:
        outfile.write(output1)
        writedf(main_df, outfile)
        outfile.write(output2)


### Run your Smith-Waterman Algorithm
if __name__ == "__main__":
    runSW(args.input, args.score, args.output, args.opengap, args.extgap)
