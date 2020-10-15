import argparse
import os

def get_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("-o","--outdir",type=str,default="output")
    parser.add_argument("-d","--path",type=str,default=None,help="path of cellranger count or aggr")
    parser.add_argument("-t","--tcr",type=str,default=None,help="path of tcr")
    parser.add_argument("-n","--n_jobs",type=int,default=8,help="number of jobs")
    args=parser.parse_args()
    return args
