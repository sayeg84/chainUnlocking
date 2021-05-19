import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from distutils.util import strtobool as sb


parser = argparse.ArgumentParser(description = "Make plots for results column")

parser.add_argument("--path",type=str,help="path of the file to plot",required=True)
parser.add_argument("-H",action="store_true")

args = parser.parse_args()

if args.H:
    df = pd.read_csv(args.path,header=0)
else:
    df = pd.read_csv(args.path,header=None)
    df = df.rename(columns = {0:"ls",1:"T_m",2:"T_std",3:"rmsd_m",4:"rmsd_std"})
print(df.head())
for s in ["T","rmsd"]:
    fig = plt.figure(figsize=(8,4))
    plt.errorbar(df["ls"],df["{0}_m".format(s)],yerr=df["{0}_std".format(s)],marker="o",markeredgecolor="black")
    plt.grid()
    plt.show()
