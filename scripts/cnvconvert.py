import argparse
import numpy as np
import pandas as pd


# all my constants
class CST:
    AMPL_REGION_TRESHOLD = 3
    AMPLF_MIN = 4
    AMPLF_MAX = 100
    # CNV_COLUMN_NAME = 'seg.mean'
    CNV_COLUMN_NAME = 'lowratio'
    CNV_COLUMN_INDEX = 7
    LEVEL_COLUMN_NAME = 'level'
    LEVEL_COLUMN_INDEX = 10
    CHR_COLUMN = 2
    ID_COLUMN = 'ID'
    LEVELS = 10
    MIN = 0
    MAX = 100


# different asymtotic functions
def f1(min, max, levels):
    l = []
    bin = (max - min) / levels
    for x in range(min, max, int(bin)):
        l.append(float(x / (x + 1)) * max)
    return l


def f2(min, max, levels):
    pass


func_dict = {"f1": f1,
             "f2": f2}


def binarysearch_left(l, value, low, high):
    if low == high:
        return low
    mid = int((low + high) / 2)
    if l[mid] == value:
        return mid
    elif l[mid] > value:
        return binarysearch_left(l, value, low, mid - 1)
    else:
        return binarysearch_left(l, value, mid + 1, high)


def compute_intervals(min=0, max=100, levels=10, func="f1"):
    """
    Generates the cnv intervals for every level.
    Assume for the start that we have 10 levels.
    The interval size follow for now an asymptotic function.
    """
    intervals = func_dict[func](min, max, levels)
    return intervals


def collapse_levels(df, levels=10):
    # count the number of bins that have same level
    l = df[[CST.ID_COLUMN, CST.LEVEL_COLUMN_NAME]].groupby(CST.LEVEL_COLUMN_NAME).count().reset_index().sort_values(
        by=CST.LEVEL_COLUMN_NAME)
    count_levels = pd.DataFrame({CST.LEVEL_COLUMN_NAME: list(range(0, levels))})
    count_levels = pd.merge(count_levels, l, how="left")
    count_levels = count_levels[CST.ID_COLUMN].fillna(0)

    last_level = 0
    for j in range(1, levels):
        if count_levels[j] > 0:
            last_level += 1
            df.loc[df[CST.LEVEL_COLUMN_NAME] == j, CST.LEVEL_COLUMN_NAME] = last_level

    return df


def compute_levels(df, min, max, levels):
    """
    Assign levels to each CNV
    """

    # initialize level
    df[CST.LEVEL_COLUMN_NAME] = 0

    levels = compute_intervals(min, max, levels)
    print(levels)
    for i in range(0, df.shape[0]):

        # consider level 0 all CNV (seg.mean) <= 3
        if df.iloc[i, CST.CNV_COLUMN_INDEX] > CST.AMPL_REGION_TRESHOLD:
            df.iloc[i, CST.LEVEL_COLUMN_INDEX] = binarysearch_left(levels, df.iloc[i, CST.CNV_COLUMN_INDEX], 0,
                                                                   len(levels))

    return df


def convert(fin, fout):
    """
    Method to convert the cnv profile into a level profile
           2
     1  |-------|  1  - levels
    -----        -----
    |                |
      30    60    29  - cnv call
    """
    print("Start converting cnv calls to levels")

    # read cnv bins
    df = pd.read_csv(fin, header=0, sep="\t")

    df = df.astype({CST.CNV_COLUMN_NAME:float})

    vmin, vmax = np.min(df[CST.CNV_COLUMN_NAME].tolist()), np.max(df[CST.CNV_COLUMN_NAME].tolist())

    print(vmin, vmax)
    # assign levels for each CNV
    df = compute_levels(df, min=CST.MIN, max=CST.MAX, levels=CST.LEVELS)
    # df = collapse_levels(df, CST.LEVELS)

    df.to_csv(fout, header=True, index_label=False, sep="\t")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fin', help='Input CNV calls (short.txt file from SmurfSeq)')
    parser.add_argument('--fout', help='Output levels')
    args = parser.parse_args()
    # convert("../data/process/AC/cnv/AC.short.txt", "../data/process/AC/cnv/AC.levels.txt")
    # convert(args.fin, args.fout)
