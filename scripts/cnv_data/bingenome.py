import pandas
import collections
import numpy as np

def bin_genome(file, binsize=5):
    fout = """hg38_bins_{}k_size.bed""".format(binsize)
    binsize = binsize * 1000  # scale to kb

    exclude_suffix = ("alt", "random", "v1", "v2")

    total_length = 0
    dict = {}

    with open(fout, "w") as g:
        with open(file, "r") as f:
            for l in f.readlines():
                chr, stop = l.strip("\n").split("\t")
                stop = int(stop)
                # include only primary chromosomes
                if chr.endswith(exclude_suffix) is False:
                    for start in range(0, stop, binsize):
                        g.write("""{}\t{}\t{}\n""".format(chr, str(start), str(int(min(start + binsize, stop)))))

                    # compute total length
                    total_length += stop
                    dict[chr] = (0, stop)


def convert_absolute2relative_genome(mappability, chromsize):
    offset = 0
    exclude_suffix = ("alt", "random", "v1", "v2")
    ordering = ["chr1",
                "chr2",
                "chr3",
                "chr4",
                "chr5",
                "chr6",
                "chr7",
                "chr8",
                "chr9",
                "chr10",
                "chr11",
                "chr12",
                "chr13",
                "chr14",
                "chr15",
                "chr16",
                "chr17",
                "chr18",
                "chr19",
                "chr20",
                "chr21",
                "chr22",
                "chrX",
                "chrY"
                ]
    dict = {}
    total = 0
    with open(chromsize, "r") as f:
        for l in f.readlines():
            chr, stop = l.strip("\n").split("\t")
            stop = int(stop)
            # include only primary chromosomes
            if chr.endswith(exclude_suffix) is False:
                dict[chr] = (0, stop)
                total += stop

    for c in ordering:
        old_offset = offset
        offset = dict[c][1] + offset
        dict[c] = old_offset

    with open(mappability + "_relative.bed", "w") as g:
        with open(mappability, "r") as f:
            for l in f.readlines():
                chr, start, stop, map = l.strip().split("\t")
                g.write("""{}\t{}\t{}\t{}\n""".format(chr,
                                                      str(int(start) - dict[chr]),
                                                      str(int(stop) - dict[chr]), map))

def collapse_gccontent(gccontentfile, nbins=1000):
    """
    Reduce number of bins for gccontenfile
    """

    with open(gccontentfile+"_reduced","w") as g:
        lstart = -1
        lstop = -1
        lchr = -1
        lgc = []
        count = 0

        for l in open(gccontentfile):
            count += 1
            print(count)
            chr, start, stop, gc = l.strip().split("\t")

            if count == nbins:
                mean = np.mean(lgc)
                g.write("""{}\t{}\t{}\t{}\n""".format(lchr, lstart, lstop, str(mean)))
                lchr = chr
                lstart = start
                lstop = stop
                lgc = []
                lgc.append(float(gc))
                count = 0
                continue

            if chr == lchr:
                lstop = stop
                lgc.append(float(gc))
                continue

            if chr != lchr and lchr == -1:
                lchr = chr
                lstart = start
                lstop = stop
                lgc.append(float(gc))
                continue

            # write to file
            if chr != lchr:
                g.write("""{}\t{}\t{}\t{}\n""".format(lchr, lstart, lstop, str(np.mean(lgc))))
                lchr = chr
                lstart = start
                lstop = stop
                lgc = []
                lgc.append(float(gc))
                count = 0

        g.write("""{}\t{}\t{}\t{}\n""".format(lchr, lstart, lstop, str(np.mean(lgc))))

def compress_gccontent(file):

    with open(file+"_reduced","w") as g:
        lchr = -1
        lgc = []
        larm = -1
        for l in open(file):
            chr, gc, arm = l.strip().split("\t")
            gc = float(gc)

            if chr != lchr and lchr == -1:
                lchr = chr
                lgc.append(gc)
                larm = arm
                continue

            if chr != lchr:
                print(lchr)
                g.write("""{}\t{}\t{}\n""".format(lchr, str(np.mean(lgc)), larm))
                lchr = chr
                lgc = []
                lgc.append(gc)
                larm = arm
        g.write("""{}\t{}\t{}\n""".format(lchr, str(np.mean(lgc)), larm))


if __name__ == "__main__":
    # bin_genome("hg38.chrom.sizes", binsize=5)
    # bin_genome("hg38.chrom.sizes", binsize=50)
    # convert_absolute2relative_genome("map_hg38_500kb.wig.bed", "hg38.chrom.sizes")
    # compute_gcfile("hg38_bins_5k_size.bed", "hg38.gc5Base.bed")
    # collapse_gccontent("hg38.gc5Base.bed")
    compress_gccontent("hg38_bins_50k_size_gccontent.txt")
    compress_gccontent("hg38_bins_5k_size_gccontent.txt")
