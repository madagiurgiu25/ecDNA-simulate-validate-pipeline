import numpy
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import *
from collections import defaultdict

import subprocess
import sys
import os

def create_genome(bedfile, reference, outfasta):
    """
    Generates a circulome reference
    """

    # read names and postions from bed file
    positions = defaultdict(list)
    circles = defaultdict(list)
    with open(bedfile) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chr, start, stop, strand,  target, coverage, composition, fragment = line.split()
            positions[chr].append((int(start), int(stop)))
            circles[target].append((chr, int(start), int(stop), strand, int(coverage), composition, fragment))

    # parse faste file and turn into dictionary
    records = SeqIO.to_dict(SeqIO.parse(open(reference), 'fasta'))

    # search for short sequences
    circleseq_records = []
    for kcirc in circles:

        sequence_circle = ""
        id = ""
        description = ""

        for (chr, start, stop, strand, coverage, composition, fragment) in circles[kcirc]:

            # get sequence
            long_seq_record = records[chr]
            long_seq = long_seq_record.seq

            # get metadata
            description = composition

            short_seq = str(long_seq)[start-1:stop]
            if strand == "-":
                short_seq = str(Seq(short_seq).reverse_complement())
                # short_seq = short_seq[::-1]
            sequence_circle += short_seq.replace("N","A") # replace by default all N's with A's

        # extend reference to close the circle
        # start, end = sequence_circle[:5000], sequence_circle[-5000:]
        # sequence_circle = end + sequence_circle + start
        circle_record = SeqRecord(Seq(sequence_circle), id=kcirc, description=description)
        circleseq_records.append(circle_record)

    # write to file
    with open(outfasta, 'w') as f:
        SeqIO.write(circleseq_records, f, 'fasta')

def mapping(fastq, reference):

    sam = "circles.sam"
    bam = "circles.bam"
    result = subprocess.call("""ngmlr --bam-fix --threads 8 --reference {} --query {} --output {} --presets ont""".format(reference, fastq, sam), shell=True)
    result = subprocess.call("""samtools sort -O BAM -o {} {}""".format(bam, sam), shell=True)
    result = subprocess.call("""samtools index {} && rm {}""".format(bam, sam), shell=True)

def simulate_reads_circular(bedfile, reference, fastqdir):

    if not os.path.exists(fastqdir):
        os.makedirs(fastqdir)
    os.chdir(fastqdir)

    file_circle_reference = "circles_reference.fa"
    file_fastq = "circles.fastq"

    print("Create circulome fasta")
    create_genome(bedfile, reference, file_circle_reference)

    # simulate reads
    print("Simulate reads ")
    result = subprocess.call("""pbsim --circular 1 --depth 100 --seed 500 --hmm_model /Users/madag/Projects/tools/pbsim2/data/R94.model {}""".format(file_circle_reference), shell=True)
    result = subprocess.call("""cat *.fastq > {} && rm sd_*""".format(file_fastq), shell=True)


# def simulate_reads(bedfile, reference, fastqdir):
#     """
#     Generates simulate reads
#     """
#     if not os.path.exists(fastqdir):
#         os.makedirs(fastqdir)
#     os.chdir(fastqdir)

#     file_circle_frame = "circles_frames.fa"
#     file_circle_reference = "circles_reference.fa"
#     file_fastq = "circles.fastq"

#     print("Create circulome fasta")
#     circleseq_records = create_genome(bedfile, reference, file_circle_reference)

#     temp_ref = []
#     for circle_ref in circleseq_records:
#         print("Roll over ", circle_ref.id)
#         # roll over the circle sequence in 1000bp window and generate frames of the circle
#         # to cover better the breakpoints
#         temp_ref.append(circle_ref)
#         sequence = circle_ref.seq
#         for i in range(0,10):
#             # put i*1000bp to end
#             temp_sequence = sequence[(i+1)*1000:] + sequence[:(i+1)*1000]
#             circle_frame = SeqRecord(Seq(temp_sequence), id=circle_ref.id+"_"+str(i), description=circle_ref.description)
#             temp_ref.append(circle_frame)
    
#     with open(file_circle_frame, 'w') as f:
#         SeqIO.write(temp_ref, f, 'fasta')

#     # simulate reads
#     print("Simulate reads ", circle_ref.id)
#     result = subprocess.call("""pbsim --depth 30 --seed 500 --hmm_model /Users/madag/Projects/tools/pbsim2/data/R94.model {}""".format(file_circle_frame), shell=True)
#     result = subprocess.call("""cat *.fastq > {} && rm sd_*""".format(file_fastq), shell=True)

#     # mapp
#     mapping(reference, file_fastq)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate circular genome')
    parser.add_argument('-i','--input', help='Circles configuration (BED format)', required=True)
    parser.add_argument('-g','--genome', help='Reference genome', required=True)
    # parser.add_argument('-f','--fasta', help='Circles reference genome', required=True)
    parser.add_argument('-o','--output', help='Output directory', required=True)
    parser.add_argument('--option', required=True, default='circref-and-reads')
    args = vars(parser.parse_args())
    
    if args['option'] == 'reference-only':
        print("Create only circular reference genome")
        create_genome(args['input'], args['genome'], args['output'])
    else:
        print("Create circular reference genome + generate reads")
        simulate_reads_circular(args['input'], args['genome'], args['output'])

# python simcircles.py -i circles.bed -g ../../reference/GRCh38.primary_assembly.genome.fa -o dir
# python3.7 simcircles.py -i /Users/madag/Projects/PhD/reconstruct/simulation/circles_new.bed -g /Users/madag/Projects/PhD/reference/GRCh38.primary_assembly.genome.fa -o dir
