import argparse
from Bio import pairwise2
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def compare_assembly2reference(infasta, reffasta):

	records_reference = SeqIO.to_dict(SeqIO.parse(open(reffasta), 'fasta'))
	records_asm = SeqIO.to_dict(SeqIO.parse(open(infasta), 'fasta'))

	for key1 in records_reference:
		for key2 in records_asm:
			list_scores = []
			seq1 = records_reference[key1].seq

			# determine shortest sequence length
			min_len = min(len(records_reference[key1].seq), len(records_asm[key2].seq))
			print(min_len)

			# compute similarity over the circular sequence
			for offset in range(0, len(records_asm[key2].seq), 1000):
				seq2 = records_asm[key2].seq[offset:] + records_asm[key2].seq[:offset]
				alignment_score = pairwise2.align.globalms(seq1, seq2, 2, -1, -.5, -.1, score_only=True, penalize_end_gaps=False, penalize_extend_when_opening=True)
				list_scores.append(alignment_score)

			print(key1, key2, max(list_scores)/min_len)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Compare circular genomes')
	parser.add_argument('-i', '--input', help='Assembly fasta', required=True)
	parser.add_argument('-r', '--reference', help='Reference genome (fasta)', required=True)
	args = parser.parse_args()

	compare_assembly2reference(args.input, args.reference)
