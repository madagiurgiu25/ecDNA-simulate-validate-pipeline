import logging
import os
import json
import pandas as pd
import pprint

from snakemake import shell
from snakemake.utils import validate

report: "../reports/workflow.rst"

###### Config file and sample sheets #####
configfile: "configs/config.yaml"

CONFIG_GENOME = config["pickgenome"]
REF = config["genomes"][CONFIG_GENOME]["genome"]
REF_MMI = config["genomes"][CONFIG_GENOME]["mmi"]
GTF = config["genomes"][CONFIG_GENOME]["gtf"]
KMER = config["genomes"][CONFIG_GENOME]["kmer"]

OUTPUTDIR = config["outputdir"]

# mapper
MAPPER = config["mapper"]

# running mode of the pipeline
RUNMODE = config["runmode"]
RUNFULL = 'simulate-mapping-sv-reconstruct-validate'
RUNRECONSTRUCT = 'reconstruct'
RUNMAPPING = 'mapping-sv-reconstruct'
RUNSV = 'sv-reconstruct'
RUNCNV = 'cnv'
RUNSIM = 'simulate'
RUNSIMTOSV = 'simulate-mapping-sv'

run_options = [RUNFULL, RUNRECONSTRUCT, RUNMAPPING, RUNSV, RUNSIM, RUNCNV, RUNSIMTOSV]

if RUNMODE not in run_options:
    s = ",".join(run_options)
    raise Exception("""ERROR: invalid runmode in snakemake pipeline; possible definition runmode={}""".format(s))

print(RUNMODE)

# run mapping
if RUNMODE in [RUNMAPPING]:
	FASTQ = config["fastq"]

# run sv
if RUNMODE in [RUNSV]:
	BW = config["bw"]
	BAM = config["bam"]

# run reconstruct
if RUNMODE in [RUNRECONSTRUCT]:
	VCF = config["vcf"]
	BW = config["bw"]
	BAM = config["bam"]

SAMPLE = []
INPUT = {}
# input bed of the circular structure
if config["batch"] is None or len(config["batch"]) == 0:
	# single ecDNA template
	print("###DEBUG: Pipeline runs in single mode")
	SAMPLE = [config["sample"]]
	if RUNMODE in [RUNFULL, RUNSIMTOSV, RUNSIM]:
		INPUT = {config["sample"]:config["input"]}
		INPUT_TYPE = config["inputtype"]
	elif RUNMODE in [RUNSV]:
		INPUT = {config["sample"]: config["bam"]}
		INPUT_TYPE = "bam"
	elif RUNMODE in [RUNMAPPING]:
		INPUT = {config["sample"]: config["fastq"]}
		INPUT_TYPE = "fastq"
	elif RUNMODE in [RUNRECONSTRUCT]:
		INPUT = {config["sample"]: config["vcf"]}
		INPUT_TYPE = "vcf"

	# save config
	config_filename = os.path.join(OUTPUTDIR,SAMPLE[0],"config.yaml")
	os.makedirs(os.path.join(OUTPUTDIR, SAMPLE[0]), exist_ok=True)
	with open(config_filename, 'w') as f:
	    json.dump(config, f, indent=4)
else:
	SAMPLE = []
	INPUT = {}
	# batch of ecDNA templates
	print("###DEBUG: Pipeline runs in batch mode")
	INPUT_TYPE = config["inputtype"]
	ROOT = config["batch"]

	for line in shell("find {} -name '*.bed' ".format(ROOT), iterable=True):
		fin = line
		fname = line.split("/")[-1].split(".")[0]
		SAMPLE.append(fname)
		INPUT[fname] = fin

		# create config per sample
		config_filename = os.path.join(OUTPUTDIR,fname,"config.yaml")
		os.makedirs(os.path.join(OUTPUTDIR,fname), exist_ok=True)
		with open(config_filename, 'w') as f:
		    json.dump(config, f, indent=4)

def get_sim_depth(id, s):
	"""
	Get simulation depth/coverage for each circular structure in the sample
	"""
	slice = CONF_CIRCLES[s][id == CONF_CIRCLES[s].structure]
	cov = slice.coverage.tolist()[0]
	return cov

CONF_CIRCLES = {}
CIRCLES = {}
COVS = {}

if INPUT and len(INPUT) > 10:
	for s in INPUT:
		CONF_CIRCLES[s] = pd.read_csv(INPUT[s], header=0, sep="\t")
		CIRCLES[s] = CONF_CIRCLES[s].structure.drop_duplicates().tolist()

		# coverage per sample
		COVS[s] = {}
		for c in CIRCLES[s]:
			COVS[s][c] = get_sim_depth(c, s)

# pprint.pprint(CIRCLES)
# pprint.pprint(COVS)
# print(INPUT,SAMPLE,OUTPUTDIR)

# simulation params
if RUNMODE in [RUNFULL, RUNSIM, RUNSIMTOSV]:
	DEPTH = config["simulation"]["depth"]
	AVGLEN = config["simulation"]["avglen"]
	ACC_MEAN = config["simulation"]["macc"]  # mean-accuracy of the model (simulate noisy reads)
	HMM_MODEL = config["simulation"]["hmm"]
	CIRCULAR = config["simulation"]["circular"]

# assembly params
if RUNMODE in [RUNFULL]:
	MINLEN = config["filtering"]["assembly"]["readlen"]
	BWA_ALING = config["params"]["bwaasm"]
	# assembly align to reference
	MAX_MEMORY = config["params"]["maxmem"]

# mapping params
READ_MINLEN = config["filtering"]["mapping"]["readlen"]
READ_QUAL = config["filtering"]["mapping"]["qual"]

# cnv params
CHRSIZE = config["cnv"]["hg38"]["chromsize"]
BINS = config["cnv"]["hg38"]["bins"]
EXCLUDE = config["cnv"]["hg38"]["exclude"]
GCCONTENT = config["cnv"]["hg38"]["gccontent"]
SCRIPT_SMURFSEQ= config["cnv"]["smurfseq"]


# reconstruct
if RUNMODE in [RUNFULL, RUNRECONSTRUCT, RUNSV]:
	DECOIL_SCRIPT = config["decoil"]

# validate assembly
if RUNMODE in [RUNFULL]:
	VALIDATE_SCRIPT = config["validate"]

# other parameters
THREADS = config["params"]["threads"]

# prevent snakemake of changing wildcards
wildcard_constraints:
    sample="|".join(SAMPLE),
    outdir="|".join([OUTPUTDIR])

