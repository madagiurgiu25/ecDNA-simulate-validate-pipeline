# stucture from https://github.com/alperyilmaz/conda-snakemake-gatk
from snakemake.utils import min_version
min_version("6.0")

import sys
print("Python version")
print(sys.executable)

# intialize all parameters
include: "rules/initialize.smk"

if RUNMODE in [RUNFULL, RUNSIM, RUNSIMTOSV]:
	include: "rules/simulate.smk"

if RUNMODE in [RUNFULL]:
	include: "rules/validate.smk"

if RUNMODE in [RUNFULL, RUNMAPPING, RUNSIMTOSV]:
	if MAPPER == "ngmlr":
		include: "rules/mapping_ngmlr.smk"
	else:
		include: "rules/mapping_winnowmap.smk"
	include: "rules/filter.smk"
	include: "rules/cov.smk"

if RUNMODE in [RUNFULL]:
	include: "rules/assembly.smk"

if RUNMODE in [RUNFULL, RUNMAPPING, RUNSV, RUNSIMTOSV]:
	include: "rules/sv.smk"

if RUNMODE in [RUNFULL, RUNCNV]:
	include: "rules/cnv.smk"

if RUNMODE in [RUNFULL, RUNRECONSTRUCT, RUNSV, RUNMAPPING]:
	include: "rules/reconstruct.smk"

# include: "rules/circlator.smk"

def get_files_simulation(run_mode):
	if run_mode in [RUNFULL, RUNSIM, RUNSIMTOSV]:
		return expand([
			"{outdir}/{sample}/reads/circle-template.bed",
			"{outdir}/{sample}/reads/circle-template.fa",
			"{outdir}/{sample}/reads/circle-template.fa.fai",
			"{outdir}/{sample}/reads/circle-template.fastq",
			"{outdir}/{sample}/reads/circle-genomeref.bed",
			"{outdir}/{sample}/reads/circle-genomeref.fa",
			"{outdir}/{sample}/reads/circle-genomeref.fa.fai",
			"{outdir}/{sample}/reads/circle-templatesplit.fa",
			"{outdir}/{sample}/reads/circle-templatesplit.fa.fai",
		],outdir=OUTPUTDIR,sample=SAMPLE,id=CIRCLES)

	return []

def get_temp_files_simulation(run_mode):
	output=[]
	if run_mode in [RUNFULL, RUNSIM, RUNSIMTOSV]:
		for s in SAMPLE:
			for i in CIRCLES[s]:
				istr = str(i)
				output.append("""{}/{}/reads/circle-template-{}.fa""".format(OUTPUTDIR,s,istr))
				output.append("""{}/{}/reads/circle-template-{}.bed""".format(OUTPUTDIR,s,istr))
				output.append("""{}/{}/reads/circle-template-{}.fastq""".format(OUTPUTDIR,s,istr))
				output.append("""{}/{}/logs/sim_fastq-{}.log""".format(OUTPUTDIR,s,istr))
				output.append("""{}/{}/logs/circle_fasta-{}.log""".format(OUTPUTDIR,s,istr))
		return output
	return []


def get_file_mapping(run_mode):
	if run_mode in [RUNFULL, RUNMAPPING, RUNSIMTOSV]:
		return expand([
			"{outdir}/{sample}/mapping/coverage.bw",
			"{outdir}/{sample}/mapping/ngmlr.bam",
			"{outdir}/{sample}/mapping/ngmlr.bam.bai",
			"{outdir}/{sample}/logs/cov.log"
		],outdir=OUTPUTDIR,sample=SAMPLE)
	return []


def get_file_assembly(run_mode):
	if run_mode in [RUNFULL]:
		return expand(["{outdir}/{sample}/assembly/Assembly.fasta",
		               "{outdir}/{sample}/assembly_unicycler/assembly.fasta",
		               "{outdir}/{sample}/logs/assembly.log"
		               ],outdir=OUTPUTDIR,sample=SAMPLE)
	return []


def get_file_sv(run_mode):
	if run_mode in [RUNFULL, RUNMAPPING, RUNSV, RUNSIMTOSV]:
		return expand([
			"{outdir}/{sample}/sv/sv.sniffles.vcf",
			"{outdir}/{sample}/sv/sv.sniffles.bedpe",
			"{outdir}/{sample}/logs/survivor.log"
		],outdir=OUTPUTDIR,sample=SAMPLE)

	return []


def get_files_reconstruction(run_mode):
	if run_mode in [RUNFULL, RUNRECONSTRUCT, RUNMAPPING, RUNSV]:
		return expand([
			"{outdir}/{sample}/sv/sv.sniffles.filtered.vcf",
			"{outdir}/{sample}/sv/sv.sniffles.filtered.bedpe",
            "{outdir}/{sample}/reconstruct/reconstruct.ecDNA.fasta",
            "{outdir}/{sample}/reconstruct/reconstruct.ecDNA.bed",
#			"{outdir}/{sample}/reconstruct_v1/reconstruct.fasta",
#			"{outdir}/{sample}/reconstruct_v1/reconstruct.bed",
#			"{outdir}/{sample}/reconstruct_v2/reconstruct.fasta",
#			"{outdir}/{sample}/reconstruct_v2/reconstruct.bed",
# 			"{outdir}/{sample}/reconstruct_v3/reconstruct.fasta",
# 			"{outdir}/{sample}/reconstruct_v3/reconstruct.bed",
			"{outdir}/{sample}/logs/filter.log",
			"{outdir}/{sample}/logs/reconstruct_v1.log",
#			"{outdir}/{sample}/logs/reconstruct_v2.log",
			# "{outdir}/{sample}/logs/reconstruct_v3.log",
			"{outdir}/{sample}/logs/survivor_filtered.log"
		],outdir=OUTPUTDIR,sample=SAMPLE)
	return []


def get_files_validation(run_mode):
	if run_mode == RUNFULL:
		return expand(["{outdir}/{sample}/assembly_qc",
		               "{outdir}/{sample}/assembly_qc_ref",
		               "{outdir}/{sample}/assembly/assembly_bwa.sam",
		               "{outdir}/{sample}/assembly/assembly_minimap2_asm5.sam",
		               "{outdir}/{sample}/assembly/assembly_minimap2_spliced.sam",
		               "{outdir}/{sample}/validate/blast_decoil_v1.txt",
		               "{outdir}/{sample}/validate/blast_decoil_v2.txt",
		               "{outdir}/{sample}/validate/blast_shasta.txt",
		               "{outdir}/{sample}/validate/blast_unicycler.txt",
		               "{outdir}/{sample}/logs/assembly_against_ref_quast.log"],outdir=OUTPUTDIR,sample=SAMPLE)

	return []


rule all:
	input:
		config=expand("{outdir}/{sample}/config.yaml", outdir=OUTPUTDIR, sample=SAMPLE),
		simulate=get_files_simulation(RUNMODE),
		simulate_temp=get_temp_files_simulation(RUNMODE),
		assembly=get_file_assembly(RUNMODE),
		mapping=get_file_mapping(RUNMODE),
		sv=get_file_sv(RUNMODE),
		reconstruct=get_files_reconstruction(RUNMODE),
		validation=get_files_validation(RUNMODE)

# snakemake --cores 32 -n --use-conda --configfile configs/config.yaml --config batch=/data/gpfs-1/users/giurgium_c/Nanopore/Simulation/sim_all/batch_test outputdir=/data/gpfs-1/users/giurgium_c/Nanopore/Simulation/process/batch_test
