# split config file into
# individual circles

def get_input(wildcards):
	return INPUT[wildcards.sample]

def get_ids_circles(wildscards):
	l = CONF_CIRCLES[wildcards.sample].structure.drop_duplicates().tolist()
	return l

def get_ids(s):
	l = CONF_CIRCLES[s].structure.drop_duplicates().tolist()
	return l

def get_merge_fastq_files(wildcards):
	ids = get_ids(wildcards.sample)
	files = expand("{{outdir}}/{{sample}}/reads/circle-template-{id}.fastq", id=ids)
	# print(files)
	return files

def get_merge_fasta_files(wildcards):
	ids = get_ids(wildcards.sample)
	files = expand("{{outdir}}/{{sample}}/reads/circle-template-{id}.fa", id=ids)
	# print(files)
	return files

rule split_bed:
	input:
		bedfile = get_input
	output:
		"{outdir}/{sample}/reads/circle-template-{id}.bed"
	conda:
		"../envs/simulate.yaml"
	params:
		id=lambda wildcards: CONF_CIRCLES[wildcards.sample].structure.drop_duplicates().tolist()
	shell:
		"""
		# use grep -P in linux but for macos is grep -Ei
		(head -1 {input.bedfile}; grep -P "\t{wildcards.id}\t" {input}) > {output}
		"""

rule get_circle_fasta:
	input:
		bed="{outdir}/{sample}/reads/circle-template-{id}.bed",
		ref=REF
	output:
		ofa="{outdir}/{sample}/reads/circle-template-{id}.fa"
	conda:
		"../envs/simulate.yaml"
	log:
		"{outdir}/{sample}/logs/circle_fasta-{id}.log"
	shell:
		"""
		mkdir -p {wildcards.outdir}
		python -u scripts/simcircles.py \
			--input {input.bed} \
			--output {output.ofa} \
			--genome {input.ref} \
			--option reference-only &> {log}
		"""

rule sim_fastq:
	input:
		fa="{outdir}/{sample}/reads/circle-template-{id}.fa",
		bed="{outdir}/{sample}/reads/circle-template-{id}.bed"
	output:
		out="{outdir}/{sample}/reads/circle-template-{id}.fastq",
		dir=temp(directory("{outdir}/tmp_{sample}_{id}"))
	shadow: "shallow"
	threads: THREADS
	params:
		id="{id}",
		depth=lambda wildcards: COVS[wildcards.sample]["{}".format(wildcards.id)],
		avglen=AVGLEN,
		macc=ACC_MEAN,
		hmm=HMM_MODEL,
		circular=CIRCULAR
	log:
		"{outdir}/{sample}/logs/sim_fastq-{id}.log"
	shell:
		"""
		mkdir -p {output.dir}
		mkdir -p {wildcards.outdir}/logs-{wildcards.id}
		pbsim --depth {params.depth} --length-mean {params.avglen} --seed 500 \
			--prefix {output.dir}/circ \
			--accuracy-mean {params.macc} \
			--hmm_model {params.hmm} \
			--circular {params.circular} \
			{input.fa} &> {log}
		cat {output.dir}/circ*.fastq | awk -v id={wildcards.id} '{{if (NR % 4 == 1 || NR % 4 == 3) {{print $0"_"id}} else {{print}} }}' > {output.out}
		"""


rule merge_fastq:
	input:
		get_merge_fastq_files
		# expand("{{outdir}}/{{sample}}/reads/circle-template-{id}.fastq",id=get_ids("{{sample}}"))
	output:
		"{outdir}/{sample}/reads/circle-template.fastq",
	shell:
		"""
		cat {input} > {output}
		"""

rule merge_fasta:
	input:
		get_merge_fasta_files
		# expand("{{outdir}}/{{sample}}/reads/circle-template-{id}.fa",id=get_ids("{{sample}}"))
	output:
		"{outdir}/{sample}/reads/circle-template.fa",
	shell:
		"""
		cat {input} > {output}
		"""

rule merge_bed:
	input:
		get_input
	output:
		"{outdir}/{sample}/reads/circle-template.bed",
	shell:
		"""
		cp {input} {output}
		"""

# get all reference overlapping regions
rule merge_regions:
	input:
		bed=get_input
	output:
		"{outdir}/{sample}/reads/circle-genomeref.bed"
	conda:
		"../envs/simulate.yaml"
	shell:
		"""
		head -1 {input.bed} > {output}
		tail -n+2 {input.bed} | sort -k1,1 -k2,2n | bedtools merge -i -  >> {output}
		"""

# get fasta regions
rule fasta_regions:
	input:
		bed="{outdir}/{sample}/reads/circle-genomeref.bed",
		ref=REF
	output:
		"{outdir}/{sample}/reads/circle-genomeref.fa",
	conda:
		"../envs/simulate.yaml"
	shell:
		"""
		bedtools getfasta -fi {input.ref} -bed {input.bed} > {output}
		"""

rule fasta_split:
	input:
		bed="{outdir}/{sample}/reads/circle-template.bed",
		ref=REF
	output:
		"{outdir}/{sample}/reads/circle-templatesplit.fa",
	conda:
		"../envs/simulate.yaml"
	shell:
		"""
		bedtools getfasta -fi {input.ref} -bed {input.bed} > {output}
		"""

rule fasta_index:
	input:
		fa="{outdir}/{sample}/reads/circle-template.fa",
		gfa="{outdir}/{sample}/reads/circle-genomeref.fa",
		fasplit="{outdir}/{sample}/reads/circle-templatesplit.fa"
	output:
		faidx="{outdir}/{sample}/reads/circle-template.fa.fai",
		gfaidx="{outdir}/{sample}/reads/circle-genomeref.fa.fai",
		fasplitidx="{outdir}/{sample}/reads/circle-templatesplit.fa.fai"
	conda:
		"../envs/simulate.yaml"
	shell:
		"""
		samtools faidx {input.fa}
		samtools faidx {input.gfa}
		samtools faidx {input.fasplit}
		"""




	# # generate the circle fasta
	# rule compose_circle_fasta:
	#     input:
	#         bed=INPUT,
	#         ref=REF
	#     output:
	#         ofa="{outdir}/{sample}/reads/circle-template.fa",
	#         ofb="{outdir}/{sample}/reads/circle-template.bed"
	#     conda:
	#         "../envs/simulate.yaml"
	#     log:
	#         "{outdir}/{sample}/logs/circle_fasta.log"
	#     shell:
	#         """
	#         mkdir -p {wildcards.outdir}/logs
	#         python -u scripts/simcircles.py \
	#             --input {input.bed} \
	#             --output {output.ofa} \
	#             --genome {input.ref} \
	#             --option reference-only &> {log}
	#         cp {input.bed} {output.ofb}
	#         """
	#
	# rule fasta_split:
	#     input:
	#         bed="{outdir}/{sample}/reads/circle-template.bed",
	#         ref=REF
	#     output:
	#         "{outdir}/{sample}/reads/circle-template-split.fa",
	#     conda:
	#         "../envs/simulate.yaml"
	#     shell:
	#         """
	#         bedtools getfasta -fi {input.ref} -bed {input.bed} > {output}
	#         """
	#
	# # get all reference overlapping regions
	# rule merge_regions:
	#     input:
	#         bed=INPUT
	#     output:
	#         "{outdir}/{sample}/reads/circle-genomeref.bed"
	#     conda:
	#         "../envs/simulate.yaml"
	#     shell:
	#         """
	#         head -1 {input.bed} > {output}
	#         tail -n+2 {input.bed} | sort -k1,1 -k2,2n | bedtools merge -i -  >> {output}
	#         """
	#
	# # get fasta regions
	# rule fasta_regions:
	#     input:
	#         bed="{outdir}/{sample}/reads/circle-genomeref.bed",
	#         ref=REF
	#     output:
	#         "{outdir}/{sample}/reads/circle-genomeref.fa",
	#     conda:
	#         "../envs/simulate.yaml"
	#     shell:
	#         """
	#         bedtools getfasta -fi {input.ref} -bed {input.bed} > {output}
	#         """
	#
	# rule fasta_index:
	#     input:
	#         fa="{outdir}/{sample}/reads/circle-template.fa",
	#         gfa="{outdir}/{sample}/reads/circle-genomeref.fa",
	#         fasplit="{outdir}/{sample}/reads/circle-template-split.fa"
	#     output:
	#         faidx="{outdir}/{sample}/reads/circle-template.fa.fai",
	#         gfaidx="{outdir}/{sample}/reads/circle-genomeref.fa.fai",
	#         fasplitidx="{outdir}/{sample}/reads/circle-template-split.fa.fai"
	#     shell:
	#         """
	#         samtools faidx {input.fa}
	#         samtools faidx {input.gfa}
	#         samtools faidx {input.fasplit}
	#         """
	#
	# # simulate fastq reads
	# rule sim_fastq:
	#     input:
	#         fa="{outdir}/{sample}/reads/circle-template.fa"
	#     output:
	#         out="{outdir}/{sample}/reads/circle-template.fastq",
	#         dir=temp(directory("{outdir}/tmp_{sample}"))
	#     shadow: "shallow"
	#     threads: THREADS
	#     params:
	#         depth=DEPTH,
	#         avglen=AVGLEN,
	#         macc=ACC_MEAN,
	#         hmm=HMM_MODEL,
	#         circular=CIRCULAR
	#     log:
	#         "{outdir}/{sample}/logs/sim_fastq.log"
	#     shell:
	#         """
	#         mkdir -p {output.dir}
	#         mkdir -p {wildcards.outdir}/logs
	#         pbsim --depth {params.depth} --length-mean {params.avglen} --seed 500 \
	#             --prefix {output.dir}/circ \
	#             --accuracy-mean {params.macc} \
	#             --hmm_model {params.hmm} \
	#             --circular {params.circular} \
	#             {input.fa} &> {log}
	#         cat {output.dir}/circ*.fastq > {output.out}
	#         """
