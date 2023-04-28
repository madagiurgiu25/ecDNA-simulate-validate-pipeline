rule filter:
	input:
		# choose as input the sv simulated calles or real vcf (given as param)
		lambda wildcards: "{outdir}/{sample}/sv/sv.sniffles.vcf" if RUNMODE in [RUNFULL, RUNMAPPING, RUNSV, RUNSIMPLE] else VCF
	output:
		"{outdir}/{sample}/sv/sv.sniffles.filtered.vcf"
#	conda:
#		"../envs/reconstruct.yaml"
	params:
		threads=THREADS,
		name="{sample}"
	log:
		"{outdir}/{sample}/logs/filter.log"
	shell:
		"""
		# cd ../decoil/ && python -m pip install -e . && cd ../simulation && \
		mkdir -p "{wildcards.outdir}/{wildcards.sample}/sv" && \
		decoil filter -i {input} -o {output}
		"""

# cat {input} | awk -F"\t" '{{if($0 ~ /#/) {{print}} else {{ if ($7 ~ /(PASS|STRANDBIAS)/) {{print}} }} }}' > {output}


rule survivor_filtered:
	priority: 1
	input:
		filt="{outdir}/{sample}/sv/sv.sniffles.filtered.vcf",
	output:
		filt="{outdir}/{sample}/sv/sv.sniffles.filtered.bedpe",
	conda:
		"../envs/sv.yaml"
	log:
		"{outdir}/{sample}/logs/survivor_filtered.log"
	shell:
		"""
		SURVIVOR vcftobed {input.filt} -1 -1 {output.filt} &> {log}
		"""

rule reconstruct_v1:
	input:
		# vcf="{outdir}/{sample}/sv/sv.sniffles.filtered.vcf",
		vcf="{outdir}/{sample}/sv/sv.sniffles.vcf",
		ref=REF,
		bigwig=lambda wildcards: "{outdir}/{sample}/mapping/coverage.bw" if RUNMODE in [RUNFULL, RUNMAPPING, RUNSV, RUNSIMPLE] else BW,
		bam=lambda wildcards: "{outdir}/{sample}/mapping/map.bam" if RUNMODE in [RUNFULL, RUNMAPPING, RUNSV, RUNSIMPLE] else BAM
	output:
		fasta="{outdir}/{sample}/reconstruct/reconstruct.ecDNA.fasta",
		bed="{outdir}/{sample}/reconstruct/reconstruct.ecDNA.bed"
#	conda:
#		"../envs/reconstruct.yaml"
	params:
		outdir="{outdir}/{sample}/reconstruct",
		name="{sample}"
	log:
		"{outdir}/{sample}/logs/reconstruct_v1.log"
	shell:
		"""
		decoil reconstruct \
			--min-vaf 0.3 --min-cov-alt 8 --min-cov 10 \
			--fragment-min-cov 10 --fragment-min-size 500 \
			--vcf {input.vcf} \
			--genome {input.ref} \
			--bam {input.bam} \
			--coverage {input.bigwig} \
			--name {params.name} \
			--output {params.outdir} &> {log}
		"""

#rule reconstruct_v2:
#	input:
#		# vcf="{outdir}/{sample}/sv/sv.sniffles.filtered.vcf",
#		vcf="{outdir}/{sample}/sv/sv.sniffles.vcf",
#		ref=REF,
#		bigwig=lambda wildcards: "{outdir}/{sample}/mapping/coverage.bw" if RUNMODE in [RUNFULL, RUNMAPPING, RUNSV] else BW,
#		bam=lambda wildcards: "{outdir}/{sample}/mapping/ngmlr.bam" if RUNMODE in [RUNFULL, RUNMAPPING, RUNSV] else BAM,
#	output:
#		fasta="{outdir}/{sample}/reconstruct_v2/reconstruct.fasta",
#		bed="{outdir}/{sample}/reconstruct_v2/reconstruct.bed"
#	conda:
#		"../envs/reconstruct.yaml"
#	params:
#		outdir="{outdir}/{sample}/reconstruct_v2",
#		decoilscript=DECOIL_SCRIPT,
#		name="{sample}"
#	log:
#		"{outdir}/{sample}/logs/reconstruct_v2.log"
#	shell:
#		"""
#		python3 {params.decoilscript} --v2 \
#			--vcf {input.vcf} \
#			--genome {input.ref} \
#			--bam {input.bam} \
#			--coverage {input.bigwig} \
#			--name {params.name} \
#			--output {params.outdir} &> {log}
#		"""

# rule reconstruct_v3:
# 	input:
# 		# vcf="{outdir}/{sample}/sv/sv.sniffles.filtered.vcf",
# 		vcf="{outdir}/{sample}/sv/sv.sniffles.vcf",
# 		ref=REF,
# 		bigwig=lambda wildcards: "{outdir}/{sample}/mapping/coverage.bw" if RUNMODE in [RUNFULL, RUNMAPPING, RUNSV] else BW
# 	output:
# 		fasta="{outdir}/{sample}/reconstruct_v3/reconstruct.fasta",
# 		bed="{outdir}/{sample}/reconstruct_v3/reconstruct.bed"
# 	conda:
# 		"../envs/reconstruct.yaml"
# 	params:
# 		outdir="{outdir}/{sample}/reconstruct_v3",
# 		decoilscript=DECOIL_SCRIPT,
# 		name="{sample}"
# 	log:
# 		"{outdir}/{sample}/logs/reconstruct_v3.log"
# 	shell:
# 		"""
# 		python3 {params.decoilscript} -2 \
# 			--vcf {input.vcf} \
# 			--genome {input.ref} \
# 			--coverage {input.bigwig} \
# 			--name {params.name} \
# 			--output {params.outdir} &> {log}
# 		"""
