configfile: "configs/ichorCNA.yaml"
configfile: "configs/config.yaml"

INPUT = config["input"]
SAMPLE = config["sample"]
OUTPUTDIR = config["outputdir"]

rule all:
	input:
		expand("{outdir}/{sample}/cnv/ichorCNA/all.cna.seg",outdir=OUTPUTDIR,sample=SAMPLE),
		expand("{outdir}/{sample}/cnv/ichorCNA/bin.wig",outdir=OUTPUTDIR,sample=SAMPLE)

rule read_counter:
	input:
		INPUT
	output:
		"{outdir}/{sample}/cnv/ichorCNA/bin.wig"
	params:
		binSize=config["binSize"],
		qual="20",
		chrs=config["chrs"]
	resources:
		mem=4
	conda:
		"../../envs/ichorCNA.yaml"
	log:
		"{outdir}/{sample}/logs/ichorDNA_bin.log"
	shell:
		"""
		# install HMMcopy to use readCounter
		readCounter {input} -c {params.chrs} -w {params.binSize} -q {params.qual} > {output} 2> {log}
		"""

rule ichorCNA:
	input:
		"{outdir}/{sample}/cnv/ichorCNA/bin.wig"
	output:
		cna="{outdir}/{sample}/cnv/ichorCNA/all.cna.seg"
	params:
		outDir="{outdir}/{sample}/cnv/ichorCNA",
		rscript=config["ichorCNA_rscript"],
		id="{sample}",
		ploidy=config["ichorCNA_ploidy"],
		normal=config["ichorCNA_normal"],
		gcwig=config["ichorCNA_gcWig"],
		mapwig=config["ichorCNA_mapWig"],
		normalpanel=config["ichorCNA_normalPanel"],
		estimateNormal=config["ichorCNA_estimateNormal"],
		estimatePloidy=config["ichorCNA_estimatePloidy"],
		estimateClonality=config["ichorCNA_estimateClonality"],
		scStates=config["ichorCNA_scStates"],
		maxCN=config["ichorCNA_maxCN"],
		includeHOMD=config["ichorCNA_includeHOMD"],
		chrs=config["ichorCNA_chrs"],
		chrTrain=config["ichorCNA_chrTrain"],
		genomeBuild=config["ichorCNA_genomeBuild"],
		genomeStyle=config["ichorCNA_genomeStyle"],
		centromere=config["ichorCNA_centromere"],
		fracReadsChrYMale=config["ichorCNA_fracReadsInChrYForMale"],
		minMapScore=config["ichorCNA_minMapScore"],
		maxFracGenomeSubclone=config["ichorCNA_maxFracGenomeSubclone"],
		maxFracCNASubclone=config["ichorCNA_maxFracCNASubclone"],
		exons=config["ichorCNA_exons"],
		txnE=config["ichorCNA_txnE"],
		txnStrength=config["ichorCNA_txnStrength"],
		plotFileType=config["ichorCNA_plotFileType"],
		plotYlim=config["ichorCNA_plotYlim"],
		libdir=config["ichorCNA_libdir"]
	resources:
		mem=4
	conda:
		"../../envs/ichorCNA.yaml"
	log:
		"{outdir}/{sample}/logs/ichorCNA_calling.log"
	shell:
		"Rscript {params.rscript} \
		--id {params.id} \
		--WIG {input} \
		--gcWig {params.gcwig} \
		--mapWig {params.mapwig} \
		--genomeStyle {params.genomeStyle} \
		--genomeBuild {params.genomeBuild} \
		--outDir {params.outDir} > {log} 2> {log}"
