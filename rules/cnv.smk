"""
CNV calling pipeline.
"""

rule getbins:
    input:
        bam="{outdir}/{sample}/mapping/ngmlr.bam",
        chrsize=CHRSIZE,
        bins=BINS
    output:
        t=temp("{outdir}/{sample}/cnv/ngmlr.sam"),
        bincounts="{outdir}/{sample}/cnv/bin_counts.bed",
        binstats="{outdir}/{sample}/cnv/bin_stats.txt"
    params:
        smurfseq=SCRIPT_SMURFSEQ
    shell:
        """
        export PATH={params.smurfseq}:$PATH
        
        samtools view -h -O SAM {input.bam} > {output.t}
        getBinCounts.py -i {output.t} \
            -c {input.chrsize} \
            -b {input.bins} \
            -o {output.bincounts} \
            -s {output.binstats}
        """

rule cnv_calling:
    input:
        bincounts="{outdir}/{sample}/cnv/bin_counts.bed",
        gccontent=GCCONTENT,
        exclude=EXCLUDE
    output:
        "{outdir}/{sample}/cnv/{sample}.pdf",
        "{outdir}/{sample}/cnv/{sample}.data.txt"
    log:
        "{outdir}/{sample}/logs/cnv.log"
    params:
        dir="{outdir}/{sample}/cnv"
    shell:
        """
        cnvAnalysis.R {input.bincounts} {params.dir}/{wildcards.sample} {input.gccontent} {input.exclude} &> {log}
        """

rule get_levels:
    input:
        "{outdir}/{sample}/cnv/{sample}.data.txt"
    output:
        "{outdir}/{sample}/cnv/{sample}.levels.txt"
    log:
        "{outdir}/{sample}/logs/cnv_levels.log"
    shell:
        """
        touch {output}
        # python scripts/cnvconvert.py --fin {input} --fout {output} &> {log}
        """
