def get_mem_gb(wildcards, attempt):
    return attempt * 3000

rule winnowmap:
    input:
        fastq="{outdir}/{sample}/filt.fastq",
#        fastq=lambda wildcards: "{outdir}/{sample}/reads/circle-template.fastq" if RUNMODE in [RUNFULL, RUNMAPPING, RUNSIMTOSV]  else INPUT,
        ref=REF,
        kmers=KMERS
    output:
        bam="{outdir}/{sample}/mapping/map.bam",
        bai="{outdir}/{sample}/mapping/map.bam.bai",
        sam=temp("{outdir}/{sample}/mapping/map.sam")
    conda:
        "../envs/winnowmap.yaml"
    params:
        threads=THREADS
    resources:
        mem_mb=get_mem_gb
    shell:
        """
        winnowmap -W {input.kmers} -ax map-ont {input.ref} {input.fastq} > {output.sam}
        samtools sort -O BAM -o {output.bam} {output.sam}
        samtools index {output.bam}
        """

