def get_mem_gb(wildcards, attempt):
    return attempt * 3000

rule ngmlr:
    input:
        fastq="{outdir}/{sample}/filt.fastq",
#        fastq=lambda wildcards: "{outdir}/{sample}/reads/circle-template.fastq" if RUNMODE in [RUNFULL, RUNMAPPING, RUNSIMTOSV]  else INPUT,
        ref=REF
    output:
        bam="{outdir}/{sample}/mapping/map.bam",
        bai="{outdir}/{sample}/mapping/map.bam.bai",
        sam=temp("{outdir}/{sample}/mapping/map.sam")
    conda:
        "../envs/mapping.yaml"
    params:
        threads=THREADS
    resources:
        mem_mb=get_mem_gb
    shell:
        """
        ngmlr --bam-fix --threads {params.threads} --reference {input.ref} --query {input.fastq} --output {output.sam} --presets ont
        samtools sort -O BAM -o {output.bam} {output.sam}
        samtools index {output.bam}
        """
