
rule nanofilt:
    input:
        # choose as input the fastq simulated or real fastq (given as param)
        lambda wildcards: "{outdir}/{sample}/reads/circle-template.fastq" if RUNMODE in [RUNFULL, RUNMAPPING, RUNSIMTOSV]  else FASTQ
    output:
        temp("{outdir}/{sample}/filt.fastq")
    params:
        rlen = READ_MINLEN,
        qual = READ_QUAL
    conda:
        "../envs/mapping.yaml"
    shell:
        "NanoFilt -l {params.rlen} -q {params.qual} --headcrop 20 --tailcrop 20 --readtype 1D {input} > {output}"
