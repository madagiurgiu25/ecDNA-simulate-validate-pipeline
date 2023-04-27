rule assembly_shasta:
    input:
        fastq=lambda wildcards: "{outdir}/{sample}/reads/circle-template.fastq" if RUNMODE in [RUNFULL,RUNMAPPING] else INPUT
    output:
        dir=directory("{outdir}/{sample}/assembly"),
        f="{outdir}/{sample}/assembly/Assembly.fasta"
    threads: THREADS
    params:
        minlen=MINLEN
    conda:
        "../envs/assembly.yaml"
    log:
        "{outdir}/{sample}/logs/assembly.log"
    shell:
        """
        echo {input}
        mkdir -p {output.dir}
        shasta \
            --threads {threads}\
            --input {input.fastq} \
            --config Nanopore-Oct2021 \
            --Reads.minReadLength {params.minlen} \
            --assemblyDirectory {output.dir} &> {log}
        """

rule assembly_unicycler:
    input:
        fastq=lambda wildcards: "{outdir}/{sample}/reads/circle-template.fastq" if RUNMODE in [RUNFULL,RUNMAPPING] else INPUT
    output:
        asm="{outdir}/{sample}/assembly_unicycler/assembly.fasta"
    conda:
        "../envs/unicycler.yaml"
    params:
        dir="{outdir}/{sample}/assembly_unicycler"
    log:
        "{outdir}/logs/{sample}/assembly_unicycler.log"
    shell:
        """
        unicycler \
            -l {input} \
            --min_fasta_length 1000 -t 2 --mode normal --keep 2 \
            --out {params.dir} &> {log}
        """

# cat {output.sam_spliced} | awk -F"\t" '{{for(i=1;i<=NF;i++) if (i==21) {{printf "zd:i:1\t"$i"\t"}} else {{printf $i"\t"}}; print ""}}' > {output.sam_spliced_zd}
