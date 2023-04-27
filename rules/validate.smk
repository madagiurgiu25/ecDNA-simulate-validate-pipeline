rule assembly_qc:
    input:
        assembly="{outdir}/{sample}/assembly/Assembly.fasta",
        ref="{outdir}/{sample}/reads/circle-template.fa",
        gtf=GTF
    output:
        directory("{outdir}/{sample}/assembly_qc")
    threads: THREADS
    conda:
        "../envs/assembly.yaml"
    log:
        "{outdir}/{sample}/logs/assembly_qc.log"
    shell:
        """
        quast {input.assembly} \
            -r {input.ref} \
            -g {input.gtf} \
            -o {output} &> {log}
        """

rule assembly_against_ref:
    input:
        assembly="{outdir}/{sample}/assembly/Assembly.fasta",
        refzoom="{outdir}/{sample}/reads/circle-genomeref.fa",
        ref=REF,
        refmmi=REF_MMI,
        gtf=GTF
    output:
        quast=directory("{outdir}/{sample}/assembly_qc_ref"),
        sam="{outdir}/{sample}/assembly/assembly_bwa.sam",
        xsam="{outdir}/{sample}/assembly/assembly_minimap2_asm5.sam",
        sam_spliced="{outdir}/{sample}/assembly/assembly_minimap2_spliced.sam"
    threads: THREADS
    params:
        BWA_ALING
    resources: mem_mb=MAX_MEMORY
    conda:
        "../envs/assembly.yaml"
    log:
        quast="{outdir}/{sample}/logs/assembly_against_ref_quast.log"
    shell:
        """
        quast {input.assembly} \
            -r {input.refzoom} \
            -g {input.gtf} \
            -o {output.quast} &> {log.quast}
        bwa mem {params} -t {threads} {input.ref} {input.assembly} | samtools sort -o {output.sam} -
        minimap2 -t 1 -a -x asm5 {input.refmmi} {input.assembly}  | samtools sort -o {output.xsam} -  
        minimap2 -t 1 -a -x splice {input.refmmi} {input.assembly} | samtools sort -o {output.sam_spliced} - 
        """

rule compare_reconstruc2denovo:
    input:
        shasta="{outdir}/{sample}/assembly/Assembly.fasta",
        unicycler="{outdir}/{sample}/assembly_unicycler/assembly.fasta",
	    ref="{outdir}/{sample}/reads/circle-template.fa"
    output:
        blast_shasta="{outdir}/{sample}/validate/blast_shasta.txt",
        blast_unicycler="{outdir}/{sample}/validate/blast_unicycler.txt"
    conda:
        "../envs/validate.yaml"
    shell:
        """
        echo -e "#qseqid\tsseqid\tpident\tlength\tevalue" > {output.blast_shasta} && \
        blastn \
            -query {input.shasta} \
            -subject {input.ref} \
            -task 'megablast' -outfmt 6 >> {output.blast_shasta}
            
        echo -e "#qseqid\tsseqid\tpident\tlength\tevalue" > {output.blast_unicycler} && \
        blastn \
            -query {input.unicycler} \
            -subject {input.ref} \
            -task 'megablast' -outfmt 6 >> {output.blast_unicycler}
        """


rule compare_reconstruct2simulation:
    input:
        ain_old="{outdir}/{sample}/reconstruct_v1/reconstruct.fasta",
        ain_new="{outdir}/{sample}/reconstruct_v2/reconstruct.fasta",
        # shasta="{outdir}/{sample}/assembly/Assembly.fasta",
        # unicycler="{outdir}/{sample}/assembly_unicycler/assembly.fasta",
	    ref="{outdir}/{sample}/reads/circle-template.fa"
    output:
        # strecher="{outdir}/{sample}/reconstruct/similarity_strecher.txt",
        blast1="{outdir}/{sample}/validate/blast_decoil_v1.txt",
        blast2="{outdir}/{sample}/validate/blast_decoil_v2.txt"
        # blast_shasta="{outdir}/{sample}/validate/blast_shasta.txt",
        # blast_unicycler="{outdir}/{sample}/validate/blast_unicycler.txt"
    conda:
        "../envs/validate.yaml"
    # params:
    #     validatescript=VALIDATE_SCRIPT
    # log:
    #     "{outdir}/{sample}/logs/validation.log"
    shell:
        """
        echo -e "#qseqid\tsseqid\tpident\tlength\tevalue" > {output.blast1} && \
        blastn \
            -query {input.ain_old} \
            -subject {input.ref} \
            -task 'megablast' -outfmt 6 >> {output.blast1}
        
        echo -e "#qseqid\tsseqid\tpident\tlength\tevalue" > {output.blast2} && \
        blastn \
            -query {input.ain_new} \
            -subject {input.ref} \
            -task 'megablast' -outfmt 6 >> {output.blast2}
        """

        # echo -e "#qseqid\tsseqid\tpident\tlength\tevalue" > {output.blast}
        # && \
        # python3 {params.validatescript} \
        #     --input {input.ain} \
        #     --reference {input.ref} &> {log}

        # stretcher \
        # -asequence {input.ref} \
        # -bsequence {input.ain} \
        # -gapopen 10 -gapextend 1 \
        # -outfile {output.strecher}

rule overlap:
    input:
        rec1="{outdir}/{sample}/reconstruct_v1/reconstruct.fasta",
        rec2="{outdir}/{sample}/reconstruct_v2/reconstruct.fasta",
        circ="{outdir}/{sample}/reads/circle-template.fa",
        circ_split="{outdir}/{sample}/reads/circle-template-split.fa"
    output:
        ref2recons1="{outdir}/{sample}/validate/ref-to-reconstruct_v1.paf",
        ref2recons2="{outdir}/{sample}/validate/ref-to-reconstruct_v2.paf",
        ref2sim="{outdir}/{sample}/validate/ref-to-sim.paf"
    conda:
        "../envs/validate.yaml"
    shell:
        """
        minimap2 -x asm5 {input.circ_split} {input.rec1} > {output.ref2recons1}
        minimap2 -x asm5 {input.circ_split} {input.rec2} > {output.ref2recons2}
        minimap2 -x asm5 {input.circ_split} {input.circ} > {output.ref2sim}
        """
