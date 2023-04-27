# rule canu:
#     input:
#         "{outdir}/{sample}/reads/circle-template.fastq"
#     output:
#         output="{outdir}/{sample}/assembly_canu/assembly.fasta"
#     conda:
#         "../envs/circlator.yaml"
#     params:
#         dir="{outdir}/{sample}/assembly_canu"
#     log:
#         "{outdir}/{sample}/logs/canu.log"
#     shell:
#         """
#         canu -p assembly \
#             -d {params.dir} \
#             genomesize=10m maxThreads=4 useGrid=false maxMemory=4g obtovlThreads=4 \
#             -nanopore-raw {input} &> {log}
#         """
rule flye:
    input:
        "{outdir}/{sample}/reads/circle-template.fastq"
    output:
        "{outdir}/{sample}/assembly_flye/assembly.fasta"
    params:
        dir="{outdir}/{sample}/assembly_flye"
    # conda:
    #     "../envs/flye.yaml"
    log:
        "{outdir}/{sample}/logs/flye.log"
    shell:
        """
        flye \
            --genome-size 3m \
            --scaffold \
            --out-dir {params.dir} \
            --nano-raw {input} &> {log}
        """


# rule get_fasta:
#     input:
#         "{outdir}/{sample}/reads/circle-template.fastq"
#     output:
#         "{outdir}/{sample}/tmp/reads.fasta"
#     conda:
#         "../envs/circlator.yaml"
#     shell:
#         """
#         seqtk seq -a {input} > {output}
#         """
#
# rule circlator_assemble_contigs:
#     input:
#         fasta="{outdir}/{sample}/tmp/reads.fasta"
#     output:
#         file="{outdir}/{sample}/assembly_circlator/contigs.fasta"
#     params:
#         dir="{outdir}/{sample}/assembly_circlator"
#     threads: THREADS
#     conda:
#         "../envs/circlator.yaml"
#     log:
#         "{outdir}/{sample}/logs/circlator_contigs.log"
#     shell:
#         """
#         mkdir -p {params.dir}
#         circlator assemble \
#             --spades_k 21,33,55,77 \
#             --threads {threads} \
#             --spades_use_first \
#             --data_type nanopore-raw {input.fasta} {params.dir} &> {log}
#         """
# #  77,21,107,97
#
# # remove very short contigs
# rule circlator_clean_contigs:
#     input:
#         "{outdir}/{sample}/assembly_circlator/contigs.fasta"
#     output:
#         "{outdir}/{sample}/assembly_circlator/clean.fasta"
#     threads: THREADS
#     conda:
#         "../envs/circlator.yaml"
#     params:
#         prefix="{outdir}/{sample}/assembly_circlator/clean"
#     log:
#         "{outdir}/{sample}/logs/circlator_clean.log"
#     shell:
#         """
#         circlator clean {input} {params.prefix} &> {log}
#         """
#
# # run circularization pipeline
# rule circlator_pipeline:
#     input:
#         contigs="{outdir}/{sample}/assembly_circlator/contigs.fasta",
#         reads="{outdir}/{sample}/tmp/reads.fasta"
#     output:
#         file="{outdir}/{sample}/assembly_circlator_all/06.fixstart.fasta"
#     params:
#         dir="{outdir}/{sample}/assembly_circlator_all"
#     threads: THREADS
#     conda:
#         "../envs/circlator.yaml"
#     log:
#         "{outdir}/{sample}/logs/circlator_pipeline.log"
#     shell:
#         """
#         mkdir -p {params.dir}
#         circlator all --assembler canu --data_type nanopore-raw --split_all_reads {input.contigs} {input.reads} {params.dir} &> {log}
#         """
