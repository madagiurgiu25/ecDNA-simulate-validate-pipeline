
rule sniffles:
    input:
        bam=lambda wildcards: "{outdir}/{sample}/mapping/map.bam" if RUNMODE in [RUNFULL,RUNMAPPING, RUNSIMTOSV, RUNSIMPLE] else BAM,
        bai=lambda wildcards: "{outdir}/{sample}/mapping/map.bam.bai" if RUNMODE in [RUNFULL,RUNMAPPING, RUNSIMTOSV, RUNSIMPLE] else BAM + ".bai",
    output:
        "{outdir}/{sample}/sv/sv.sniffles.vcf"
#    conda:
#        "../envs/sv.yaml"
    params:
        threads = THREADS
    log:
        "{outdir}/{sample}/logs/sniffles.log"
    shell:
        """
        sniffles -t {params.threads} -m {input.bam} \
        -v {output} --min_homo_af 0.7 \
        --min_het_af 0.1 --min_length 50 \
        --cluster --genotype --min_support 4 --report-seq
        """


#         """
#	 # sniffles2
#         sniffles -t {params.threads} -i {input.bam} -v {output} \
#            --minsvlen 50 \
#            --cluster-merge-bnd 300 --minsupport 4 --non-germline &> {log}
#        """

rule survivor:
    priority: 1
    input:
        unfilt="{outdir}/{sample}/sv/sv.sniffles.vcf"
    output:
        unfilt="{outdir}/{sample}/sv/sv.sniffles.bedpe"
#    conda:
#        "../envs/sv.yaml"
    log:
        "{outdir}/{sample}/logs/survivor.log"
    shell:
        """
        SURVIVOR vcftobed {input.unfilt} -1 -1 {output.unfilt} &> {log}
        """

# rule svim_reads:
#     input:
#         reads="{outdir}/{sample}/reads/circle-template.fastq",
#         reference=REF
#     output:
#         final="{outdir}/{sample}/sv/sv.svim_reads.vcf",
#         dir=directory("{outdir}/{sample}/sv/svim_reads")
#     conda:
#         "../envs/sv.yaml"
#     params:
#         qual = READ_QUAL
#     log:
#         "{outdir}/{sample}/logs/svim_reads.log"
#     shell:
#          """
#          mkdir -p {output.dir} && \
#          svim reads --read_names --nanopore --aligner ngmlr \
#             --min_mapq {params.qual} --max_sv_size 1000000 \
#             --heterozygous_threshold 0.1 --homozygous_threshold 0.7 \
#             --sample {wildcards.sample} {output.dir} {input.reads} {input.reference}  &> {log} && \
#          mv {output.dir}/variants.vcf {output.final}
#          """

# rule svim_align:
#     input:
#         bam="{outdir}/{sample}/mapping/ngmlr.bam",
#         bai="{outdir}/{sample}/mapping/ngmlr.bam.bai",
#         reference=REF
#     output:
#         final="{outdir}/{sample}/sv/sv.svim.vcf",
#         dir=directory("{outdir}/{sample}/sv/svim_align")
#     conda:
#         "../envs/sv.yaml"
#     log:
#         "{outdir}/{sample}/logs/svim.log"
#     shell:
#          """
#          mkdir -p {output.dir} && \
#          svim alignment --read_names --sample {wildcards.sample} \
#             --max_sv_size 1000000 \
#             --heterozygous_threshold 0.1 --homozygous_threshold 0.7 \
#             {output.dir} {input.bam} {input.reference} &> {log} && \
#          cp {output.dir}/variants.vcf {output.final}
#          """
