rule coverage:
    input:
        bam="{outdir}/{sample}/mapping/map.bam",
        bai="{outdir}/{sample}/mapping/map.bam.bai"
    output:
        "{outdir}/{sample}/mapping/coverage.bw"
    conda:
       "../envs/cov.yaml"
    log:
        cov="{outdir}/{sample}/logs/cov.log"
    params:
        threads=THREADS
    shell:
        "bamCoverage --bam {input.bam} -o {output} --binSize 50 -p {params.threads} &> {log.cov}"
