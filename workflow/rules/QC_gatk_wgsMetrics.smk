rule QC_gatk_CollectWgsMetrics:
    input:
        ref = hg38_reference_genome,
        bam = rules.pb_germline_FE.output.outBam,
        bai = rules.pb_germline_FE.output.outBai,
    output:
        wgs_metrics = config['result_folder'] + "/QC_gatk_CollectWgsMetrics/{prefix}_bamMetrics.txt",
    params:
        DOCKER = config['QC_gatk_CollectWgsMetrics']['docker'],
        MEM = config['QC_gatk_CollectWgsMetrics']['mem'],
        LSF_LOG = config['logs_folder'] + "/QC_gatk_CollectWgsMetrics/LSF_{prefix}.log",
        JOBNAME = "QC_gatk_CollectWgsMetrics_{prefix}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['QC_gatk_CollectWgsMetrics']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
        mem_gatk = str(config['QC_gatk_CollectWgsMetrics']['mem']).split('B')[0],
    log: config['logs_folder'] + "/QC_gatk_CollectWgsMetrics/{prefix}.log"
    threads: config['QC_gatk_CollectWgsMetrics']['threads']
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk --java-options '-Xms{params.mem_gatk} -Xmx{params.mem_gatk} -XX:ParallelGCThreads={threads}' CollectWgsMetrics "
        "-R {input.ref} "
        "-I {input.bam} "
        "-O {output.wgs_metrics} "
        "> {log} 2>&1 "


# --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms48G -Xmx48G -XX:ParallelGCThreads=8"
# gatk CollectWgsMetrics -R Ref.fa -I wgs.bam -O metrics.txt
# TEST: time /gatk/gatk --java-options "-Xms10G -Xmx10G -XX:ParallelGCThreads=10" CollectWgsMetrics -R /storage1/fs1/jin810/Active/references/Homo_sapiens_assembly38.fasta -I data/germline_output/WGS_001_father_UDN121697.bam -O data/germline_output/WGS_001_father_UDN121697_CollectWGSMetrics.txt

## Cram also runable:
#/gatk/gatk CollectWgsMetrics -I /storage1/fs1/jin810/Active/fup/Kruer_CP_WGS/snakemake_WGS_FE_pipeline_call_denovo/results/bam2cram/CP-F0309/CP-F0309-001U_FE_Parabricks_Germline_bam.cram -O /storage1/fs1/jin810/Active/fup/Kruer_CP_WGS/test_CollectWgsMetrics_on_CP-F0309-001U.txt -R /storage1/fs1/jin810/Active/references/Homo_sapiens_assembly38.fasta