rule pb_collectmultiplemetrics:
    input:
        ref = hg38_reference_genome,
        bam = rules.pb_germline_FE.output.outBam,
        bai = rules.pb_germline_FE.output.outBai,
    output:
        qc_metrics_dir = directory(config['result_folder'] + "/pb_collectmultiplemetrics/QC_Metrics_{prefix}"),
    params:
        tmpDir = config['SCRATCH1_pbTemp'] + "/pb_collectmultiplemetrics/{prefix}",
        DOCKER = config['pb_collectmultiplemetrics']['docker'],
        MEM = config['pb_collectmultiplemetrics']['mem'],
        LSF_LOG = config['logs_folder'] + "/pb_collectmultiplemetrics/LSF_{prefix}.log",
        JOBNAME = "pb_collectmultiplemetrics_{prefix}",
        GPU = "-gpu 'num=" + str(config['pb_collectmultiplemetrics']['nvidia_gpu']) + ":gmodel=TeslaV100_SXM2_32GB:j_exclusive=yes'",
        RESOURCE = "-R 'gpuhost rusage[mem=" + config['pb_collectmultiplemetrics']['mem'] + "] span[hosts=1]'",
    log: config['logs_folder'] + "/pb_collectmultiplemetrics/{prefix}.log"
    threads: config['pb_collectmultiplemetrics']['threads']
    resources:
        nvidia_gpu = config['pb_collectmultiplemetrics']['nvidia_gpu'],
        gpu_model = config['gpuModel']
    shell:
        "mkdir -p {params.tmpDir} && "
        "pbrun collectmultiplemetrics --x4 "
        "--ref {input.ref} "
        "--bam {input.bam} "
        "--out-qc-metrics-dir {output.qc_metrics_dir} "
        "--tmp-dir {params.tmpDir} "
        "--processor-threads {threads} "
        "--bam-decompressor-threads {threads} "
        "> {log} 2>&1"
        
# --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms48G -Xmx48G -XX:ParallelGCThreads=8"
# gatk CollectMultipleMetrics --REFERENCE_SEQUENCE Ref.fa -I wgs.bam -O metrics \
# --PROGRAM CollectAlignmentSummaryMetrics \
# --PROGRAM CollectInsertSizeMetrics \
# --PROGRAM QualityScoreDistribution \
# --PROGRAM MeanQualityByCycle \
# --PROGRAM CollectBaseDistributionByCycle \
# --PROGRAM CollectGcBiasMetrics \
# --PROGRAM CollectSequencingArtifactMetrics \
# --PROGRAM CollectQualityYieldMetrics