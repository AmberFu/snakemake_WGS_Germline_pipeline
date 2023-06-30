rule pb_bammetrics:
    input:
        ref = hg38_reference_genome,
        bam = config['result_folder']+"/pb_germline_FE/{prefix}.bam",
        bai = config['result_folder']+"/pb_germline_FE/{prefix}.bam.bai",
        #bam = rules.pb_germline_FE.output.outBam,
        #bai = rules.pb_germline_FE.output.outBai,
    output:
        bam_metrics = config['result_folder'] + "/pb_bammetrics/{prefix}_bamMetrics.txt",
    params:
        tmpDir = config['SCRATCH1_pbTemp'] + "/pb_bammetrics/{prefix}",
        DOCKER = config['pb_bammetrics']['docker'],
        MEM = config['pb_bammetrics']['mem'],
        LSF_LOG = config['logs_folder'] + "/pb_bammetrics/LSF_{prefix}.log",
        JOBNAME = "pb_bammetrics_{prefix}",
        GPU = "-gpu 'num=" + str(config['pb_bammetrics']['nvidia_gpu']) + ":gmodel=TeslaV100_SXM2_32GB:j_exclusive=yes'",
        RESOURCE = "-R 'gpuhost rusage[mem=" + config['pb_bammetrics']['mem'] + "] span[hosts=1]'",
        QUEUE = config['subscription_queue'],
        GROUP = config['subscription_compute_group'] + " -sla " + config['subscription_sla'],
    log: config['logs_folder'] + "/pb_bammetrics/{prefix}.log"
    threads: config['pb_bammetrics']['threads']
    resources:
        nvidia_gpu = config['pb_bammetrics']['nvidia_gpu'],
        gpu_model = config['gpuModel']
    shell:
        "mkdir -p {params.tmpDir} && "
        "pbrun bammetrics --x4 "
        "--ref {input.ref} "
        "--bam {input.bam} "
        "--out-metrics-file {output.bam_metrics} "
        "--num-threads {threads} "
        "--tmp-dir {params.tmpDir} "
        "> {log} 2>&1 "
        
# --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms48G -Xmx48G -XX:ParallelGCThreads=8"
# gatk CollectWgsMetrics -R Ref.fa -I wgs.bam -O metrics.txt
# TEST: time /gatk/gatk --java-options "-Xms10G -Xmx10G -XX:ParallelGCThreads=10" CollectWgsMetrics -R /storage1/fs1/jin810/Active/references/Homo_sapiens_assembly38.fasta -I data/germline_output/WGS_001_father_UDN121697.bam -O data/germline_output/WGS_001_father_UDN121697_CollectWGSMetrics.txt