'''
SAMPLE_DICT = {'PREFIX_01': {'R1': 'xxx_r1.fastq.gz', 'R2': 'xxx_r2.fastq.gz'},
               'PREFIX_02': {'R1': 'xox_r1.fastq.gz', 'R2': 'xox_r2.fastq.gz'},
               'PREFIX_03': {'R1': 'xxo_r1.fastq.gz', 'R2': 'xxo_r2.fastq.gz'}}
'''

rule pb_germline_FE:
    input:
        ref = hg38_reference_genome,
        fqR1 = lambda wildcards: SAMPLE_DICT[wildcards.prefix]['R1'],
        fqR2 = lambda wildcards: SAMPLE_DICT[wildcards.prefix]['R2'],
        knownSite19443 = KnownSite19443,
        knownSite20267 = KnownSite20267,
        knownSite20211 = KnownSite20211,
    output:
        outVar = config['result_folder'] + "/pb_germline_FE/{prefix}.g.vcf",
        outVarIdx = config['result_folder'] + "/pb_germline_FE/{prefix}.g.vcf.idx",
        outBam = config['result_folder'] + "/pb_germline_FE/{prefix}.bam",
        outBai = config['result_folder'] + "/pb_germline_FE/{prefix}.bam.bai",
        outRecal = config['result_folder'] + "/pb_germline_FE/{prefix}_report.txt",
        outDupMetrics = config['result_folder'] + "/pb_germline_FE/{prefix}_dup_metrics.txt"
    params:
        DOCKER = config['pb_germline_FE']['docker'],
        MEM = config['pb_germline_FE']['mem'],
        MEM_LIMIT = config['pb_germline_FE']['mem_limit'],
        LSF_LOG = config['logs_folder'] + "/pb_germline_FE/LSF_{prefix}.log",
        JOBNAME = "pb_germline_{prefix}",
        GPU = "-gpu 'num=" + str(config['pb_germline_FE']['nvidia_gpu']) + ":gmodel=TeslaV100_SXM2_32GB:j_exclusive=yes'",
        RESOURCE = "-R 'gpuhost rusage[mem=" + config['pb_germline_FE']['mem'] + "] span[hosts=1]'",
        QUEUE = config['subscription_queue'],
        GROUP = config['subscription_compute_group'] + " -sla " + config['subscription_sla'],
        tmpDir = config['SCRATCH1_pbTemp'] + "/pb_germline_FE/{prefix}",
        readGroupPL = config['platform'],
        readGroupSM = "{prefix}"
    log: config['logs_folder'] + "/pb_germline_FE/{prefix}.log"
    threads: config['pb_germline_FE']['threads']
    resources:
        nvidia_gpu = config['pb_germline_FE']['nvidia_gpu'],
        gpu_model = config['gpuModel']
    shell:
        "mkdir -p {params.tmpDir} && "
        "pbrun germline --x4 "
        "--memory-limit {params.MEM_LIMIT} "
        "--ref {input.ref} "
        "--in-fq {input.fqR1} {input.fqR2} "
        "--knownSites {input.knownSite19443} "
        "--knownSites {input.knownSite20267} "
        "--knownSites {input.knownSite20211} "
        "--num-gpus {resources.nvidia_gpu} "
        "--out-variants {output.outVar} "
        "--out-bam {output.outBam} "
        "--out-recal-file {output.outRecal} "
        "--out-duplicate-metrics {output.outDupMetrics} "
        "--tmp-dir {params.tmpDir} "
        "--read-group-pl {params.readGroupPL} "
        "--read-group-sm {params.readGroupSM} "
        "--static-quantized-quals 10 "
        "--static-quantized-quals 20 "
        "--static-quantized-quals 30 "
        "--bwa-options='-Y -K 100000000' "
        "--gvcf "
        "> {log} 2>&1 "

