### Prepare a Dict:
prefix_dict = {}
for index, row in SAMPLES.iterrows():
    prefix_dict[row.prefix] = row.BAMCRAM_STEMNAME
    
rule pb_germline_FE:
    input:
        ref = hg38_reference_genome,
        fqR1 = lambda wildcards: config['result_folder'] + "/gzip_fastq/" + prefix_dict[wildcards.prefix] + "_r1.fastq.gz",
        fqR2 = lambda wildcards: config['result_folder'] + "/gzip_fastq/" + prefix_dict[wildcards.prefix] + "_r2.fastq.gz",
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
        LSF_LOG = config['logs_folder'] + "/pb_germline_FE/LSF_{prefix}.log",
        JOBNAME = "pb_germline_{prefix}",
        GPU = "-gpu 'num=" + str(config['pb_germline_FE']['nvidia_gpu']) + ":gmodel=TeslaV100_SXM2_32GB:j_exclusive=yes'",
        RESOURCE = "-R 'gpuhost rusage[mem=" + config['pb_germline_FE']['mem'] + "] span[hosts=1]'",
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

