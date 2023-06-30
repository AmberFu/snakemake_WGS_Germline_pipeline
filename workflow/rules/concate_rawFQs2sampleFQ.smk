'''
SAMPLE_DICT = {'HG002_NA24385_son': RAWFQs_SON_LIST,
               'HG003_NA24149_father': RAWFQs_DAD_LIST,
               'HG004_NA24143_mother': RAWFQs_MOM_LIST}
'''

rule rawFQs_to_sampleFQ_R1:
    input:
        fqR1_paired_zip = lambda wildcards: expand(config['result_folder']+"/gzip_repaired_fastq/{raw_prefix}_R1.fastq.paired.fq.gz", raw_prefix=SAMPLE_DICT[wildcards.prefix]),
    output:
        sampleFQ_R1 = config['result_folder']+"/rawFQs_to_sampleFQ/{prefix}_R1.fastq.gz",
    params:
        DOCKER = config['rawFQs_to_sampleFQ_R1']['docker'],
        MEM = config['rawFQs_to_sampleFQ_R1']['mem'],
        LSF_LOG = config['logs_folder'] + "/rawFQs_to_sampleFQ_R1/LSF_{prefix}.log",
        JOBNAME = "rawFQs_to_sampleFQ_R1_{prefix}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['rawFQs_to_sampleFQ_R1']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    log: config['logs_folder'] + "/rawFQs_to_sampleFQ_R1/{prefix}.log"
    threads: config['rawFQs_to_sampleFQ_R1']['threads']
    shell: 
        "cat {input.fqR1_paired_zip} > {output.sampleFQ_R1} 2>{log}"
        

rule rawFQs_to_sampleFQ_R2:
    input:
        fqR2_paired_zip = lambda wildcards: expand(config['result_folder']+"/gzip_repaired_fastq/{raw_prefix}_R2.fastq.paired.fq.gz", raw_prefix=SAMPLE_DICT[wildcards.prefix]),
    output:
        sampleFQ_R2 = config['result_folder']+"/rawFQs_to_sampleFQ/{prefix}_R2.fastq.gz",
    params:
        DOCKER = config['rawFQs_to_sampleFQ_R2']['docker'],
        MEM = config['rawFQs_to_sampleFQ_R2']['mem'],
        LSF_LOG = config['logs_folder'] + "/rawFQs_to_sampleFQ_R2/LSF_{prefix}.log",
        JOBNAME = "rawFQs_to_sampleFQ_R2_{prefix}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['rawFQs_to_sampleFQ_R2']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    log: config['logs_folder'] + "/rawFQs_to_sampleFQ_R2/{prefix}.log"
    threads: config['rawFQs_to_sampleFQ_R2']['threads']
    shell: 
        "cat {input.fqR2_paired_zip} > {output.sampleFQ_R2} 2>{log}"