'''
Re-pair each Raw FASTQ:

ie, RAWFQs_DICT = 
{'HG002_NA24385_son_D1_S1_L001_001': 
 {'r1': '/storage1/fs1/jin810/Active/datasets/NIST_AshkenazimTrio_Illumina_2_250bps/HG002_NA24385_son/D1_S1_L001_R1_001.fastq.gz',
  'r2': '/storage1/fs1/jin810/Active/datasets/NIST_AshkenazimTrio_Illumina_2_250bps/HG002_NA24385_son/D1_S1_L001_R2_001.fastq.gz'},
 'HG002_NA24385_son_D1_S1_L001_002': 
 {'r1': '/storage1/fs1/jin810/Active/datasets/NIST_AshkenazimTrio_Illumina_2_250bps/HG002_NA24385_son/D1_S1_L001_R1_002.fastq.gz',
  'r2': '/storage1/fs1/jin810/Active/datasets/NIST_AshkenazimTrio_Illumina_2_250bps/HG002_NA24385_son/D1_S1_L001_R2_002.fastq.gz'},
  ...}
  
=> RAW_PREFIX: ie, HG002_NA24385_son_D1_S1_L001_001
'''


rule unzip_fastq_R1:
    input:
        fqR1_gz = lambda wildcards: RAWFQs_DICT[wildcards.raw_prefix]['r1'],
    output:
        fqR1 = temp(config['result_folder'] + "/fastq_pair/{raw_prefix}_R1.fastq"),
    params:
        DOCKER = config['unzip_fastq_R1']['docker'],
        MEM = config['unzip_fastq_R1']['mem'],
        LSF_LOG = config['logs_folder'] + "/unzip_fastq_R1/LSF_{raw_prefix}.log",
        JOBNAME = "unzip_fastq_R1_{raw_prefix}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['unzip_fastq_R1']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    log: config['logs_folder'] + "/unzip_fastq_R1/{raw_prefix}.log"
    threads: config['unzip_fastq_R1']['threads']
    shell:
        "gunzip -c {input.fqR1_gz} > {output.fqR1} 2>{log}"

rule unzip_fastq_R2:
    input:
        fqR2_gz = lambda wildcards: RAWFQs_DICT[wildcards.raw_prefix]['r2'],
    output:
        fqR2 = temp(config['result_folder'] + "/fastq_pair/{raw_prefix}_R2.fastq"),
    params:
        DOCKER = config['unzip_fastq_R2']['docker'],
        MEM = config['unzip_fastq_R2']['mem'],
        LSF_LOG = config['logs_folder'] + "/unzip_fastq_R2/LSF_{raw_prefix}.log",
        JOBNAME = "unzip_fastq_R2_{raw_prefix}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['unzip_fastq_R2']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    log: config['logs_folder'] + "/unzip_fastq_R2/{raw_prefix}.log"
    threads: config['unzip_fastq_R2']['threads']
    shell:
        "gunzip -c {input.fqR2_gz} > {output.fqR2} 2>{log}"


rule fastq_pair:
    input:
        fqR1 = rules.unzip_fastq_R1.output.fqR1,
        fqR2 = rules.unzip_fastq_R2.output.fqR2,
    output:
        fqR1_paired = temp(config['result_folder'] + "/fastq_pair/{raw_prefix}_R1.fastq.paired.fq"),
        fqR2_paired = temp(config['result_folder'] + "/fastq_pair/{raw_prefix}_R2.fastq.paired.fq"),
    params:
        t_number = "160000000",
        DOCKER = config['fastq_pair']['docker'],
        MEM = config['fastq_pair']['mem'],
        LSF_LOG = config['logs_folder'] + "/fastq_pair/LSF_{raw_prefix}.log",
        JOBNAME = "fastq_pair_{raw_prefix}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['fastq_pair']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    log: config['logs_folder'] + "/fastq_pair/{raw_prefix}.log"
    threads: config['fastq_pair']['threads']
    shell:
        "fastq_pair -t {params.t_number} {input.fqR1} {input.fqR2}"
        "> {log} 2>&1"

rule gzip_repaired_fastq_R1:
    input:
        fqR1_paired = rules.fastq_pair.output.fqR1_paired,
    output:
        fqR1_paired_zip = config['result_folder']+"/gzip_repaired_fastq/{raw_prefix}_R1.fastq.paired.fq.gz",
    params:
        DOCKER = config['gzip_repaired_fastq_R1']['docker'],
        MEM = config['gzip_repaired_fastq_R1']['mem'],
        LSF_LOG = config['logs_folder'] + "/gzip_repaired_fastq_R1/LSF_{raw_prefix}.log",
        JOBNAME = "gzip_repaired_fastq_R1_{raw_prefix}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['gzip_repaired_fastq_R1']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    log: config['logs_folder'] + "/gzip_repaired_fastq_R1/{raw_prefix}.log"
    threads: config['gzip_repaired_fastq_R1']['threads']
    shell: 
        "gzip -c {input.fqR1_paired} > {output.fqR1_paired_zip} 2>{log}"
        
        
rule gzip_repaired_fastq_R2:
    input:
        fqR2_paired = rules.fastq_pair.output.fqR2_paired,
    output:
        fqR2_paired_zip = config['result_folder']+"/gzip_repaired_fastq/{raw_prefix}_R2.fastq.paired.fq.gz",
    params:
        DOCKER = config['gzip_repaired_fastq_R2']['docker'],
        MEM = config['gzip_repaired_fastq_R2']['mem'],
        LSF_LOG = config['logs_folder'] + "/gzip_repaired_fastq_R2/LSF_{raw_prefix}.log",
        JOBNAME = "gzip_repaired_fastq_R2_{raw_prefix}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['gzip_repaired_fastq_R2']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    log: config['logs_folder'] + "/gzip_repaired_fastq_R2/{raw_prefix}.log"
    threads: config['gzip_repaired_fastq_R2']['threads']
    shell: 
        "gzip -c {input.fqR2_paired} > {output.fqR2_paired_zip} 2>{log} "
