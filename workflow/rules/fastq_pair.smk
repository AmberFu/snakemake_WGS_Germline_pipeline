### Follow /storage1/fs1/jin810/Active/fup/UDN_CP_download_for_Nahyun/repair_FQs_BSUBs/generate_BSUB_repairFQ_germline.sh
rule unzip_fastq:
    input:
        fqR1_gz = lambda wildcards: config['fastq_folder'] + prefix_dict[wildcards.prefix] + "_r1.fastq.gz",
        fqR2_gz = lambda wildcards: config['fastq_folder'] + prefix_dict[wildcards.prefix] + "_r2.fastq.gz",
    output:
        fqR1 = temp(config['result_folder'] + "/unzipped_FQs/{prefix}_r1.fastq"),
        fqR2 = temp(config['result_folder'] + "/unzipped_FQs/{prefix}_r2.fastq"),
    params:
        DOCKER = config['unzip_fastq']['docker'],
        MEM = config['unzip_fastq']['mem'],
        LSF_LOG = config['logs_folder'] + "/unzip_fastq/LSF_{prefix}.log",
        JOBNAME = "unzip_fastq_{prefix}",
    log: config['logs_folder'] + "/unzip_fastq/{prefix}.log"
    threads: config['unzip_fastq']['threads']
    shell:
        "gunzip -c {input.fqR1_gz} > {output.fqR1} 2>{log} && "
        "gunzip -c {input.fqR2_gz} > {output.fqR2} 2>>{log}"

### Run very long time...
rule fastq_pair:
    input:
        fqR1 = rules.unzip_fastq.output.fqR1,
        fqR2 = rules.unzip_fastq.output.fqR2,
    output:
        fqR1_paired = temp(config['result_folder'] + "/fastq_pair/{prefix}_r1.fastq.paired.fq"),
        fqR2_paired = temp(config['result_folder'] + "/fastq_pair/{prefix}_r2.fastq.paired.fq"),
    params:
        t_number = "160000000"
        DOCKER = config['fastq_pair']['docker'],
        MEM = config['fastq_pair']['mem'],
        LSF_LOG = config['logs_folder'] + "/fastq_pair/LSF_{prefix}.log",
        JOBNAME = "fastq_pair_{prefix}",
    log: config['logs_folder'] + "/fastq_pair/{prefix}.log"
    threads: config['fastq_pair']['threads']
    shell:
        "fastq_pair -t {params.t_number} {input.fqR1} {input.fqR2}"
        "> {log} 2>&1"

rule gzip_repaired_fastq:
    input:
        fqR1_paired = rules.fastq_pair.output.fqR1_paired,
        fqR2_paired = rules.fastq_pair.output.fqR2_paired,
    output:
        fqR1_paired_zip = config['result_folder'] + "/gzip_repaired_fastq/{prefix}_r1.fastq.paired.fq.gz",
        fqR2_paired_zip = config['result_folder'] + "/gzip_repaired_fastq/{prefix}_r2.fastq.paired.fq.gz",
    params:
        DOCKER = config['gzip_repaired_fastq']['docker'],
        MEM = config['gzip_repaired_fastq']['mem'],
        LSF_LOG = config['logs_folder'] + "/gzip_repaired_fastq/LSF_{prefix}.log",
        JOBNAME = "gzip_repaired_fastq_{prefix}",
    log: config['logs_folder'] + "/gzip_repaired_fastq/{prefix}.log"
    threads: config['gzip_repaired_fastq']['threads']
    shell: 
        "gzip -c {input.fqR1_paired} > {output.fqR1_paired_zip} 2>{log} && "
        "gzip -c {input.fqR2_paired} > {output.fqR2_paired_zip} 2>>{log} "