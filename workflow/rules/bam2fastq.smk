#############
## BAM
#############
## Sort BAM by read name
rule collate_bam:
    input:
        bam = config['bam_cram_folder']+"/{filename}.bam",
    output:
        sorted_bam = temp(config['result_folder'] + "/collate_bam/{filename}_collate.bam"),
    params:
        prefix = config['result_folder'] + "/collate_bam/{filename}_collate",
        DOCKER = config['collate_bam']['docker'],
        MEM = config['collate_bam']['mem'],
        LSF_LOG = config['logs_folder'] + "/collate_bam/LSF_{filename}.log",
        JOBNAME = "collate_bam_{filename}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['collate_bam']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    log: config['logs_folder'] + "/collate_bam/{filename}.log"
    threads: config['collate_bam']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    shell:
        "samtools collate -@ {threads} -o {output.sorted_bam} {input.bam} {params.prefix} >> {log} 2>&1"


rule bam2fastq:
    input:
        sorted_bam = rules.collate_bam.output.sorted_bam,
    output:
        r1 = temp(config['result_folder']+"/bam2fastq/{filename}_r1.fastq"),
        r2 = temp(config['result_folder']+"/bam2fastq/{filename}_r2.fastq"),
    log: config['logs_folder']+"/bam2fastq/{filename}.log"
    params:
        DOCKER = config['bam2fastq']['docker'],
        MEM = config['bam2fastq']['mem'],
        LSF_LOG = config['logs_folder'] + "/bam2fastq/LSF_{filename}.log",
        JOBNAME = "bam2fastq_{filename}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['bam2fastq']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    threads: config['bam2fastq']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    shell:
        "samtools fastq -@ {threads} -1 {output.r1} -2 {output.r2} "
        "-0 /dev/null -s /dev/null -n "
        "{input.sorted_bam} > {log} 2>&1"

rule gzip_r1_fastq_bam:
    input:
        r1 = rules.bam2fastq.output.r1,
    output:
        r1gz = config['result_folder']+"/gzip_fastq/{filename}_r1.fastq.gz",
    log: config['logs_folder']+"/gzip_r1_fastq_bam/{filename}.log"
    params:
        DOCKER = config['gzip_fastq_bam']['docker'],
        MEM = config['gzip_fastq_bam']['mem'],
        LSF_LOG = config['logs_folder'] + "/gzip_r1_fastq_bam/LSF_{filename}.log",
        JOBNAME = "gzip_r1_fastq_bam_{filename}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['gzip_fastq_bam']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    threads: config['gzip_fastq_bam']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    shell:
        "gzip -c {input.r1} > {output.r1gz} 2>{log}"
        
        
rule gzip_r2_fastq_bam:
    input:
        r2 = rules.bam2fastq.output.r2,
    output:
        r2gz = config['result_folder']+"/gzip_fastq/{filename}_r2.fastq.gz",
    log: config['logs_folder']+"/gzip_r2_fastq_bam/{filename}.log"
    params:
        DOCKER = config['gzip_fastq_bam']['docker'],
        MEM = config['gzip_fastq_bam']['mem'],
        LSF_LOG = config['logs_folder'] + "/gzip_r2_fastq_bam/LSF_{filename}.log",
        JOBNAME = "gzip_r2_fastq_bam_{filename}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['gzip_fastq_bam']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    threads: config['gzip_fastq_bam']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    shell:
        "gzip -c {input.r2} > {output.r2gz} 2>{log}"



#############
## CRAM
#############
rule collate_cram:
    input:
        cram = config['bam_cram_folder']+"/{filename}.cram",
    output:
        sorted_cram = temp(config['result_folder'] + "/collate_cram/{filename}_collate.bam"),
    params:
        DOCKER = config['collate_cram']['docker'],
        MEM = config['collate_cram']['mem'],
        LSF_LOG = config['logs_folder'] + "/collate_cram/LSF_{filename}.log",
        JOBNAME = "collate_cram_{filename}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['collate_cram']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    log: config['logs_folder'] + "/collate_cram/{filename}.log"
    threads: config['collate_cram']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    shell:
        "samtools collate -@ {threads} -o {output.sorted_cram} {input.cram} > {log} 2>&1"

rule cram2fastq:
    input:
        sorted_cram = rules.collate_cram.output.sorted_cram,
    output:
        r1 = temp(config['result_folder']+"/cram2fastq/{filename}_r1.fastq"),
        r2 = temp(config['result_folder']+"/cram2fastq/{filename}_r2.fastq"),
    log:
        config['logs_folder']+"/cram2fastq/{filename}.log"
    params:
        DOCKER = config['cram2fastq']['docker'],
        MEM = config['cram2fastq']['mem'],
        LSF_LOG = config['logs_folder'] + "/cram2fastq/LSF_{filename}.log",
        JOBNAME = "cram2fastq_{filename}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['cram2fastq']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    threads: config['cram2fastq']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    shell:
        "samtools fastq -@ {threads} -1 {output.r1} -2 {output.r2} "
        "-0 /dev/null -s /dev/null -n "
        "{input.sorted_cram} > {log} 2>&1"


rule gzip_r1_fastq_cram:
    input:
        r1 = rules.cram2fastq.output.r1,
    output:
        r1gz = config['result_folder']+"/gzip_fastq/{filename}_r1.fastq.gz",
    log: config['logs_folder']+"/gzip_r1_fastq_cram/{filename}.log"
    params:
        DOCKER = config['gzip_fastq_cram']['docker'],
        MEM = config['gzip_fastq_cram']['mem'],
        LSF_LOG = config['logs_folder'] + "/gzip_r1_fastq_cram/LSF_{filename}.log",
        JOBNAME = "gzip_r1_fastq_cram_{filename}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['gzip_fastq_cram']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    threads: config['gzip_fastq_cram']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    shell:
        "gzip -c {input.r1} > {output.r1gz} 2>{log}"


rule gzip_r2_fastq_cram:
    input:
        r2 = rules.cram2fastq.output.r2,
    output:
        r2gz = config['result_folder']+"/gzip_fastq/{filename}_r2.fastq.gz",
    log: config['logs_folder']+"/gzip_r2_fastq_cram/{filename}.log"
    params:
        DOCKER = config['gzip_fastq_cram']['docker'],
        MEM = config['gzip_fastq_cram']['mem'],
        LSF_LOG = config['logs_folder'] + "/gzip_r2_fastq_cram/LSF_{filename}.log",
        JOBNAME = "gzip_r2_fastq_cram_{filename}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['gzip_fastq_cram']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    threads: config['gzip_fastq_cram']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    shell:
        "gzip -c {input.r2} > {output.r2gz} 2>{log}"
        
        
        

### SAMTOOLS doc: http://www.htslib.org/doc/samtools-fasta.html
# samtools collate -u -O in_pos.bam | \\
# samtools fastq -1 paired1.fq -2 paired2.fq -0 /dev/null -s /dev/null -n