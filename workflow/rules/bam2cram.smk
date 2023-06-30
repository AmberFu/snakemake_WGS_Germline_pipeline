rule bam2cram:
    input:
        ref = hg38_reference_genome,
        #bam = config['result_folder']+"/pb_germline_FE/{prefix}.bam",
        #bai = config['result_folder']+"/pb_germline_FE/{prefix}.bam.bai",
        bam = rules.pb_germline_FE.output.outBam,
        bai = rules.pb_germline_FE.output.outBai,
    output:
        cram = config['result_folder'] + "/bam2cram/{prefix}_bam.cram",
    params:
        DOCKER = config['bam2cram']['docker'],
        MEM = config['bam2cram']['mem'],
        LSF_LOG = config['logs_folder'] + "/bam2cram/LSF_{prefix}.log",
        JOBNAME = "bam2cram_{prefix}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['bam2cram']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    threads: config['bam2cram']['threads']
    log: config['logs_folder']+"/bam2cram/{prefix}.log"
    resources:
        tmpdir=config['ACTIVE_TMP']
    shell: 
        "samtools view -C -@ {threads} -T {input.ref} "
        "-o {output.cram} {input.bam} > {log} 2>&1"


rule index_cram:
    input: rules.bam2cram.output.cram
    output: 
        crai = config['result_folder'] + "/bam2cram/{prefix}_bam.cram.crai"
    params:
        DOCKER = config['index_cram']['docker'],
        MEM = config['index_cram']['mem'],
        LSF_LOG = config['logs_folder'] + "/index_cram/LSF_{prefix}.log",
        JOBNAME = "index_cram_{prefix}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['index_cram']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    threads: config['index_cram']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    log: config['logs_folder']+"/index_cram/{prefix}.log"
    shell: 
        "samtools index -@ {threads} -c "
        "{input} > {log} 2>&1"

### 
PREFIX = list(SAMPLE_DICT.keys())


rule md5CRAM:
    input:
        cram = expand(config['result_folder'] + "/bam2cram/{prefix}_bam.cram", prefix=PREFIX),
        crai = expand(config['result_folder'] + "/bam2cram/{prefix}_bam.cram.crai", prefix=PREFIX),
    output: config['result_folder'] + "/bam2cram/crams.md5"
    params:
       DOCKER = config['md5CRAM']['docker'],
       MEM = config['md5CRAM']['mem'],
       LSF_LOG = config['logs_folder'] + "/md5CRAM/LSF_{prefix}.log",
       JOBNAME = "md5CRAM_{prefix}",
       GPU = "",
       RESOURCE = "-R 'rusage[mem=" + config['md5CRAM']['mem'] + "] span[hosts=1]'",
       QUEUE = config['general_queue'],
       GROUP = config['general_compute_group'],
    threads: config['md5CRAM']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    shell: "md5sum {input.cram} {input.crai} > {output}"


rule tarCRAM:
    input: 
        cram = expand(config['result_folder'] + "/bam2cram/{prefix}_bam.cram", prefix=PREFIX),
        crai = expand(config['result_folder'] + "/bam2cram/{prefix}_bam.cram.crai", prefix=PREFIX),
        md5 = rules.md5CRAM.output,
    output: config['result_folder'] + "/" + config['project_name'] + "_CRAMs.tar.gz"
    params:
       DOCKER = config['tarCRAM']['docker'],
       MEM = config['tarCRAM']['mem'],
       LSF_LOG = config['logs_folder'] + "/tarCRAM/LSF_Tar_CRAMs.log",
       JOBNAME = "tarCRAM_{prefix}",
       GPU = "",
       RESOURCE = "-R 'rusage[mem=" + config['tarCRAM']['mem'] + "] span[hosts=1]'",
       QUEUE = config['general_queue'],
       GROUP = config['general_compute_group'],
    threads: config['tarCRAM']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    log: config['logs_folder']+"/tarCRAM/Tar_CRAMs.log"
    shell: "tar zcvf {output} {input.cram} {input.crai} {input.md5} 2> {log}"

