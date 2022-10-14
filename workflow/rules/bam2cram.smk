rule bam2cram:
    input:
        ref = hg38_reference_genome,
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
    threads: config['bam2cram']['threads']
    log: config['logs_folder']+"/bam2cram/{prefix}.log"
    shell: 
        "samtools view -C -@ {threads} -T {input.ref} "
        "-o {output.cram} {input.bam} > {log} 2>&1"


rule index_cram:
    input: rules.bam2cram.output.cram
    output: config['result_folder'] + "/bam2cram/{prefix}_bam.cram.crai"
    params:
       DOCKER = config['index_cram']['docker'],
       MEM = config['index_cram']['mem'],
       LSF_LOG = config['logs_folder'] + "/index_cram/LSF_{prefix}.log",
       JOBNAME = "index_cram_{prefix}",
       GPU = "",
       RESOURCE = "-R 'rusage[mem=" + config['index_cram']['mem'] + "] span[hosts=1]'",
    threads: config['index_cram']['threads']
    log: config['logs_folder']+"/index_cram/{prefix}.log"
    shell: 
        "samtools index -@ {threads} -c "
        "{input} > {log} 2>&1"

### 
PREFIX = []
for index, row in SAMPLES.iterrows():
    PREFIX.append(row['prefix'])


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
    threads: config['md5CRAM']['threads']
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
    threads: config['tarCRAM']['threads']
    log: config['logs_folder']+"/tarCRAM/Tar_CRAMs.log"
    shell: "tar zcvf {output} {input.cram} {input.crai} {input.md5} 2> {log}"

