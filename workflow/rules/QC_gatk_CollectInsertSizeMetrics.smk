'''
This step is going to get insertion size infomation. (Cut-off: Median insert size: < 250)

GATK4 CollectInsertSizeMetrics: 
https://gatk.broadinstitute.org/hc/en-us/articles/360037055772-CollectInsertSizeMetrics-Picard-
'''

rule QC_gatk_CollectInsertSizeMetrics:
    input:
        ref = hg38_reference_genome,
        bam = rules.pb_germline_FE.output.outBam,
        bai = rules.pb_germline_FE.output.outBai,
        #bam = config['result_folder']+"/pb_germline_FE/{prefix}.bam",
        #bai = config['result_folder']+"/pb_germline_FE/{prefix}.bam.bai",
        #cram = rules.bam2cram.output.cram,
        #crai = rules.index_cram.output.crai,
    output:
        insertion_size_metrics = config['result_folder'] + "/QC_gatk_CollectInsertSizeMetrics/{prefix}_InsertSizeMetrics.txt",
        insert_size_histogram_pdf = config['result_folder'] + "/QC_gatk_CollectInsertSizeMetrics/{prefix}_InsertSizeHistogram.pdf",
    params:
        DOCKER = config['QC_gatk_CollectInsertSizeMetrics']['docker'],
        MEM = config['QC_gatk_CollectInsertSizeMetrics']['mem'],
        LSF_LOG = config['logs_folder'] + "/QC_gatk_CollectInsertSizeMetrics/LSF_{prefix}.log",
        JOBNAME = "QC_gatk_CollectInsertSizeMetrics_{prefix}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['QC_gatk_CollectInsertSizeMetrics']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
        mem_gatk = str(config['QC_gatk_CollectInsertSizeMetrics']['mem']).split('B')[0],
    log: config['logs_folder'] + "/QC_gatk_CollectInsertSizeMetrics/{prefix}.log"
    threads: config['QC_gatk_CollectInsertSizeMetrics']['threads']
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk --java-options '-Xms{params.mem_gatk} -Xmx{params.mem_gatk} -XX:ParallelGCThreads={threads}' CollectInsertSizeMetrics "
        "-R {input.ref} "
        "-I {input.bam} "
        "-O {output.insertion_size_metrics} "
        "-H {output.insert_size_histogram_pdf} "
        "> {log} 2>&1 "






'''
CollectInsertSizeMetrics

java -jar picard.jar CollectInsertSizeMetrics \
      I=input.bam \
      O=insert_size_metrics.txt \
      H=insert_size_histogram.pdf \
      M=0.5 # if run small file
'''