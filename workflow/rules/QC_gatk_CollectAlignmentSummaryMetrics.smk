'''
This step is going to generate Rates of chimeric reads (filter out > 0.05)

GATK4 CollectAlignmentSummaryMetrics: 
https://gatk.broadinstitute.org/hc/en-us/articles/360036883111-CollectAlignmentSummaryMetrics-Picard-
'''

rule QC_gatk_CollectAlignmentSummaryMetrics:
    input:
        ref = hg38_reference_genome,
        bam = rules.pb_germline_FE.output.outBam,
        bai = rules.pb_germline_FE.output.outBai,
        #bam = config['result_folder']+"/pb_germline_FE/{prefix}.bam",
        #bai = config['result_folder']+"/pb_germline_FE/{prefix}.bam.bai",
        #cram = rules.bam2cram.output.cram,
        #crai = rules.index_cram.output.crai,
    output:
        alignment_metrics = config['result_folder'] + "/QC_gatk_CollectAlignmentSummaryMetrics/{prefix}_AlignmentSummaryMetrics.txt",
    params:
        DOCKER = config['QC_gatk_CollectAlignmentSummaryMetrics']['docker'],
        MEM = config['QC_gatk_CollectAlignmentSummaryMetrics']['mem'],
        LSF_LOG = config['logs_folder'] + "/QC_gatk_CollectAlignmentSummaryMetrics/LSF_{prefix}.log",
        JOBNAME = "QC_gatk_CollectAlignmentSummaryMetrics_{prefix}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['QC_gatk_CollectAlignmentSummaryMetrics']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
        mem_gatk = str(config['QC_gatk_CollectAlignmentSummaryMetrics']['mem']).split('B')[0],
    log: config['logs_folder'] + "/QC_gatk_CollectAlignmentSummaryMetrics/{prefix}.log"
    threads: config['QC_gatk_CollectAlignmentSummaryMetrics']['threads']
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk --java-options '-Xms{params.mem_gatk} -Xmx{params.mem_gatk} -XX:ParallelGCThreads={threads}' CollectAlignmentSummaryMetrics "
        "-R {input.ref} "
        "-I {input.bam} "
        "-O {output.alignment_metrics} "
        "> {log} 2>&1 "


'''
java -jar picard.jar CollectAlignmentSummaryMetrics \
          R=reference_sequence.fasta \
          I=input.bam \
          O=output.txt

'''