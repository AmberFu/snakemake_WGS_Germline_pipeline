rule bgzip_gvcf:
    input:
        gvcf = rules.pb_germline_FE.output.outVar,
        #gvcf = config['result_folder']+"/pb_germline_FE/{prefix}.g.vcf",
    output:
        bgzip_gvcf = config['result_folder'] + "/bgzip_gvcf/{prefix}.g.vcf.bgz",
    params:
        DOCKER = config['bgzip_gvcf']['docker'],
        MEM = config['bgzip_gvcf']['mem'],
        LSF_LOG = config['logs_folder'] + "/bgzip_gvcf/LSF_{prefix}.log",
        JOBNAME = "bgzip_gvcf_{prefix}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['bgzip_gvcf']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    threads: config['bgzip_gvcf']['threads']
    log: config['logs_folder']+"/bgzip_gvcf/{prefix}.log"
    resources:
        tmpdir=config['ACTIVE_TMP']
    shell: 
        "bgzip -@ {threads} -c {input.gvcf} > {output.bgzip_gvcf} 2> {log}"

rule tabix_bgzip_gvcf:
    input:
        bgzip_gvcf = rules.bgzip_gvcf.output.bgzip_gvcf,
    output:
        tabix_bgzip_gvcf = config['result_folder'] + "/bgzip_gvcf/{prefix}.g.vcf.bgz.tbi",
    params:
        DOCKER = config['tabix_bgzip_gvcf']['docker'],
        MEM = config['tabix_bgzip_gvcf']['mem'],
        LSF_LOG = config['logs_folder'] + "/tabix_bgzip_gvcf/LSF_{prefix}.log",
        JOBNAME = "tabix_bgzip_gvcf_{prefix}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['tabix_bgzip_gvcf']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    threads: config['tabix_bgzip_gvcf']['threads']
    log: config['logs_folder']+"/tabix_bgzip_gvcf/{prefix}.log"
    resources:
        tmpdir=config['ACTIVE_TMP']
    shell: 
        "tabix {input.bgzip_gvcf} 2> {log}"

