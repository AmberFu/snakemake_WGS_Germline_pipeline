# check contamination:
rule verify_bam2:
    input:
        #bam = config['result_folder']+"/pb_germline_FE/{prefix}.bam",
        #bai = config['result_folder']+"/pb_germline_FE/{prefix}.bam.bai",
        bam = rules.pb_germline_FE.output.outBam,
        bai = rules.pb_germline_FE.output.outBai,
        ref = hg38_reference_genome,
    output:
        selfSM = config['result_folder'] + "/verify_bam2/{prefix}.selfSM",
        Ancestry = config['result_folder'] + "/verify_bam2/{prefix}.Ancestry",
    threads: config['verify_bam2']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    log: config['logs_folder']+"/verify_bam2/{prefix}.log"
    params:
        DOCKER = config['verify_bam2']['docker'],
        MEM = config['verify_bam2']['mem'],
        LSF_LOG = config['logs_folder'] + "/verify_bam2/LSF_{prefix}.log",
        JOBNAME = "verify_bam2_{prefix}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['verify_bam2']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
        prefix = config['result_folder'] + "/verify_bam2/{prefix}",
        SVDPrefix = VERIFYBAMID2_SVDPREFIX,
    shell: 
        "verifybamid2 --Reference {input.ref} --BamFile {input.bam} --NumThread {threads} "
        "--SVDPrefix {params.SVDPrefix} --Output {params.prefix} 2>&1"
