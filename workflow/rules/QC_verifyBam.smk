# check contamination:
rule verify_bam:
    input:
        bam = rules.pb_germline_FE.output.outBam,
        bai = rules.pb_germline_FE.output.outBai,
        vcf = ONEKGENOME_SNP,
    output:
        bestSM = config['result_folder'] + "/verify_bam/{prefix}_verify_bam.bestSM",
        depthSM = config['result_folder'] + "/verify_bam/{prefix}_verify_bam.depthSM",
        selfSM = config['result_folder'] + "/verify_bam/{prefix}_verify_bam.selfSM"
    threads: config['verify_bam']['threads']
    log: config['logs_folder']+"/verify_bam/{prefix}.log"
    params:
        DOCKER = config['verify_bam']['docker'],
        MEM = config['verify_bam']['mem'],
        LSF_LOG = config['logs_folder'] + "/verify_bam/LSF_{prefix}.log",
        JOBNAME = "verify_bam_{filename}",
        prefix = config['result_folder'] + "/verify_bam/{prefix}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['verify_bam']['mem'] + "] span[hosts=1]'",
    shell: 
        "verifyBamID --NumThread {threads} --vcf {input.vcf} --bam {input.bam} "
        "--best --ignoreRG --out {params.prefix} 2>&1 {log}"


rule verify_bam2:
    input:
        bam = rules.pb_germline_FE.output.outBam,
        bai = rules.pb_germline_FE.output.outBai,
        ref = hg38_reference_genome,
    output:
        selfSM = config['result_folder'] + "/verify_bam2/{prefix}.selfSM",
        Ancestry = config['result_folder'] + "/verify_bam2/{prefix}.Ancestry",
    threads: config['verify_bam2']['threads']
    log: config['logs_folder']+"/verify_bam2/{prefix}.log"
    params:
        DOCKER = config['verify_bam2']['docker'],
        MEM = config['verify_bam2']['mem'],
        LSF_LOG = config['logs_folder'] + "/verify_bam2/LSF_{prefix}.log",
        JOBNAME = "verify_bam2_{prefix}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['verify_bam2']['mem'] + "] span[hosts=1]'",
        prefix = config['result_folder'] + "/verify_bam2/{prefix}",
        SVDPrefix = VERIFYBAMID2_SVDPREFIX,
    shell: 
        "verifybamid2 --Reference {input.ref} --BamFile {input.bam} --NumThread {threads} "
        "--SVDPrefix {params.SVDPrefix} --Output {params.prefix} 2>&1 {log}"
