
rule hail_joint_calling:
    input:
        bgz_gvcfs = expand(config['result_folder'] + "/bgzip_gvcf/{prefix}.g.vcf.bgz", prefix=PREFIX),
        bgz_gvcfs_index = expand(config['result_folder'] + "/bgzip_gvcf/{prefix}.g.vcf.bgz.tbi", prefix=PREFIX),
    output:
        vds = directory(config['result_folder'] + "/hail_joint_calling/"+config['project_name']+".vds"),
    threads: config['hail_joint_calling']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    log: config['logs_folder']+"/hail_joint_calling.log"
    params:
        DOCKER = config['hail_joint_calling']['docker'],
        MEM = config['hail_joint_calling']['mem'],
        LSF_LOG = config['logs_folder'] + "/LSF_hail_joint_calling.log",
        JOBNAME = "hail_joint_calling",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['hail_joint_calling']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    script:
        "../scripts/hail_jointcalling_VDS.py"