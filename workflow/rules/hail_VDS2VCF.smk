rule hail_VDS2VCF:
    input:
        vds = rules.hail_joint_calling.output.vds,
    output:
        denseMT = temp(directory(config['result_folder'] + "/hail_VDS2VCF/"+config['project_name']+"_dense.mt")),
        vcf = config['result_folder'] + "/hail_VDS2VCF/"+config['project_name']+".vcf",
    threads: config['hail_VDS2VCF']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    log: config['logs_folder']+"/hail_VDS2VCF.log"
    params:
        DOCKER = config['hail_VDS2VCF']['docker'],
        MEM = config['hail_VDS2VCF']['mem'],
        LSF_LOG = config['logs_folder'] + "/LSF_hail_VDS2VCF.log",
        JOBNAME = "hail_VDS2VCF",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['hail_VDS2VCF']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    script:
        "../scripts/hail_vds2vcf.py"