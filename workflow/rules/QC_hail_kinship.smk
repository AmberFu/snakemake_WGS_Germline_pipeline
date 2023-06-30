rule hail_vqsrVcf2mt:
    input:
        VQSR_VCF_BGZ = rules.bgz_VQSR_VCF.output.vcf_VQSR_BGZ,
        VQSR_VCF_BGZ_tabix = rules.tabix_bgzip_VQSR_VCF.output.tabix_vcf_VQSR_BGZ,
    output:
        mt = config['result_folder'] + "/" + config['project_name']+"_VQSR.mt",
    threads: config['hail_vqsrVcf2mt']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    log: config['logs_folder']+"/hail_vqsrVcf2mt.log"
    params:
        DOCKER = config['hail_vqsrVcf2mt']['docker'],
        MEM = config['hail_vqsrVcf2mt']['mem'],
        LSF_LOG = config['logs_folder'] + "/LSF_hail_vqsrVcf2mt.log",
        JOBNAME = "hail_vqsrVcf2mt",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['hail_vqsrVcf2mt']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    script:
        "../scripts/hail_VCF2MT.py"


rule hail_kinship:
    input:
        mt = rules.hail_vqsrVcf2mt.output.mt,
    output:
        #scores_table = config['result_folder'] + "/hail_kinship/"+config['project_name']+"_scores_table.ht",
        pc_rel = config['result_folder'] + "/hail_kinship/"+config['project_name']+"_pc_rel.ht",
        related_samples_to_remove = config['result_folder'] + "/hail_kinship/"+config['project_name']+"_related_samples_to_remove.tsv",
    threads: config['hail_kinship']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    log: config['logs_folder']+"/hail_kinship.log"
    params:
        DOCKER = config['hail_kinship']['docker'],
        MEM = config['hail_kinship']['mem'],
        LSF_LOG = config['logs_folder'] + "/LSF_hail_kinship.log",
        JOBNAME = "hail_kinship",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['hail_kinship']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
        maf = 0.001,
        k = 2,
        kin_cutoff = 0.125,
    script:
        "../scripts/hail_Kinship_analysis.py"
