rule gatk_indexVCF:
    input:
        vcf = rules.hail_VDS2VCF.output.vcf,
    output:
        vcf_index = config['result_folder'] + "/hail_VDS2VCF/"+config['project_name']+".vcf.idx",
    threads: config['gatk_indexVCF']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    log: config['logs_folder']+"/gatk_indexVCF.log"
    params:
        DOCKER = config['gatk_indexVCF']['docker'],
        MEM = config['gatk_indexVCF']['mem'],
        LSF_LOG = config['logs_folder'] + "/LSF_gatk_indexVCF.log",
        JOBNAME = "gatk_indexVCF",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['gatk_indexVCF']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
        mem_gatk = str(config['QC_gatk_CollectWgsMetrics']['mem']).split('B')[0],
    shell:
        "gatk --java-options '-Xms{params.mem_gatk} -Xmx{params.mem_gatk} -XX:ParallelGCThreads={threads}' IndexFeatureFile "
        "--input {input.vcf}"



rule gatk_VQSR:
    input:
        vcf = rules.hail_VDS2VCF.output.vcf,
        vcf_idx = rules.gatk_indexVCF.output.vcf_index,
    output:
        recal = config['result_folder'] + "/gatk_VQSR/"+config['project_name']+".recal",
        tranches = config['result_folder'] + "/gatk_VQSR/"+config['project_name']+".tranches",
    threads: config['gatk_VQSR']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    log: config['logs_folder']+"/gatk_VQSR.log"
    params:
        DOCKER = config['gatk_VQSR']['docker'],
        MEM = config['gatk_VQSR']['mem'],
        LSF_LOG = config['logs_folder'] + "/LSF_gatk_VQSR.log",
        JOBNAME = "gatk_VQSR",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['gatk_VQSR']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
        mem_gatk = str(config['gatk_VQSR']['mem']).split('B')[0],
        OMNI = config['vcf']['storage1']['omni_vcf'],
        ONEKG = config['vcf']['storage1']['onekgsnp_vcf'],
        HAPMAP = config['vcf']['storage1']['hapmap_vcf'],
        MILLS = config['vcf']['storage1']['mills'],
    shell:
        "gatk --java-options '-Xms{params.mem_gatk} -Xmx{params.mem_gatk} -XX:ParallelGCThreads={threads}' VariantRecalibrator "
        "-V {input.vcf} "
        "-O {output.recal} "
        "--tranches-file {output.tranches} "
        "--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.HAPMAP} "
        "--resource:omni,known=false,training=true,truth=true,prior=12.0 {params.OMNI} "
        "--resource:1000g,known=false,training=true,truth=false,prior=10.0 {params.ONEKG} "
        "--resource:mills,known=false,training=true,truth=true,prior=12.0 {params.MILLS} "
        "-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP -an MQ "
        "-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 "
        "-tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 "
        "-tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 "
        "--mode BOTH "
        
rule gatk_ApplyVQSR:
    input:
        vcf = rules.hail_VDS2VCF.output.vcf,
        vcf_idx = rules.gatk_indexVCF.output.vcf_index,
        recal = rules.gatk_VQSR.output.recal,
        tranches = rules.gatk_VQSR.output.tranches,
    output:
        vcf_VQSR = temp(config['result_folder'] + "/gatk_VQSR/"+config['project_name']+"_VQSR.vcf"),
    threads: config['gatk_ApplyVQSR']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    log: config['logs_folder']+"/gatk_ApplyVQSR.log"
    params:
        DOCKER = config['gatk_ApplyVQSR']['docker'],
        MEM = config['gatk_ApplyVQSR']['mem'],
        LSF_LOG = config['logs_folder'] + "/LSF_gatk_ApplyVQSR.log",
        JOBNAME = "gatk_ApplyVQSR",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['gatk_ApplyVQSR']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
        mem_gatk = str(config['gatk_ApplyVQSR']['mem']).split('B')[0],
    shell:
        "gatk --java-options '-Xms{params.mem_gatk} -Xmx{params.mem_gatk} -XX:ParallelGCThreads={threads}' ApplyVQSR "
        "-V {input.vcf} "
        "--recal-file {input.recal} "
        "--tranches-file {input.tranches} "
        "-O {output.vcf_VQSR} "
        "--mode BOTH "


rule bgz_VQSR_VCF:
    input:
        vcf_VQSR = rules.gatk_ApplyVQSR.output.vcf_VQSR,
    output:
        vcf_VQSR_BGZ = config['result_folder'] + "/gatk_VQSR/"+config['project_name']+"_VQSR.vcf.bgz",
    threads: config['bgz_VQSR_VCF']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    log: config['logs_folder']+"/bgz_VQSR_VCF.log"
    params:
        DOCKER = config['bgz_VQSR_VCF']['docker'],
        MEM = config['bgz_VQSR_VCF']['mem'],
        LSF_LOG = config['logs_folder'] + "/LSF_bgz_VQSR_VCF.log",
        JOBNAME = "bgz_VQSR_VCF",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['bgz_VQSR_VCF']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    shell:
        "bgzip -c {input.vcf_VQSR} > {output.vcf_VQSR_BGZ}"


rule tabix_bgzip_VQSR_VCF:
    input:
        vcf_VQSR_BGZ = rules.bgz_VQSR_VCF.output.vcf_VQSR_BGZ,
    output:
        tabix_vcf_VQSR_BGZ = config['result_folder'] + "/gatk_VQSR/"+config['project_name']+"_VQSR.vcf.bgz.tbi",
    params:
        DOCKER = config['tabix_bgzip_VQSR_VCF']['docker'],
        MEM = config['tabix_bgzip_VQSR_VCF']['mem'],
        LSF_LOG = config['logs_folder'] + "/LSF_tabix_bgzip_VQSR_VCF.log",
        JOBNAME = "tabix_bgzip_VQSR_VCF",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['tabix_bgzip_VQSR_VCF']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    threads: config['tabix_bgzip_VQSR_VCF']['threads']
    log: config['logs_folder']+"/tabix_bgzip_VQSR_VCF.log"
    resources:
        tmpdir=config['ACTIVE_TMP']
    shell: 
        "tabix {input.vcf_VQSR_BGZ} 2> {log}"


'''
### Index VCF:
gatk --java-options "-Xmx50g -Xms50g" IndexFeatureFile --input $OUTPUT_CONCAT_VCF


### RUN VQSR
# Follow gnomAD v3.1 post: https://gnomad.broadinstitute.org/news/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/#variant-qc
# Follow GATK4 post - https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering#1
recal_outpath=$out_vcf_dir/$out_vcf_basename".recal"
tranches_outpath=$out_vcf_dir/$out_vcf_basename".tranches"
resource_omni="/storage1/fs1/jin810/Active/known_sites/1000G_omni2.5.hg38.vcf"
resource_1kg="/storage1/fs1/jin810/Active/known_sites/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
resource_hapmap="/storage1/fs1/jin810/Active/known_sites/hapmap_3.3.hg38.vcf.gz"
resource_dbsnp="/storage1/fs1/jin810/Active/known_sites/resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf"
resource_mills="/storage1/fs1/jin810/Active/known_sites/resources-broad-hg38-v0-Mills_and_1000G_gold_standard.indels.hg38.vcf"
resource_axiomPoly="/storage1/fs1/jin810/Active/known_sites/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"

time /gatk/gatk --java-options "-Xmx50g -Xms50g" VariantRecalibrator \
    -V $in_vcf \
    -O $recal_outpath \
    --tranches-file $tranches_outpath \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $resource_hapmap \
    --resource:omni,known=false,training=true,truth=true,prior=12.0 $resource_omni \
    --resource:1000g,known=false,training=true,truth=false,prior=10.0 $resource_1kg \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 $resource_mills \
    -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP -an MQ \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 \
    -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 \
    -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    --mode BOTH
'''