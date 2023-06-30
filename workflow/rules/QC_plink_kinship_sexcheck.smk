#############
## KINSHIP
#############
rule plink_kinship:
    input:
        vqsr_vcf = rules.bgz_VQSR_VCF.output.vcf_VQSR_BGZ,
        tabix_vqsr_vcf = rules.tabix_bgzip_VQSR_VCF.output.tabix_vcf_VQSR_BGZ,
    output:
        ibd_report = config['result_folder'] + "/plink_kinship/" + config['project_name'] + "_kinship.txt",
    params:
        DOCKER = config['plink_kinship']['docker'],
        MEM = config['plink_kinship']['mem'],
        LSF_LOG = config['logs_folder'] + "/LSF_plink_kinship.log",
        JOBNAME = "plink_kinship",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['plink_kinship']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    threads: config['plink_kinship']['threads']
    log: config['logs_folder']+"/plink_kinship.log"
    resources:
        tmpdir=config['ACTIVE_TMP']
    shell: 
        "/bin/plink --geno 0.01 --genome --hwe 0.001 --maf 0.05 "
        "--snps-only --allow-extra-chr --vcf-half-call m "
        "--double-id --vcf {input.vqsr_vcf} "
        "--out {output.ibd_report} "

#############
## SEX CHECK
#############
## Note: VCF is quicker than VCF.bgz!
rule plink_sexcheck1:
    input:
        vqsr_vcf = rules.gatk_ApplyVQSR.output.vcf_VQSR,
    output:
        sexcheck1_bed = temp(config['result_folder'] + "/plink_sexcheck/" + config['project_name'] + "_plink_sexcheck1.bed"),
        sexcheck1_bim = temp(config['result_folder'] + "/plink_sexcheck/" + config['project_name'] + "_plink_sexcheck1.bim"),
        sexcheck1_fam = temp(config['result_folder'] + "/plink_sexcheck/" + config['project_name'] + "_plink_sexcheck1.fam"),
        sexcheck1_log = config['result_folder'] + "/plink_sexcheck/" + config['project_name'] + "_plink_sexcheck1.log",
    threads: config['plink_sexcheck1']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    log: config['logs_folder']+"/plink_sexcheck1.log"
    params:
        prefix = config['project_name'] + "_plink_sexcheck1",
        DOCKER = config['plink_sexcheck1']['docker'],
        MEM = config['plink_sexcheck1']['mem'],
        LSF_LOG = config['logs_folder'] + "/LSF_plink_sexcheck1.log",
        JOBNAME = "plink_sexcheck1",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['plink_sexcheck1']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    shell:
        "/bin/plink --double-id --allow-extra-chr --vcf-half-call m "
        "--vcf {input.vqsr_vcf} "
        "--out {params.prefix} > {log} 2>&1"


rule plink_sexcheck2:
    input:
        sexcheck1_prefix = rules.plink_sexcheck1.params.prefix,
    output:
        sexcheck2_bed = temp(config['result_folder'] + "/plink_sexcheck/" + config['project_name'] + "_plink_sexcheck2.bed"),
        sexcheck2_bim = temp(config['result_folder'] + "/plink_sexcheck/" + config['project_name'] + "_plink_sexcheck2.bim"),
        sexcheck2_fam = temp(config['result_folder'] + "/plink_sexcheck/" + config['project_name'] + "_plink_sexcheck2.fam"),
        sexcheck2_nosex = temp(config['result_folder'] + "/plink_sexcheck/" + config['project_name'] + "_plink_sexcheck2.nosex"),
        sexcheck2_log = config['result_folder'] + "/plink_sexcheck/" + config['project_name'] + "_plink_sexcheck2.log",
    threads: config['plink_sexcheck2']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    log: config['logs_folder']+"/plink_sexcheck2.log"
    params:
        prefix = config['project_name'] + "_plink_sexcheck2",
        DOCKER = config['plink_sexcheck2']['docker'],
        MEM = config['plink_sexcheck2']['mem'],
        LSF_LOG = config['logs_folder'] + "/LSF_plink_sexcheck2.log",
        JOBNAME = "plink_sexcheck2",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['plink_sexcheck2']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    shell:
        "/bin/plink -split-x hg38 --make-bed --allow-extra-chr "
        "--bfile {input.sexcheck1_prefix} "
        "--out {params.prefix} > {log} 2>&1"


rule plink_sexcheck3:
    input:
        sexcheck2_prefix = rules.plink_sexcheck2.params.prefix,
    output:
        sexcheck3_bed = temp(config['result_folder'] + "/plink_sexcheck/" + config['project_name'] + "_plink_sexcheck3.bed"),
        sexcheck3_bim = temp(config['result_folder'] + "/plink_sexcheck/" + config['project_name'] + "_plink_sexcheck3.bim"),
        sexcheck3_fam = temp(config['result_folder'] + "/plink_sexcheck/" + config['project_name'] + "_plink_sexcheck3.fam"),
        sexcheck3_nosex = temp(config['result_folder'] + "/plink_sexcheck/" + config['project_name'] + "_plink_sexcheck3.nosex"),
        sexcheck3_sexcheck = config['result_folder'] + "/plink_sexcheck/" + config['project_name'] + "_plink_sexcheck3.sexcheck",
        sexcheck3_log = config['result_folder'] + "/plink_sexcheck/" + config['project_name'] + "_plink_sexcheck3.log",
    threads: config['plink_sexcheck3']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    log: config['logs_folder']+"/plink_sexcheck3.log"
    params:
        prefix = config['project_name'] + "_plink_sexcheck3",
        DOCKER = config['plink_sexcheck3']['docker'],
        MEM = config['plink_sexcheck3']['mem'],
        LSF_LOG = config['logs_folder'] + "/LSF_plink_sexcheck3.log",
        JOBNAME = "plink_sexcheck3",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['plink_sexcheck3']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    shell:
        "/bin/plink --impute-sex --make-bed --allow-extra-chr "
        "--bfile {input.sexcheck2_prefix} "
        "--out {params.prefix} > {log} 2>&1"



'''
DOCKERIMG_1="broadinstitute/gatk:latest"    ##-----gatk select variants
DOCKERIMG_2_1="spashleyfu/hail_0.2.78:latest"    ##-----hail joint
DOCKERIMG_2_2="spashleyfu/hail_vep_gnomad:latest"    ##-----hail mt to vcf
#DOCKERIMG_3="biocontainers/bcftools:v1.5_cv3"   ##-----bcf tool-bcf2vcf
DOCKERIMG_4="sam16711/plink:latest"   ##-----Kindship Sexcheck
DOCKERIMG_5="spashleyfu/ubuntu18_vep104:hail_gsutil"   ##-----PCA liftover
DOCKERIMG_6="dreammaerd/laser_trace_v2.04:latest"     ##-----Laser PCA
DOCKERIMG_7="dreammaerd/r-jupyter:4.2.1"    ##-----R script

# kinship Commend
kinshipCmd="/bin/plink --geno 0.01 --genome --hwe 0.001 --maf 0.05 "
kinshipCmd+="--snps-only --allow-extra-chr --vcf-half-call m "
kinshipCmd+="--vcf $VCF_LOC_k "
kinshipCmd+="--out $Kin_out"

# Sexcheck Commend
sexcheck1_out=$Sexcheck_OUT/"$batch_name"_all_s1

sexcheck1="/bin/plink --double-id --allow-extra-chr --vcf-half-call m "
sexcheck1+="--vcf $VCF_LOC_s "
sexcheck1+="--out $sexcheck1_out"

# Step 2-split X
sexcheck2_out=$Sexcheck_OUT/"$batch_name"_all_s2

sexcheck2="/bin/plink -split-x hg38 --make-bed --allow-extra-chr "
sexcheck2+="--bfile $sexcheck1_out "
sexcheck2+="--out $sexcheck2_out"

# Step 3-sex check
sexcheck3_out=$Sexcheck_OUT/"$batch_name"_all_s3

sexcheck3="/bin/plink --impute-sex --make-bed --allow-extra-chr "
sexcheck3+="--bfile $sexcheck2_out "
sexcheck3+="--out $sexcheck3_out"



### Actual running: VCF is quicker than VCF.bgz
fup@compute1-exec-130:/storage1/fs1/jin810/Active/fup/UDN_CP_download_for_Nahyun$ ls -lh /storage1/fs1/jin810/Active/fup/UDN_CP_download_for_Nahyun/WGS_50PilotTrios_VQSR.vcf
-rw-------. 1 fup domain users 104G Jan 23 03:35 

fup@compute1-exec-130:/storage1/fs1/jin810/Active/fup/UDN_CP_download_for_Nahyun$ time /bin/plink  --double-id --allow-extra-chr --vcf-half-call m --vcf $VCF_LOC_s --out WGS_50PilotTrios_VQSR_plink_sexcheck1
PLINK v1.90b6.18 64-bit (16 Jun 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to WGS_50PilotTrios_VQSR_plink_sexcheck1.log.
Options in effect:
  --allow-extra-chr
  --double-id
  --out WGS_50PilotTrios_VQSR_plink_sexcheck1
  --vcf /storage1/fs1/jin810/Active/fup/UDN_CP_download_for_Nahyun/WGS_50PilotTrios_VQSR.vcf.bgz
  --vcf-half-call m

385377 MB RAM detected; reserving 192688 MB for main workspace.
--vcf: WGS_50PilotTrios_VQSR_plink_sexcheck1.bed +
WGS_50PilotTrios_VQSR_plink_sexcheck1.bim +
WGS_50PilotTrios_VQSR_plink_sexcheck1.fam written.

real	7m29.035s
user	7m21.776s
sys	0m5.088s

fup@compute1-exec-130:/storage1/fs1/jin810/Active/fup/UDN_CP_download_for_Nahyun$ VCF_LOC_s=$PWD"/WGS_50PilotTrios_VQSR.vcf"
fup@compute1-exec-130:/storage1/fs1/jin810/Active/fup/UDN_CP_download_for_Nahyun$ time /bin/plink  --double-id --allow-extra-chr --vcf-half-call m --vcf $VCF_LOC_s --out WGS_50PilotTrios_VQSR_plink_sexcheck1.txt
PLINK v1.90b6.18 64-bit (16 Jun 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to WGS_50PilotTrios_VQSR_plink_sexcheck1.txt.log.
Options in effect:
  --allow-extra-chr
  --double-id
  --out WGS_50PilotTrios_VQSR_plink_sexcheck1.txt
  --vcf /storage1/fs1/jin810/Active/fup/UDN_CP_download_for_Nahyun/WGS_50PilotTrios_VQSR.vcf
  --vcf-half-call m

385377 MB RAM detected; reserving 192688 MB for main workspace.
--vcf: WGS_50PilotTrios_VQSR_plink_sexcheck1.txt.bed +
WGS_50PilotTrios_VQSR_plink_sexcheck1.txt.bim +
WGS_50PilotTrios_VQSR_plink_sexcheck1.txt.fam written.

real	1m29.683s
user	1m10.682s
sys	0m17.794s


fup@compute1-exec-130:/storage1/fs1/jin810/Active/fup/UDN_CP_download_for_Nahyun$ ls -lh WGS_50PilotTrios_VQSR_plink_sexcheck1.txt*   
-rw-------. 1 fup domain users 1.4G Feb  6 23:50 WGS_50PilotTrios_VQSR_plink_sexcheck1.txt.bed
-rw-------. 1 fup domain users 749M Feb  6 23:50 WGS_50PilotTrios_VQSR_plink_sexcheck1.txt.bim
-rw-------. 1 fup domain users 9.1K Feb  6 23:48 WGS_50PilotTrios_VQSR_plink_sexcheck1.txt.fam
-rw-------. 1 fup domain users  702 Feb  6 23:50 WGS_50PilotTrios_VQSR_plink_sexcheck1.txt.log

fup@compute1-exec-130:/storage1/fs1/jin810/Active/fup/UDN_CP_download_for_Nahyun$  time /bin/plink -split-x hg38 --make-bed --allow-extra-chr --bfile WGS_50PilotTrios_VQSR_plink_sexcheck1 --out WGS_50PilotTrios_VQSR_plink_sexcheck2
PLINK v1.90b6.18 64-bit (16 Jun 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to WGS_50PilotTrios_VQSR_plink_sexcheck2.log.
Options in effect:
  --allow-extra-chr
  --bfile WGS_50PilotTrios_VQSR_plink_sexcheck1
  --make-bed
  --out WGS_50PilotTrios_VQSR_plink_sexcheck2
  --split-x hg38

385377 MB RAM detected; reserving 192688 MB for main workspace.
31995139 variants loaded from .bim file.
177 people (0 males, 0 females, 177 ambiguous) loaded from .fam.
Ambiguous sex IDs written to WGS_50PilotTrios_VQSR_plink_sexcheck2.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 177 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate is 0.445533.
31995139 variants and 177 people pass filters and QC.
Note: No phenotypes present.
--split-x: 63277 chromosome codes changed.
--make-bed to WGS_50PilotTrios_VQSR_plink_sexcheck2.bed +
WGS_50PilotTrios_VQSR_plink_sexcheck2.bim +
WGS_50PilotTrios_VQSR_plink_sexcheck2.fam ... done.

real	0m16.673s
user	0m13.715s
sys	0m2.287s


fup@compute1-exec-130:/storage1/fs1/jin810/Active/fup/UDN_CP_download_for_Nahyun$ time /bin/plink --impute-sex --make-bed --allow-extra-chr --bfile WGS_50PilotTrios_VQSR_plink_sexcheck2 --out WGS_50PilotTrios_VQSR_plink_sexcheck3
PLINK v1.90b6.18 64-bit (16 Jun 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to WGS_50PilotTrios_VQSR_plink_sexcheck3.log.
Options in effect:
  --allow-extra-chr
  --bfile WGS_50PilotTrios_VQSR_plink_sexcheck2
  --impute-sex
  --make-bed
  --out WGS_50PilotTrios_VQSR_plink_sexcheck3

385377 MB RAM detected; reserving 192688 MB for main workspace.
31995139 variants loaded from .bim file.
177 people (0 males, 0 females, 177 ambiguous) loaded from .fam.
Ambiguous sex IDs written to WGS_50PilotTrios_VQSR_plink_sexcheck3.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 177 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate is 0.445533.
31995139 variants and 177 people pass filters and QC.
Note: No phenotypes present.
--impute-sex: 1161727 Xchr and 0 Ychr variant(s) scanned, 107/177 sexes
imputed. Report written to WGS_50PilotTrios_VQSR_plink_sexcheck3.sexcheck .
--make-bed to WGS_50PilotTrios_VQSR_plink_sexcheck3.bed +
WGS_50PilotTrios_VQSR_plink_sexcheck3.bim +
WGS_50PilotTrios_VQSR_plink_sexcheck3.fam ... done.

real	0m13.933s
user	0m11.363s
sys	0m2.191s


fup@compute1-exec-130:/storage1/fs1/jin810/Active/fup/UDN_CP_download_for_Nahyun$ time /bin/plink --geno 0.01 --genome --hwe 0.001 --maf 0.05 --double-id --snps-only --allow-extra-chr --vcf-half-call m --vcf WGS_50PilotTrios_VQSR.vcf --out WGS_50PilotTrios_VQSR_plink_kinship
PLINK v1.90b6.18 64-bit (16 Jun 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to WGS_50PilotTrios_VQSR_plink_kinship.log.
Options in effect:
  --allow-extra-chr
  --double-id
  --geno 0.01
  --genome
  --hwe 0.001
  --maf 0.05
  --out WGS_50PilotTrios_VQSR_plink_kinship
  --snps-only
  --vcf WGS_50PilotTrios_VQSR.vcf
  --vcf-half-call m

385377 MB RAM detected; reserving 192688 MB for main workspace.
--vcf: WGS_50PilotTrios_VQSR_plink_kinship-temporary.bed +
WGS_50PilotTrios_VQSR_plink_kinship-temporary.bim +
WGS_50PilotTrios_VQSR_plink_kinship-temporary.fam written.
22579944 out of 31995139 variants loaded from .bim file.
177 people (0 males, 0 females, 177 ambiguous) loaded from .fam.
Ambiguous sex IDs written to WGS_50PilotTrios_VQSR_plink_kinship.nosex .
Using up to 63 threads (change this with --threads).
Before main variant filters, 177 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate is 0.341042.
22059016 variants removed due to missing genotype data (--geno).
--hwe: 66667 variants removed due to Hardy-Weinberg exact test.
313641 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
140620 variants and 177 people pass filters and QC.
Note: No phenotypes present.
Excluding 672 variants on non-autosomes from IBD calculation.
IBD calculations complete.  
Finished writing WGS_50PilotTrios_VQSR_plink_kinship.genome .

real	1m39.471s
user	1m18.456s
sys	0m19.978s
'''