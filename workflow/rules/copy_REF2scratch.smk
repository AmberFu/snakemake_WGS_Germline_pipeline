### This SMK is for prepare REF on Scratch space:
rule prepare_hg38_FE_ref:
    input: config['hg38FE']
    output: config['hg38FE_scratch']
    params:
        hg38FE_prefix = config['hg38FE_prefix'],
        scratch_dir = config['hg38FE_scratch_dir'],
    shell: "cp {params.hg38FE_prefix}* {params.scratch_dir}"

rule prepare_dbsnp_vcf:
    input: config['vcf']['storage1']['dbsnp_vcf']
    output: config['vcf']['scratch1']['dbsnp_vcf']
    shell: "cp {input} {output}"

rule prepare_omni_vcf:
    input: config['vcf']['storage1']['omni_vcf']
    output: config['vcf']['scratch1']['omni_vcf']
    shell: "cp {input} {output}"


rule prepare_onekgsnp_vcf:
    input: config['vcf']['storage1']['onekgsnp_vcf']
    output: config['vcf']['scratch1']['onekgsnp_vcf']
    shell: "cp {input} {output}"


rule prepare_hapmap_vcf:
    input: config['vcf']['storage1']['hapmap_vcf']
    output: config['vcf']['scratch1']['hapmap_vcf']
    shell: "cp {input} {output}"


rule prepare_ks19443dbsnp:
    input: config['vcf']['storage1']['ks19443dbsnp']
    output: config['vcf']['scratch1']['ks19443dbsnp']
    shell: "cp {input} {output}"

rule prepare_ks20211mills:
    input: config['vcf']['storage1']['ks20211mills']
    output: config['vcf']['scratch1']['ks20211mills']
    shell: "cp {input} {output}"


rule prepare_ks20267indel:
    input: config['vcf']['storage1']['ks20267indel']
    output: config['vcf']['scratch1']['ks20267indel']
    shell: "cp {input} {output}"