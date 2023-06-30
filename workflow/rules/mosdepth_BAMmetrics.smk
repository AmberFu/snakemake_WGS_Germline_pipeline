## https://github.com/brentp/mosdepth
rule mosdepth_bammetrics:
    input: 
        ref = hg38_reference_genome,
        capture_bed = config['capture_bed'],
        bam = rules.pb_germline_FE.output.outBam,
        bai = rules.pb_germline_FE.output.outBai,
        #bam = config['result_folder']+"/pb_germline_FE/{prefix}.bam",
        #bai = config['result_folder']+"/pb_germline_FE/{prefix}.bam.bai",
        #cram = rules.bam2cram.output.cram,
        #crai = rules.index_cram.output.crai,
    output:
        global_dist = config['result_folder'] + "/mosdepth_bammetrics/{prefix}.mosdepth.global.dist.txt",
        summary = config['result_folder'] + "/mosdepth_bammetrics/{prefix}.mosdepth.summary.txt",
        perbase_bed_gz = temp(config['result_folder'] + "/mosdepth_bammetrics/{prefix}.per-base.bed.gz"),
        regions_bed_gz = temp(config['result_folder'] + "/mosdepth_bammetrics/{prefix}.regions.bed.gz"),
        thresholds_bed_gz = temp(config['result_folder'] + "/mosdepth_bammetrics/{prefix}.thresholds.bed.gz"),
    params:
        out_prefix = config['result_folder'] + "/mosdepth_bammetrics/{prefix}",
        DOCKER = config['mosdepth_bammetrics']['docker'],
        MEM = config['mosdepth_bammetrics']['mem'],
        LSF_LOG = config['logs_folder'] + "/mosdepth_bammetrics/LSF_{prefix}.log",
        JOBNAME = "mosdepth_bammetrics_{prefix}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['mosdepth_bammetrics']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    log: config['logs_folder'] + "/mosdepth_bammetrics/{prefix}.log"
    threads: config['mosdepth_bammetrics']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    conda: "../envs/mosdepth.yaml"
    shell: 
        "mosdepth --threads {threads} "
        "--thresholds 1,8,10,15,20,30,50 --mapq 20 --by {input.capture_bed} "
        "--fasta {input.ref} {params.out_prefix} {input.bam} >> {log} 2>&1"
        

### Unzip:
rule unzip_regions_bed:
    input: rules.mosdepth_bammetrics.output.regions_bed_gz
    output: config['result_folder'] + "/mosdepth_bammetrics/{prefix}.regions.bed"
    params:
        DOCKER = config['unzip_regions_bed']['docker'],
        MEM = config['unzip_regions_bed']['mem'],
        LSF_LOG = config['logs_folder'] + "/unzip_regions_bed/LSF_{prefix}.log",
        JOBNAME = "unzip_regions_bed_{prefix}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['unzip_regions_bed']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    log: config['logs_folder'] + "/unzip_regions_bed/{prefix}.log"
    threads: config['unzip_regions_bed']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    shell: 
        "gunzip {input} > {log} 2>&1"


rule unzip_thresholds_bed:
    input: rules.mosdepth_bammetrics.output.thresholds_bed_gz
    output: config['result_folder'] + "/mosdepth_bammetrics/{prefix}.thresholds.bed"
    params:
        DOCKER = config['unzip_thresholds_bed']['docker'],
        MEM = config['unzip_thresholds_bed']['mem'],
        LSF_LOG = config['logs_folder'] + "/unzip_thresholds_bed/LSF_{prefix}.log",
        JOBNAME = "unzip_thresholds_bed_{prefix}",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['unzip_thresholds_bed']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    log: config['logs_folder'] + "/unzip_thresholds_bed/{prefix}.log"
    threads: config['unzip_thresholds_bed']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    shell: 
        "gunzip {input} > {log} 2>&1"
        
        
### 
PREFIX = list(SAMPLE_DICT.keys())

rule get_info_from_regions_bed:
    input: expand(config['result_folder'] + "/mosdepth_bammetrics/{prefix}.regions.bed", prefix=PREFIX)
    output: config['result_folder'] + "/" +  config['project_name'] + "_all_samples_coverage.tsv",
    params:
        DOCKER = config['get_info_from_regions_bed']['docker'],
        MEM = config['get_info_from_regions_bed']['mem'],
        LSF_LOG = config['logs_folder'] + "/get_info_from_regions_bed/LSF_"+config['project_name']+".log",
        JOBNAME = "get_info_from_regions_bed",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['get_info_from_regions_bed']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    log: config['logs_folder'] + "/get_info_from_regions_bed/"+config['project_name']+".log"
    threads: config['get_info_from_regions_bed']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    run: 
        import pandas as pd
        import os
        info = {}
        for f in input:
            sample_name = os.path.basename(f).replace(".regions.bed","")
            tmp_df = pd.read_csv(f, sep="\t", header=None, names=["chr","stat","end","coverage"], skipfooter=3)
            info[sample_name] = tmp_df["coverage"].mean()
        df = pd.DataFrame(data=info, index=[0])
        df.to_csv(output[0], sep="\t", index=False)


rule get_info_from_thresholds_bed:
    input: expand(config['result_folder'] + "/mosdepth_bammetrics/{prefix}.thresholds.bed", prefix=PREFIX)
    output: config['result_folder'] + "/" +  config['project_name'] + "_all_samples_thresholds.tsv",
    params:
        DOCKER = config['get_info_from_thresholds_bed']['docker'],
        MEM = config['get_info_from_thresholds_bed']['mem'],
        LSF_LOG = config['logs_folder'] + "/get_info_from_thresholds_bed/LSF_"+config['project_name']+".log",
        JOBNAME = "get_info_from_thresholds_bed",
        GPU = "",
        RESOURCE = "-R 'rusage[mem=" + config['get_info_from_thresholds_bed']['mem'] + "] span[hosts=1]'",
        QUEUE = config['general_queue'],
        GROUP = config['general_compute_group'],
    log: config['logs_folder'] + "/get_info_from_thresholds_bed/"+config['project_name']+".log"
    threads: config['get_info_from_thresholds_bed']['threads']
    resources:
        tmpdir=config['ACTIVE_TMP']
    run: 
        import pandas as pd
        import os
        info = []
        names = []
        for f in input:
            sample_name = os.path.basename(f).replace(".thresholds.bed","")
            tmp_df = pd.read_csv(f, sep="\t", header=0, skipfooter=3)
            tmp_df["1X_percent"] = tmp_df.apply(lambda row: row['1X']/row['end'], axis=1)
            tmp_df["8X_percent"] = tmp_df.apply(lambda row: row['8X']/row['end'], axis=1)
            tmp_df["10X_percent"] = tmp_df.apply(lambda row: row['10X']/row['end'], axis=1)
            tmp_df["15X_percent"] = tmp_df.apply(lambda row: row['15X']/row['end'], axis=1)
            tmp_df["20X_percent"] = tmp_df.apply(lambda row: row['20X']/row['end'], axis=1)
            tmp_df["30X_percent"] = tmp_df.apply(lambda row: row['30X']/row['end'], axis=1)
            tmp_df["50X_percent"] = tmp_df.apply(lambda row: row['50X']/row['end'], axis=1)
            mean_df = tmp_df[["1X_percent","8X_percent","10X_percent","15X_percent","20X_percent","30X_percent","50X_percent"]].mean().to_frame()
            mean_df = mean_df.rename(columns={0: sample_name})
            info.append(mean_df)
        df = pd.concat(info, axis=1)
        df.to_csv(output[0], sep="\t")
