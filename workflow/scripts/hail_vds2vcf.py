'''
https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#python

def do_something(data_path, out_path, threads, myparam):
    # python code

do_something(snakemake.input[0], snakemake.output[0], snakemake.threads, snakemake.config["myparam"])
'''

import argparse
import os
import math
import subprocess
import time
from datetime import datetime
import pyspark
import hail as hl
from gnomad.utils.annotations import (
    get_adj_expr,
    get_lowqual_expr,
    bi_allelic_site_inbreeding_expr,
)
from gnomad.utils.sparse_mt import (
    get_site_info_expr,
    split_info_annotation,
    # split_lowqual_annotation,
)
from gnomad.utils.file_utils import file_exists
from gnomad.utils.vcf import adjust_vcf_incompatible_types

########################################
today = datetime.now().strftime("%Y-%m-%d")
print("Program Starting Date: {}".format(today))

# HAIL_DIR:
# HAIL_HOME=$(pip3 show hail | grep Location | awk -F' ' '{print $2 "/hail"}')
# https://hail.is/docs/0.2/install/other-cluster.html
HAIL_HOME = subprocess.getoutput("pip3 show hail | grep Location | awk -F' ' '{print $2 \"/hail\"}'")
os.environ["HAIL_DIR"] = HAIL_HOME + "/backend"
hail_jars = HAIL_HOME + "/backend/hail-all-spark.jar"
print("HAIL_DIR = {}/backend".format(HAIL_HOME))
# JAVA_HOME:
# for spashleyfu/ubuntu18_vep104:hail_gsutil
# os.environ["JAVA_HOME"] = "/opt/conda"
# for spashleyfu/hail_0.2.78:latest and spashleyfu/hail_vep_gnomad:latest
os.environ["JAVA_HOME"] = "/usr/lib/jvm/java-8-openjdk-amd64/jre"
print("JAVA_HOME = /usr/lib/jvm/java-8-openjdk-amd64/jre")

# Allocate resources:
threads = snakemake.threads
memory = snakemake.params.MEM.lower().split('b')[0] # snakemake.params.MEM is like "24GB"
print("Run with {} CPUs, and spark.driver.memory={}, spark.driver.maxResultSize={}".format(threads, memory, memory))

# CONFIG:
conf = pyspark.SparkConf().setAll([
    ('spark.master', 'local[{}]'.format(threads)),
    ('spark.app.name', 'Hail'),
    ('spark.jars', str(hail_jars)),
    ('spark.driver.extraClassPath', str(hail_jars)),
    ('spark.executor.extraClassPath', './hail-all-spark.jar'),
    ('spark.serializer', 'org.apache.spark.serializer.KryoSerializer'),
    ('spark.kryo.registrator', 'is.hail.kryo.HailKryoRegistrator'),
    ### https://discuss.hail.is/t/turning-run-combiner-performance-for-hail-local-mode/2318
    ('spark.driver.memory', str(memory)),
    ('spark.executor.memory', str(memory)),
    ('spark.driver.maxResultSize', str(memory)),
    ### For IOException: No space left on device
    # https://spark.apache.org/docs/2.3.0/configuration.html
    ('spark.local.dir', '/storage1/fs1/jin810/Active/fup/hail_tmp'),
    ])
### Using sc:
sc = pyspark.SparkContext(conf=conf)
project_name = snakemake.config["project_name"]
logfile = "/storage1/fs1/jin810/Active/fup/hailLogs/hail_VDS2VCF_{}_{}.log".format(project_name, today)
### Hail Start...
hl.init(default_reference='GRCh38',sc=sc,log=logfile)

####################################

INPUT_VDS = snakemake.input.vds
OUTPUT_denseMT = snakemake.output.denseMT
OUTPUT_VCF = snakemake.output.vcf

## 1-1. Read-in VDS
vds = hl.vds.read_vds(INPUT_VDS)

## 1-2. VDS to Merged Sparse MT
mt = hl.vds.to_merged_sparse_mt(vds)

# 2. Keyed by locus and alleles
mt = mt.key_rows_by(mt.locus, mt.alleles)
# 3. Filter out no variants positions:
mt = mt.filter_rows((hl.len(mt.alleles) > 1))
# 4. Annotate "alt_alleles_range_array"
mt = mt.annotate_rows(alt_alleles_range_array=hl.range(1, hl.len(mt.alleles)))
# 5. Annotate entries -  QUALapprox by Sum of PL[0] values!
mt = mt.annotate_entries(QUALapprox=
                         hl.int32(hl.if_else(hl.is_missing(mt.LPL),hl.null('int32'),mt.LPL[0])))
# 6. Annotate entries - VarDP as Depth over variant genotypes (does not include depth of reference samples)
mt = mt.annotate_entries(VarDP=
                         hl.int32(hl.if_else(mt.LGT.is_non_ref(), mt.DP, 0)))

# 7. Compute site level INFO expr
# https://github.com/broadinstitute/gnomad_qc/blob/e9e28ff69632e8a7965708dbc4d9ce7d89999fcf/gnomad_qc/v3/annotations/generate_qc_annotations.py#L47
info_expr = get_site_info_expr(
    mt,
    sum_agg_fields=["QUALapprox"],
    int32_sum_agg_fields=["VarDP"],
    array_sum_agg_fields=["SB","RAW_MQandDP"],
    median_agg_fields=[]
)

# 8. Add AC and AC_raw: (remove AC_raw at this point)
# First compute ACs for each non-ref allele, grouped by adj
# https://github.com/broadinstitute/gnomad_methods/blob/978b57762108d6ccbe979ed7a867ffe361d146e6/gnomad/utils/annotations.py#L701
grp_ac_expr = hl.agg.array_agg(
    lambda ai: hl.agg.filter(
        mt.LA.contains(ai),
        hl.agg.group_by(
            get_adj_expr(mt.LGT, mt.GQ, mt.DP, mt.LAD),
            hl.agg.sum(
                mt.LGT.one_hot_alleles(mt.LA.map(lambda x: hl.str(x)))[
                    mt.LA.index(ai)
                ]
            ),
        ),
    ),
    mt.alt_alleles_range_array,
)
# Then, for each non-ref allele, compute AC as the adj group
# AC_raw as the sum of adj and non-adj groups
info_expr = info_expr.annotate(
    AC_raw=grp_ac_expr.map(lambda i: hl.int32(i.get(True, 0) + i.get(False, 0))),
    AC=grp_ac_expr.map(lambda i: hl.int32(i.get(True, 0))),
)

# 9. Annotate "AN": 
# Total number of alleles in called genotypes (AN)
# https://gatk.broadinstitute.org/hc/en-us/articles/4409917348891-ChromosomeCounts
AN_from_MT = mt.cols().count() * 2
info_expr = info_expr.annotate(
    AN = AN_from_MT
)

# 10. Annotate "AF" = AC/AN
# Frequency of each ALT allele, in the same order as listed (AF)
# info_expr = info_expr.annotate(
#     AF = info_expr.AC.map(lambda i: i/AN_from_MT)
# )

# 11. Create INFO_HT:
info_ht = mt.select_rows(info=info_expr).rows()

# 12. Add lowqual flag
info_ht = info_ht.annotate(
    LowQual=get_lowqual_expr(
        info_ht.alleles,
        info_ht.info.QUALapprox,
        # The indel het prior used for gnomad v3 was 1/10k bases (phred=40).
        # This value is usually 1/8k bases (phred=39).
        indel_phred_het_prior=40,
    )
)

# 13. Adjust Types:
# https://github.com/broadinstitute/gnomad_methods/blob/c673ea3b3829f4d2ee19a6112b221d752aa22d92/gnomad/utils/vcf.py#L326
vcf_ht = adjust_vcf_incompatible_types(info_ht)

# 14. Split Sparse MT's Multi-alleles:
"""
hail.experimental.sparse_split_multi(sparse_mt, *, filter_changed_loci=False)
    The split multi logic handles the following entry fields:
    struct {
        LGT: call
        LAD: array<int32>
        DP: int32
        GQ: int32
        LPL: array<int32>
        RGQ: int32
        LPGT: call
        LA: array<int32>
        END: int32
    }
"""
split_mt = hl.experimental.sparse_split_multi(mt)

# 15. Split INFO_HT's Multi-alleles:
############# SPLIT MULTI-ALLELES Functions start...
def split_info_annotation_fup(
    info_expr: hl.expr.StructExpression, a_index: hl.expr.Int32Expression
) -> hl.expr.StructExpression:
    """
    Split multi-allelic allele-specific info fields.
    :param info_expr: Field containing info struct.
    :param a_index: Allele index. Output by hl.split_multi or hl.split_multi_hts.
    :return: Info struct with split annotations.
    * Modified from https://github.com/broadinstitute/gnomad_methods/blob/978b57762108d6ccbe979ed7a867ffe361d146e6/gnomad/utils/sparse_mt.py#L510
    """
    # Index AS annotations
    info_expr = info_expr.annotate(
        **{
            f: info_expr[f][a_index - 1]
            for f in info_expr
            if f.startswith("AC") or (f.startswith("AS_") and not f == "AS_SB_TABLE")
        },
        # AS_SB_TABLE=info_expr.AS_SB_TABLE[0].extend(info_expr.AS_SB_TABLE[a_index]),
    )
    return info_expr

def split_info_fup(info_ht: hl.Table) -> hl.Table:
    """
    Generates an info table that splits multi-allelic sites from the multi-allelic info table.
    :return: Info table with split multi-allelics
    :rtype: Table
    * Modified from https://github.com/broadinstitute/gnomad_qc/blob/8b2342133e6186ce7056820ca5ec1b66032173bd/gnomad_qc/v3/annotations/generate_qc_annotations.py#L141
    """
    # Create split version
    info_ht = hl.split_multi(info_ht)
    info_ht = info_ht.annotate(
        info=info_ht.info.annotate(
            **split_info_annotation_fup(info_ht.info, info_ht.a_index),
        ),
        # AS_lowqual=split_lowqual_annotation(info_ht.AS_lowqual, info_ht.a_index),
    )
    return info_ht


############# SPLIT MULTI-ALLELES Functions end!
vcf_ht = split_info_fup(vcf_ht)

# 16. Annotate splited INFO_HT onto MT:
split_mt = split_mt.annotate_rows(info = vcf_ht[split_mt.locus, split_mt.alleles].info)

# 17. Annotate "DP": Sum-up
split_mt = split_mt.annotate_rows(info=split_mt.info.annotate(
    DP = hl.int32(hl.agg.sum(split_mt.DP))
))

# 10. Annotate "AF": Moved here...calculate after split
split_mt = split_mt.annotate_rows(info=split_mt.info.annotate(
    AF = hl.float64(split_mt.info.AC_raw/split_mt.info.AN)
))

##### RankSum:
# Monkol suggest: 
# A. take median
# B. take the min and max, and then take the value with the highest absolute value 
#    (compare their unsigned values, then take the highest absolute value but put back the sign)
#    (keep max or min signed as original)
# 1. ABSOLUTE VALUE - hail.expr.functions.abs(x)
#    https://hail.is/docs/0.2/functions/numeric.html#hail.expr.functions.abs
# 2. MAX - hail.expr.aggregators.max(expr)
#    https://hail.is/docs/0.2/aggregators.html#hail.expr.aggregators.max
# 3. MIN - hail.expr.aggregators.min(expr)
#    https://hail.is/docs/0.2/aggregators.html#hail.expr.aggregators.min
##############

# 18. Annotate "BaseQRankSum":
# split_mt = split_mt.annotate_rows(info=split_mt.info.annotate(
#     BaseQRankSum = hl.float64(hl.agg.approx_median(split_mt.gvcf_info.BaseQRankSum))
# ))
# 18-1. Annotate "BaseQRankSum_highest_abs_max_min":
split_mt = split_mt.annotate_rows(info=split_mt.info.annotate(
    BaseQRankSum = hl.float64(
        hl.if_else(
            hl.abs(hl.agg.min(split_mt.gvcf_info.BaseQRankSum)) > hl.abs(hl.agg.max(split_mt.gvcf_info.BaseQRankSum)),
            hl.agg.min(split_mt.gvcf_info.BaseQRankSum),
            hl.agg.max(split_mt.gvcf_info.BaseQRankSum)
        ))))

# 19. Annotate "ReadPosRankSum":
# split_mt = split_mt.annotate_rows(info=split_mt.info.annotate(
#     ReadPosRankSum = hl.float64(hl.agg.approx_median(split_mt.gvcf_info.ReadPosRankSum))
# ))
# 19-1. Annotate "ReadPosRankSum_highest_abs_max_min":
split_mt = split_mt.annotate_rows(info=split_mt.info.annotate(
    ReadPosRankSum = hl.float64(
        hl.if_else(
            hl.abs(hl.agg.min(split_mt.gvcf_info.ReadPosRankSum)) > hl.abs(hl.agg.max(split_mt.gvcf_info.ReadPosRankSum)),
            hl.agg.min(split_mt.gvcf_info.ReadPosRankSum),
            hl.agg.max(split_mt.gvcf_info.ReadPosRankSum)
        ))))

# 20. Annotate "MQRankSum":
# split_mt = split_mt.annotate_rows(info=split_mt.info.annotate(
#     MQRankSum = hl.float64(hl.agg.approx_median(split_mt.gvcf_info.MQRankSum))
# ))
# 20-1. Annotate "MQRankSum_highest_abs_max_min":
split_mt = split_mt.annotate_rows(info=split_mt.info.annotate(
    MQRankSum = hl.float64(
        hl.if_else(
            hl.abs(hl.agg.min(split_mt.gvcf_info.MQRankSum)) > hl.abs(hl.agg.max(split_mt.gvcf_info.MQRankSum)),
            hl.agg.min(split_mt.gvcf_info.MQRankSum),
            hl.agg.max(split_mt.gvcf_info.MQRankSum)
        ))))

# 21. Annotate "ExcessHet": 
# Using hardy_weinberg_test()
# then, calculate phredPval = -10.0 * Math.log10(pval)
# IN JAVA:
# if (pval < 10e-60) {
#     return Pair.of(sampleCount, PHRED_SCALED_MIN_P_VALUE);
# }
# final double phredPval = -10.0 * Math.log10(pval);
# private static final double MIN_NEEDED_VALUE = 1.0E-16;
# private static final boolean ROUND_GENOTYPE_COUNTS = true;
# public static final double PHRED_SCALED_MIN_P_VALUE = -10.0 * Math.log10(MIN_NEEDED_VALUE);
MIN_NEEDED_VALUE = 1.0E-16
PHRED_SCALED_MIN_P_VALUE = -10.0 * math.log10(MIN_NEEDED_VALUE)

split_mt = split_mt.annotate_rows(
    hwe = hl.agg.hardy_weinberg_test(split_mt.GT)
)
split_mt = split_mt.annotate_rows(
    info=split_mt.info.annotate(
        ExcessHet = hl.if_else(
            split_mt.hwe.het_freq_hwe < 10e-60,
            PHRED_SCALED_MIN_P_VALUE,
            (-10 * hl.log10(split_mt.hwe.het_freq_hwe))
        )
    )
)

# 22. Annotate "InbreedingCoeff":
# Monkol's recommandation: 
# https://broadinstitute.github.io/gnomad_methods/_modules/gnomad/utils/annotations.html#bi_allelic_site_inbreeding_expr
# ##FILTER=<ID=InbreedingCoeff,Description="GATK InbreedingCoeff < -0.3">
split_mt = split_mt.annotate_rows(
    info=split_mt.info.annotate(
        InbreedingCoeff = bi_allelic_site_inbreeding_expr(split_mt.GT)
    )
)
# # 22-1. Phard the value:
# split_mt = split_mt.annotate_rows(
#     info=split_mt.info.annotate(
#         Phard_InbreedingCoeff = hl.if_else(
#             split_mt.info.InbreedingCoeff < 10e-60,
#             PHRED_SCALED_MIN_P_VALUE,
#             (-10 * hl.log10(split_mt.info.InbreedingCoeff))
#         )
#     )
# )

# 23. Densify MT:
split_mt_dense = hl.experimental.densify(split_mt)

# 24. Drop fields:
drop_fields = ["gvcf_info","QUALapprox","VarDP","hwe"]
split_mt_dense = split_mt_dense.drop(*drop_fields)
split_mt_dense = split_mt_dense.annotate_rows(
    info=split_mt_dense.info.drop(*["SB"])
)

# 25. OUTPUT MT: output_MT_path
split_mt_dense.write(OUTPUT_denseMT, overwrite=True)
split_mt_dense = hl.read_matrix_table(OUTPUT_denseMT)

# 26. Filter out ExcessHet > 54.69:
split_mt_dense = split_mt_dense.filter_rows(split_mt_dense.info.ExcessHet <= 54.69)

# 27. Create Mata-data:
meta = \
{
    'filter': {'LowQual': {'Description': 'Low quality'}},
    'info': {'AC': {'Description': 'Allele count in genotypes, for each ALT allele, in the same order as listed',
                    'Number': 'A',
                    'Type': 'Integer'},
             'AC_raw': {'Description': 'Allele count in genotypes, for each ALT allele, in the same order as listed (AC_raw as the sum of adj and non-adj groups)',
                    'Number': 'A',
                    'Type': 'Integer'},
             'SOR': {'Description': 'Symmetric Odds Ratio of 2x2 contingency table to detect strand bias',
                     'Number': '1',
                     'Type': 'Float'},
             'ReadPosRankSum': {'Description': 'Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias. (Get highest abs(max/min) of max/min in cohort.)',
                                'Number': '1',
                                'Type': 'Float'},
             # 'RAW_MQandDP': {'Description': 'Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.',
             #                 'Number': '2',
             #                 'Type': 'Integer'},
             'AN': {'Description': 'Total number of alleles in called genotypes. (sample number * 2)',
                    'Number': '1',
                    'Type': 'Integer'},
             'InbreedingCoeff': {'Description': 'Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation',
                                 'Number': '1',
                                 'Type': 'Float'},
             'AF': {'Description': 'Allele Frequency, for each ALT allele, in the same order as listed. (AC_raw/AN)',
                    'Number': 'A',
                    'Type': 'Float'},
             'FS': {'Description': "Phred-scaled p-value using Fisher's exact test to detect strand bias",
                    'Number': '1',
                    'Type': 'Float'},
             'DP': {'Description': 'Approximate read depth; some reads may have been filtered',
                    'Number': '1',
                    'Type': 'Integer'},
             'BaseQRankSum': {'Description': 'Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities. (Get highest abs(max/min) of max/min in cohort.)',
                              'Number': '1',
                              'Type': 'Float'},
             # 'MLEAF': {'Description': 'Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed',
             #           'Number': 'A',
             #           'Type': 'Float'},
             # 'MLEAC': {'Description': 'Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed',
             #           'Number': 'A',
             #           'Type': 'Integer'},
             'MQ': {'Description': 'RMS Mapping Quality', 'Number': '1', 'Type': 'Float'},
             'QD': {'Description': 'Variant Confidence/Quality by Depth',
                    'Number': '1',
                    'Type': 'Float'},
             # 'SB': {'Description': 'Strand Bias',
             #         'Number': '4',
             #         'Type': 'Integer'},
             'MQRankSum': {'Description': 'Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities. (Get highest abs(max/min) of max/min in cohort.)',
                           'Number': '1',
                           'Type': 'Float'},
             'ExcessHet': {'Description': 'Phred-scaled p-value for exact test of excess heterozygosity',
                           'Number': '1',
                           'Type': 'Float'}
            },
    'format': {'GQ': {'Description': 'Genotype Quality',
                      'Number': '1',
                      'Type': 'Integer'},
               'SB': {'Description': "Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.",
                      'Number': '4',
                      'Type': 'Integer'},
               'AD': {'Description': 'Allelic depths for the ref and alt alleles in the order listed',
                      'Number': 'R',
                      'Type': 'Integer'},
               'PID': {'Description': 'Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group',
                       'Number': '1',
                       'Type': 'String'},
               'GT': {'Description': 'Genotype', 'Number': '1', 'Type': 'String'},
               'MIN_DP': {'Description': 'Minimum DP observed within the GVCF block',
                          'Number': '1',
                          'Type': 'Integer'},
               'PGT': {'Description': 'Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another',
                       'Number': '1',
                       'Type': 'String'},
               'PL': {'Description': 'Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification',
                      'Number': 'G',
                      'Type': 'Integer'},
               'DP': {'Description': 'Approximate read depth (reads with MQ=255 or with bad mates are filtered)',
                      'Number': '1',
                      'Type': 'Integer'},
               'RGQ': {'Description': 'Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)',
                       'Number': '1',
                       'Type': 'Integer'},
               'PS': {'Description': 'Phasing set (typically the position of the first variant in the set)',
                      'Number': '1',
                      'Type': 'Integer'}}}

# 5. OUTPUT VCF: out name can .bgz and add ", tabix=True"
hl.export_vcf(split_mt_dense, OUTPUT_VCF, metadata=meta)
