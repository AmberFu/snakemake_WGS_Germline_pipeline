'''
https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#python

def do_something(data_path, out_path, threads, myparam):
    # python code

do_something(snakemake.input[0], snakemake.output[0], snakemake.threads, snakemake.config["myparam"])
'''

import os
import math
import subprocess
import time
from datetime import datetime
import pyspark
import hail as hl

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
logfile = "/storage1/fs1/jin810/Active/fup/hailLogs/hail_jointcalling_{}_{}.log".format(project_name, today)
### Hail Start...
hl.init(default_reference='GRCh38',sc=sc,log=logfile)

####################################

INPUT_GVCFs = snakemake.input.bgz_gvcfs
OUTPUT_VDS = snakemake.output.vds

combiner = hl.vds.new_combiner(
    output_path=OUTPUT_VDS,
    temp_path='/storage1/fs1/jin810/Active/fup/hail_tmp',
    gvcf_paths=INPUT_GVCFs,
    # vds_paths=vdses,
    use_genome_default_intervals=True,
)

combiner.run()

print(">>> Finished Hail Joint-calling!\n")