
""" Snake pipeline for running cellranger with CITE-seq data """

# Configure shell for all rules 
shell.executable("/bin/bash")
shell.prefix("set -o nounset -o pipefail -o errexit -x; ")
import subprocess # check if necessary
import glob
import os 
import re
import pandas as pd


# Parameters from config.yaml
SAMPLE_TABLE    = config["SAMPLE_TABLE"]
GENOME          = config["GENOME"]
COMBINED_GENOME = config["COMBINED_GENOME"]
GTF             = config["GTF"]
COMBINED_GTF    = config["COMBINED_GTF"]
RESULTS         = config["RESULTS"]
ADAPTORS        = config["ADAPTORS"]
TRIM_METHOD     = config["TRIM_METHOD"]
PROJECT         = config["PROJECT"]

# Pull out sample names and fastq files
SAMPLE_LIST = pd.read_table(config["SAMPLE_TABLE"]).set_index("sample", drop = False)
SAMPLES = SAMPLE_LIST.index.values
IS_PAIRED = "fastq2" in SAMPLE_LIST.columns

# Pull out if there are any samples with spike-ins
if "spike_in" in SAMPLE_LIST.columns:
    SPIKE_SAMPLES = SAMPLE_LIST[SAMPLE_LIST.spike_in].index.values
    NON_SPIKE_SAMPLES = SAMPLE_LIST[~SAMPLE_LIST.spike_in].index.values
else:
    NON_SPIKE_SAMPLES = SAMPLE_LIST
    SPIKE_SAMPLES = []

if IS_PAIRED:
    CMD_PARAMS = config["PE"]
else:
    CMD_PARAMS = config["SE"]

# Function to check paths for input files/directories
def _check_path(path):
    if os.path.exists(path):
        return os.path.abspath(path)
    else:
        sys.exit("ERROR: " + path + " does not exist.")



# Set sample/group names
SAMPLES = [x.strip() for x in SAMPLES]


# Check/set directory/file paths
RESULTS  = _check_path(RESULTS)
GENOME   = _check_path(GENOME)

# Only check the combined genome if combined samples are present in the sample sheet
if len(SPIKE_SAMPLES) > 0:
    COMBINED_GENOME = _check_path(COMBINED_GENOME)

# Final output files
rule all:
    input:
    	# QC output
    	expand(
    		"{results}/fastqc/fastqc_{sample}.txt",
    		results = RESULTS, sample = SAMPLES),
        # Output of adaptor trimming
        expand(
            "{results}/logs/{trim_method}_trim_{sample}.txt",
            results = RESULTS, sample = SAMPLES,
            trim_method = TRIM_METHOD),
        # STAR output
        expand(
            "{results}/star_{trim_method}_trim/{sample}_Aligned.sortedByCoord.out.bam",
            results = RESULTS, sample = SAMPLES,
            trim_method = TRIM_METHOD
            ),
        # featureCounts output
        expand(
            "{results}/featureCount_{trim_method}_trim/{sample}_countsOutput",
            results = RESULTS, sample = SAMPLES,
            trim_method = TRIM_METHOD
            ),

        # Count table
        expand(
            "{results}/{project}_countTable_{trim_method}.txt",
            results = RESULTS, project = PROJECT,
            trim_method = TRIM_METHOD
            )

# Snakes to run
include: "src/rules/fastqc.snake"
include: "src/rules/trimming.snake"
include: "src/rules/star.snake"
include: "src/rules/featurecounts.snake"