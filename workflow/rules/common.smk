import pandas as pd
import re
from snakemake.utils import validate


def get_assembly_name(assembly):
    rgx = re.compile(r".*\/(?P<FN>.*)\.((fasta|fa){1})((\.gz)?)")
    return rgx.search(assembly).group('FN')


def get_reads_names(reads):
    rgx = re.compile(r".*\/(?P<FN>.*)\.((fastq|fq){1})((\.gz)?)")
    return pd.Series([rgx.search(read).group('FN') for read in reads])

# path to assemly stored in config
FASTA = config["assembly"]

# assembly name used in output filenames
ASSEMBLY = get_assembly_name(FASTA)

# dataframe of info on samples #TODO deal not with only PE, but also SE
SAMPLES = pd.read_table(config["samples"], dtype=str)
# extract sample name from path of forward read if not filled
SAMPLES['sample_id'].fillna(get_reads_names(SAMPLES['forward_read']), inplace=True)
# make column "sample" as identificator of dataframe
SAMPLES.set_index("sample_id", drop=False)
validate(SAMPLES, schema="../schemas/samples.schema.yaml")

# window size for coverage
SIZE=config["coverage_stats_window_size_list"]


# def get_fasta(wildcards):
#     return SAMPLES[SAMPLES["assembly_name"]==wildcards.sample].assembly


def get_forward_fastq(wildcards):
    return SAMPLES[SAMPLES["sample_id"]==wildcards.sample].forward_read


def get_reverse_fastq(wildcards):
    return SAMPLES[SAMPLES["sample_id"]==wildcards.sample].reverse_read


def get_platform(wildcards):
    return SAMPLES[SAMPLES["sample"]==wildcards.sample].platform


def aggregate_file_names(wildcards, pattern):
     checkpoint_output = checkpoints.bcftools_vcf_subset.get().output[0]
     return expand(str(pattern),
                    subset=glob_wildcards(os.path.join(checkpoint_output, "{subset}/"+ASSEMBLY+".vcf.gz")).subset)

