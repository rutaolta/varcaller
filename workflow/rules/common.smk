import pandas as pd
import re
from snakemake.utils import validate


def get_assembly_name(assemblies):
    rgx = re.compile(r".*\/(?P<FN>.*)\.((fasta|fa){1})((\.gz)?)")
    return [rgx.search(assembly).group('FN') for assembly in assemblies]


def get_reads_names(reads):
    rgx = re.compile(r".*\/(?P<FN>.*)\.((fastq|fq){1})((\.gz)?)")
    return [rgx.search(read).group('FN') for read in reads]


SAMPLES = pd.read_table(config["samples"], dtype=str)
SAMPLES.insert(2, "assembly_name", get_assembly_name(SAMPLES.assembly), True)
SAMPLES.set_index("assembly_name", drop=False)
SAMPLES.insert(4, "forward_read_name", get_reads_names(SAMPLES.forward_read), True)
SAMPLES.insert(5, "reverse_read_name", get_reads_names(SAMPLES.reverse_read), True)
validate(SAMPLES, schema="../schemas/samples.schema.yaml")
#TODO make index for sample unique identification, in case of several different reads for one assembly
#TODO deal not with only PE, but also SE

def get_fasta(wildcards):
    return SAMPLES[SAMPLES["assembly_name"]==wildcards.sample].assembly


def get_forward_fastq(wildcards):
    return SAMPLES[SAMPLES["assembly_name"]==wildcards.sample].forward_read


def get_reverse_fastq(wildcards):
    return SAMPLES[SAMPLES["assembly_name"]==wildcards.sample].reverse_read


def get_platform(wildcards):
    return #TODO for alignment bwa_map
