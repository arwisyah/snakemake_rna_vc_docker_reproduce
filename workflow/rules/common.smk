#######################
import errno
import pandas as pd
import os
import multiprocessing
from snakemake.utils import validate


validate(config, schema="../schemas/config.schema.yaml")


samples = pd.read_table(config.get("samples"), index_col="sample")
units = pd.read_table(config.get("units"), index_col=["unit"], dtype=str)


## pipeline-related functions
def get_unit_fastqs(wildcards, samples, label='units',read_pair='fq'):
    for unit_set in samples.loc[wildcards.sample,[label]]:
        wildcards.sample
    return [expand_filepath(units.loc[x,[read_pair]].dropna()[0]) for x in unit_set.split(',')]


## filepath functions
def resolve_results_filepath(basepath, outname):
    return os.path.join(basepath, outname)


def expand_filepath(filepath):
    filepath = os.path.expandvars(os.path.expanduser(filepath))
    if not os.path.isabs(filepath):
        raise FileNotFoundError(
            errno.ENOENT,
            os.strerror(errno.ENOENT) + " (path must be absolute)",
            filepath,
        )
    return filepath


def resolve_single_filepath(basepath, filename):
    return os.path.join(basepath, filename)


def return_res_dir(output,elements=None):
    path_list=output[0].split(os.sep)
    path_list2=path_list[0:(len(path_list)-elements)]
    res_path=os.path.join("/",os.path.join(*(path_list2)))
    return res_path

## functions for system resources
def cpu_count():
    return multiprocessing.cpu_count()


def conservative_cpu_count(reserve_cores=1, max_cores=5):
    cores = max_cores if cpu_count() > max_cores else cpu_count()
    return max(cores - reserve_cores, 1)