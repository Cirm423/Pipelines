__author__ = "Maarten van der Sande"
__copyright__ = "Copyright 2020, Maarten van der Sande"
__email__ = "M.vanderSande@science.ru.nl"
__license__ = "MIT"


from snakemake.shell import shell
import os

# Optional parameters
provider = snakemake.params.get("provider", "UCSC")
assembly = snakemake.params.get("assembly")

annotation = ""
if any(["annotation" in out for out in snakemake.output]):
    annotation = "--annotation"

# parse the genome dir
genome_dir = "./"
if snakemake.output[0].count("/") > 1:
    genome_dir = "/".join(snakemake.output[0].split("/")[:-2])
    #genome_dir = os.path.dirname(os.path.dirname(snakemake.output[0]))

log = snakemake.log

# Finally execute genomepy
shell(
    """
    # install the genome
    genomepy install {assembly} \
    -p {provider} {annotation} -g {genome_dir} >> {log} 2>&1
    """
)