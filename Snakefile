from paths import cdf_mamp_location
from paths import cdf_tswf_e_location


# Define input and output directories for padding
padding_datadir = cdf_tswf_e_location
padding_statsdir = "998_generated/stats/"


# Define input and output directories for mamp
mamp_datadir = cdf_mamp_location
mamp_statsdir = "998_generated/mamp_processed/"


#find all the inputs
padding_data, = glob_wildcards(padding_datadir+"{sample}.cdf")
mamp_data, = glob_wildcards(mamp_datadir+"{sample}.cdf")


# Define the final rule that requests all the outputs
rule all:
    input:
        expand(padding_statsdir+"{output}.txt", output=padding_data),
        expand(mamp_statsdir+"{output}.pkl", output=mamp_data)
    conda:
	"environment.yml"
    shell:
	"""
	python aggregate.py
	"""


#rule to do padding
rule padding_generate_stats:
    input:
        file = padding_datadir+"{sample}.cdf"
    output:
        file = padding_statsdir+"{sample}.txt"
    conda:
        "environment.yml"
    shell:
        """
        python nano_produce_waveforms.py "{input.file}" "{output.file}"
        """


#rule to analyze mamp
rule mamp_generate_stats:
    input:
        file = mamp_datadir+"{sample}.cdf"
    output:
        file = mamp_statsdir+"{sample}.pkl"
    conda:
        "environment.yml"
    shell:
        """
        python nano_mamp.py "{input.file}" "{output.file}"
        """
