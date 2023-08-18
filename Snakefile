# This is a makefile for batch analyzing solar orbiter waveforms ideantified by 
# CNN as dust impacts. Every file of input produces a stats file and some plots 
# as a byproduct. 

# Define input and output directories
datadir = "997_data/solo_amplitude_data"
statsdir = "997_data/solo_features"
plotsdir = "998_generated/solo_statistics"

# Define the input and output file patterns
data_pattern = datadir + "/{sample}.txt"
stats_pattern = statsdir + "/{sample}.npz"
plots_pattern = plotsdir + "/*.png"

# so all the input data are like this
datafiles = glob_wildcards(data_pattern).sample

# Define the final rule to aggregate all outputs
rule all:
    input:
        expand(stats_pattern, sample=datafiles)
    shell:
        """
        python 003_solar_orbiter/solo_delay_plots.py
        """

# Define the rule for generating statistics
rule generate_stats:
    input:
        data = data_pattern
    output:
        stats = stats_pattern
    conda: 'environment.yml'
    shell:
        """
        python 003_solar_orbiter/solo_feature_extraction.py {input.data} {output.stats}
        """
