# TBD

# Define input and output directories
datadir = "C:/Users/skoci/Disk Google/000 Škola/UIT/getting data/solo/rpw/tds_wf_e"
statsdir = "998_generated/stats"

# Define the input and output file patterns
data_pattern = datadir + "/{sample}.cdf"
stats_pattern = statsdir + "/{sample}.txt"

# so all the input data are like this
datafiles = glob_wildcards(data_pattern).sample

# Define the final rule to aggregate all outputs
rule all:
    input:
        expand(stats_pattern, sample=datafiles)
    shell:
        """
        echo %time%
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
        python nano_produce_waveforms.py "{input.data}" "{output.stats}"
        """
