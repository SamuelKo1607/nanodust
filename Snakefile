# TBD

# Define input and output directories
datadir = "C:/Users/skoci/Disk Google/000 Škola/UIT/getting data/solo/rpw/tds_wf_e"
statsdir = "998_generated/stats"
mamp_cdf_dir = "C:/Users/skoci/Disk Google/000 Škola/UIT/getting data/solo/rpw/mamp"
mamp_pkl_dir = "998_generated/mamp_processed"

# Define the input and output file patterns
data_pattern = datadir + "/{sample}.cdf"
stats_pattern = statsdir + "/{sample}.txt"
data_mamp_pattern = mamp_cdf_dir + "/{sample}.cdf"
stat_mamp_pattern = mamp_pkl_dir + "/{sample}.pkl"

# so all the input data are like this
datafiles = glob_wildcards(data_pattern).sample
mamp_datafiles = glob_wildcards(data_mamp_pattern).sample

# Define the final rule to aggregate all outputs
rule all1:
    input:
        expand(stats_pattern, sample=datafiles)
    shell:
        """
        echo %time%
        """

rule all2:
    input:
        expand(stat_mamp_pattern, sample=mamp_datafiles)
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

rule mamp_generate_stats:
    input:
        data = data_mamp_pattern
    output:
        stats = stat_mamp_pattern
    conda: 'environment.yml'
    shell:
        """
        python -c 'from nano_mamp import main' 
	python main "{input.data}" "{output.stats}"
        """
