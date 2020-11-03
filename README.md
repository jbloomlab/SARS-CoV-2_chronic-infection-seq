## SARS-CoV-2 Chronic Infection Sequencing
*2020-11-03*

This repo contains code for analyzing deep-sequencing data from [Dr. Jonathan Li's group](https://jonathanlilab.bwh.harvard.edu/) describing an immunocompromised patient with SARS-CoV-2. There is longitudinal sequencing from this patient that demonstrates significant evolution in the RBD. 

There is Illumina deep-sequencing at nine timepoints performed on days 18, 25, 75, 81, 128, 130, 143, 146, and 152. Interestingly, there is one non-synonymous mutation that arises, Q493K, that is a reported Regeneron REGN10933 escape mutation.

We have a collaboration set with the Li group and they'll be sending over the FASTQ files hopefully this week. The general story seems to revolve around evidence that there is selection in the RBD against the immune system or the REGN antibody cocktail. 

### Notes on the Paper

An interesting facet of the results is the extent of amino acid changes observed in the Spike gene. Amin said that this is very unlike what others have seen. Amino acid changes were predominantly in the spike gene and receptor binding domain, which comprises 13% and 2% of the viral genome, but harbored 57% and 38% of the 93 observed changes, respectively.

The viral consensus sequences for these isolates can be located with GISAID accession numbers EPI_ISL_593480, EPI_ISL_593557, EPI_ISL_593558, EPI_ISL_593555, EPI_ISL_593478, EPI_ISL_593479, EPI_ISL_593556, EPI_ISL_593553, EPI_ISL_593554

Assmebly and variant calling was done with a SARS-CoV-2 specific module of [IRMA](https://wonder.cdc.gov/amd/flu/irma/). The variants were called down to 5%. We can probably refine this more. 


