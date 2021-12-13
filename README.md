# NovFamilies
Scripts used for generating the results presented in Functional and evolutionary significance of unknown genes from uncultivated taxa.

## Database building

- create_families_json.py: generate json including all the gene family information (members, taxonomic distribution, ...)

## Figures

- ``generate_data_fig1d.py``: Script generating the proportion of uncultivated taxa and number of novel families per genome in each order of the GTDB tree. This data was later used to generate the figure in itol (https://itol.embl.de/) 
- ``figure2_A_matrix.py``: Script generating the KO presence / absence matrix in figure 2A.
- ``statistical_analysis_and_figures.r``: Script including the code for generating the remaining figures and some statistical tests. 

## Analysis

- 
- ``get_score_per_pos.py``: Calculate the genomic context conservation each novel family.
- ``genomic_context_conservation_score2json.py``: Get genomic context conservation score to json. 
- ``synapos.py``: Calculate synapomorphic families. 
