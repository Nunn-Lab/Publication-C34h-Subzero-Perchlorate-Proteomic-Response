# C34h-Subzero-Perchlorate-Proteomic-Response
Supporting documentation for publication.
The documents included in this repository will allow you to recreate the analysis and figures in the paper.

**Figure 1A. NMDS of all time points between Control and WCL medium.**\
Folder: NMDS\
R script: NMDS.R\
Input files:\
    - ABACUS_output.csv\
    - metadata.csv

**Figure 1B. Weighted gene correlation network analysis (WGCNA) of different culture conditions (Time by day, WCL Treatment, Temperature, and estimated growth curve phase of lag, log, stationary or death)**\
Folder: WGCNA\
R script: WGCNA Prep.R must be run first to produce the necessary files to create the plot in WGCNA Plots.R\
Input files:\
    - ABACUS_output.csv\
    - metadata.csv\
Additional Output files: Protein IDs within each WGCNA module\
    - module_blue_prot.txt\
    - module_brown_prot.txt\
    - module_turquoise_prot.txt

**Figure 2. Abundance Analysis (QSPEC) of Control vs. WCL medium.**\
Folder: Volcano Plot\
R script: Volcano Plot.R\
Input files:\
    - 34h_annotations.csv\
    - ABACUS_output.csv\
    - QSPEC_CW.txt_qspec_fdr\
    - module_blue_prot.txt\
    - module_brown_prot.txt\
    - module_turquoise_prot.txt

**Figure 3. Enrichment results for blue (blue module) and turquoise (turquoise module) WGCNA module in terms of protein ratio**\
Folder: Enrichment\
R script: Module Enrichment.R\
Input file:\
    - Ratio Comp Summary.txt\
    - blue_uniprot_geneID.tsv\
    - turquoise_uniprot_geneID.tsv\
    - module_blue_prot.txt\
    - module_turquoise_prot.txt\
    - module_brown_prot.txt

**Supplemental Figure 1. Bacterial growth curves of C34h at -1C and -5C under different media conditions with estimated growth curves and loess smoothing**\
Folder: Cell Count\
R script: Cell Count & OD.R\
Input files:\
    - Direct_Count_and_OD.xlsx

**Other files included:**
biostats.R and documentation (source package used for some of the plots)
plainstat.R (code containing commonly used personal functions)
