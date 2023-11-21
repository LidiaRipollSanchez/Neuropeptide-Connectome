# [The neuropeptidergic connectome of Caenorhabditis elegans](https://www.cell.com/action/showPdf?pii=S0896-6273%2823%2900756-0)
Lidia Ripoll-Sánchez, Jan Watteyne,  HaoSheng Sun, Robert Fernandez,  Seth R Taylor, Alexis Weinreb, Mark Hammarlund, David M Miller III,  Oliver Hobert,  Isabel Beets, Petra E Vértes,  William R Schafer

<img width="748" alt="Screenshot 2022-11-09 at 17 25 27" src="https://user-images.githubusercontent.com/86192587/200898826-29f869e6-137e-45d5-a99d-5e8a5f7edeec.png">

**Abstract**

Efforts are currently ongoing to map synaptic wiring diagrams or connectomes in order to understand the neural basis of brain function. However, chemical synapses represent only one type of functionally important neuronal connection; in particular, extrasynaptic, “wireless” signaling by neuropeptides is widespread and plays essential roles in all nervous systems. By integrating single-cell anatomical and gene expression datasets with a biochemical analysis of receptor-ligand interactions, we have generated a draft connectome of neuropeptide signaling in the C. elegans nervous system. This connectome is characterized by a high connection density, extended signaling cascades, autocrine foci, and a decentralized topology, with a large, highly interconnected core containing three constituent communities sharing similar patterns of input connectivity. Intriguingly, several of the most important nodes in this connectome are little-studied neurons that are specialized for peptidergic neuromodulation. We anticipate that the C. elegans neuropeptidergic connectome will serve as a prototype to understand basic organizational principles of neuroendocrine signaling networks.

**About this repository**

This repository contains scripts, input data, processed data files, figures and supplementary data for all analysis presented in [Ripoll-Sánchez et al., 2023](https://www.cell.com/action/showPdf?pii=S0896-6273%2823%2900756-0)

Original scRNAseq data used for analysis comes from [the CeNGEN consortia](https://www.sciencedirect.com/science/article/pii/S0092867421007583?via%3Dihub#fig3) and neuropeptide-receptor deorphanisation data comes from [Beets et al., 2023](https://www.cell.com/cell-reports/pdf/S2211-1247(23)01069-0.pdf). Synaptic connectivity and anatomical structure of the nervous system of C.elegans come from a list of resources including [White et al., 1986](https://doi.org/10.1098/rstb.1986.0056), [Witvliet et al., 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8756380/pdf/nihms-1766209.pdf), [Varshney et al., 2011](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1001066)

All main figures analysis is done at CeNGEN threshold 4 and EC50 deorphanisation data thresholded to 500nM, to explore the networks with alternative thresholds check the Sensitivity Analysis data. 
All scripts for running the analysis shown in the paper, and all the input files required to run it are in the Scripts & data FOLDER.

**Figure 1**

* 17072023_Main_figures (1-4).pdf (Figure 1): The figure
* 14062022_visualisation_short_range_cytoscape_colorpalete.pdf.cys (FOLDER Adjacency matrices for networks): Cytoscape file for Figure 1
* 14062022_visualisation_short_range_cytoscape_colorpalete.pdf (FOLDER Adjacency matrices for networks): Network representation from Figure 1
* 08062023_neuropeptide_connectome_short_range_model.csv (FOLDER Adjacency matrices for networks): Dataset used to build Cytoscape figure, network expressed as edge list

**Figure 2**

*17072023_Main_figures (1-4).pdf (Figure 2): The figure

**Figure 3**

* 17072023_Main_figures (1-4).pdf (Figure 3)(FOLDER Figures): The figure
* Worm_head_midbody_tail3_horizontal.pdf (FOLDER Figures): worm outline figure
* GPCR_with_one_ligand (FOLDER Models Comparison): identities of GPCRs that only have 1 ligand
* cognate_Neuropeptide_to_GPCR_with_one_ligand (FOLDER Models Comparison): the identities of the cognate neuropeptides
* Neuron_expression_GPCR_1_ligand (FOLDER Models Comparison): Expression of 1 ligand GPCRs by neuron
* Neuron_expression_cognate_NPP_to_GPCR_1_ligand (FOLDER Models Comparison): Expression of cognate NPPs that specifically bind only GPCRs that bind only 1 ligand

**Figure 4**

* 17072023_Main_figures (1-4).pdf (Figure 4)(FOLDER Figures): The figure
* Barplot_Neuron_with_GPCR.png (FOLDER Expression analysis): Number of neurons that express each GPCR, used for the assortativity analysis
* Barplot_Neuron_with_neuropeptide.png (FOLDER Expression analysis): Number of neurons that express each NPP, used for the assortativity analysis
* FOLDER Individual_net_short_r (FOLDER Adjacency matrices for networks): Adjacency matrices of all 92 NPP-GPCR networks short-range used to build the graphs in the main figure and figures S3 and S4
* FOLDER Individual_net_mid_r (FOLDER Adjacency matrices for networks): Adjacency matrices of all 92 NPP-GPCR networks mid range used to build supplementary figures S5 and S6
* neuropeptide_pairs.csv (FOLDER Adjacency matrices for networks): Order of the 92 NPP-GPCR networks in list format used for supplementary figures S3-S6 and for the order of networks given in the individual_net folders

**Figure 5**

* 17072023_Main_figures (4-8).pdf (Figure 5) (FOLDER Figures): The figure
* 08062023_neuropeptide_connectome_short_range_model.csv (FOLDER Adjacency matrices for networks): Dataset used to build matrix depiction in colour in figure 
* 08062023_neuropeptide_connectome_mid_range_model.csv (FOLDER Adjacency matrices for networks): Dataset used to build matrix depiction in grey in figure 
* Barplot_NPP_per_neuron.png (FOLDER Expression analysis): number of NPPs expressed per neuron, neurons sorted in alphabetical order
* Barplot_GPCR_per_neuron.png (FOLDER Expression analysis): number of GPCRs expressed per neuron, neurons sorted in alphabetical order
* 26012022_num_neuronID.txt (FOLDER Scripts & data): neurons that correspond to each row and column in the adjacency matrices
* 072022_anatomical_class.csv (FOLDER Scripts & data): classification of each neuron in a neuronal type and matrix of anatomical proximity of neuronal processes that was used to define the short-range constraints for the short range network. 

**Figure 6**

* 17072023_Main_figures (4-8).pdf (Figure 6) (FOLDER Figures): The figure
* Degree_table_networks_comparison.xlsx (FOLDER Network Topology Measures /SUBFOLDER Degree) : Excel file with data comparing the number of connections in short, long and mid range neuropeptide connectomes (also groups)
* All grouped neuropeptide networks short range rich club.csv (FOLDER Network Topology Measures /SUBFOLDER Rich Club): short range rich club data shown in figure

**Figure 7**

* 17072023_Main_figures (4-8).pdf (Figure 7) (FOLDER Figures): The figure
* 08062023_Groups_mid_range_network_table.csv (FOLDER Mesoscale structure): table with the t-SNE values, the neuron names, their group allocation and their indegree used in the figure

**Figure 8**

* 17072023_Main_figures (4-8).pdf (Figure 8) (FOLDER Figures): The figure

**Supplementary data**

* 07072022_NPP_GPCR_networks_long_range_model (FOLDER Adjacency matrices for networks): 91 individual NPP-GPCR networks long-range model
* 07072022_neuropeptide_connectome_long_range_model (FOLDER Adjacency matrices for networks): all grouped NPP-GPCR networks long-range model

**Extended sensitivity data (further analysis than what show in the paper, All files in the Sensitivity Analysis FOLDER)**

* CSV FOLDERS include the adjacency matrices found in the Adjacency matrices for networks but at the threshold stated (T4 = CeNGEN threshold 4, T3 = CeNGEN threshold 3, the concentration at the end of the name indicates the EC50 threshold for deorphanisation, the neuropeptide_pairs.csv file on each folder indicates the identities of the individual networks )
* Degree FOLDER includes the degree distributions and comparisons of all the networks in the CSV FOLDERS and their comparison to the standard threshold 4 EC50 500nM
* Rich Club folder includes the rich club of all the threshold 4 networks in the CSV FOLDERS and their comparison to the standard threshold 4 EC50 500nM
* The other files include the images of all the networks in CSV FOLDERS, the models comparison sensitivity analysis with T4 and T3, the 
