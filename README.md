# [The neuropeptidergic connectome of Caenorhabditis elegans]([https://www.biorxiv.org/content/10.1101/2022.10.30.514396v2.full](https://www.cell.com/action/showPdf?pii=S0896-6273%2823%2900756-0))
Lidia Ripoll-Sánchez, Jan Watteyne,  HaoSheng Sun, Robert Fernandez,  Seth R Taylor, Alexis Weinreb, Mark Hammarlund, David M Miller III,  Oliver Hobert,  Isabel Beets, Petra E Vértes,  William R Schafer

<img width="748" alt="Screenshot 2022-11-09 at 17 25 27" src="https://user-images.githubusercontent.com/86192587/200898826-29f869e6-137e-45d5-a99d-5e8a5f7edeec.png">

**Abstract**

Efforts are currently ongoing to map synaptic wiring diagrams or connectomes in order to understand the neural basis of brain function. However, chemical synapses represent only one type of functionally important neuronal connection; in particular, extrasynaptic, “wireless” signaling by neuropeptides is widespread and plays essential roles in all nervous systems. By integrating single-cell anatomical and gene expression datasets with a biochemical analysis of receptor-ligand interactions, we have generated a draft connectome of neuropeptide signaling in the C. elegans nervous system. This connectome is characterized by a high connection density, extended signaling cascades, autocrine foci, and a decentralized topology, with a large, highly interconnected core containing three constituent communities sharing similar patterns of input connectivity. Intriguingly, several of the most important nodes in this connectome are little-studied neurons that are specialized for peptidergic neuromodulation. We anticipate that the C. elegans neuropeptidergic connectome will serve as a prototype to understand basic organizational principles of neuroendocrine signaling networks.

**About this repository**

This repository contains codes and processed data files for analyses presented in this [bioRxiv paper](https://www.biorxiv.org/content/10.1101/2022.10.30.514396v2.full)

Original scRNAseq data used in this paper comes from [the CeNGEN consortia](https://www.sciencedirect.com/science/article/pii/S0092867421007583?via%3Dihub#fig3) and original deorphanisation data used in this paper comes from [this paper](https://www.biorxiv.org/content/10.1101/2022.10.30.514428v1.full)

**Figure 1**

* 14062022_visualisation_short_range_cytoscape_colorpalete.pdf.cys (FOLDER Adjacency matrices for networks): Cytoscape file for Figure 1
* 14062022_visualisation_short_range_cytoscape_colorpalete.pdf (FOLDER Adjacency matrices for networks): Network representation from Figure 1
* 07072022_NPP_GPCR_networks_short_range_model.csv (FOLDER Adjacency matrices for networks): Dataset used to build Cytoscape figure, network expressed as edge list


**Figure 3**

* Worm_head_midbody_tail3_horizontal.pdf (FOLDER Figures): worm outline figure
* GPCR_with_one_ligand (FOLDER Models Comparison): identities of GPCRs that only have 1 ligand
* cognate_Neuropeptide_to_GPCR_with_one_ligand (FOLDER Models Comparison): the identities of the cognate neuropeptides
* Neuron_expression_GPCR_1_ligand (FOLDER Models Comparison): Expression of 1 ligand GPCRs by neuron
* Neuron_expression_cognate_NPP_to_GPCR_1_ligand: Expression of cognate NPPs that specifically bind only GPCRs that bind only 1 ligand

**Figure 4**

* Barplot_Neuron_with_GPCR.png : Number of neurons that express each GPCR, used for the assortativity analysis
* Barplot_Neuron_with_neuropeptide.png : Number of neurons that express each NPP, used for the assortativity analysis
* 07072022_NPP_GPCR_networks_short_range_model (FOLDER Adjacency matrices for networks): Adjacency matrices of all 91 NPP-GPCR networks short-range used to build the graphs in the figure and figure S2
* 07072022_NPP_GPCR_networks_mid_range_model (FOLDER Adjacency matrices for networks): Adjacency matrices of all 91 NPP-GPCR networks mid range used to build supplementary figure S3
* 20220103_Order of 91 NPP/GPCR networks.csv : Order of the 91 NPP-GPCR networks in list format used for supplementary figure S2 and S3

**Figure 5**

* 07072022_neuropeptide_connectome_short_range_model.csv (FOLDER Adjacency matrices for networks): Dataset used to build matrix depiction in colour in figure 
* 07072022_neuropeptide_connectome_mid_range_model.csv (FOLDER Adjacency matrices for networks): Dataset used to build matrix depiction in grey in figure 
* Barplot_NPP_per_neuron.png : number of NPPs expressed per neuron, neurons sorted in alphabetical order
* Barplot_GPCR_per_neuron.png : number of GPCRs expressed per neuron, neurons sorted in alphabetical order
* 26012022_num_neuronID.txt : neurons that correspond to each row and column in the adjancency matrices
* 072022_anatomical_class.csv (FOLDER Main): classification of each neuron in a neuronal type and matrix of anatomical proximity of neuronal processes that was used to define the short-range constraints for the short range network. 

**Figure 6**

* 07072022_degrees (comparison between different connection range).xlsx : Excel file with data and comparisons between the short, long and mid range neuropeptide connectomes and the synaptic, gap junction and monoamine connectomes. 
* All grouped neuropeptide networks short range rich club.csv : short range rich club data shown in figure
* All grouped neuropeptide networks mid range rich club.csv : mid range rich club data

**Figure 7**

* 07072022_neuropeptide_connectome_mid_range_sorted_by_groups.csv (FOLDER Adjacency matrices for networks): Dataset used to build matrix depiction in the figure, resorted mid range matrix
* 07072022_Groups_mid_range_network_table.csv (FOLDER Adjacency matrices for networks): table with the t-SNE values, the neuron names, their group allocation and their indegree used in the figure

**Figure 8**

*

**Supplementary data**

* 07072022_NPP_GPCR_networks_long_range_model (FOLDER Adjacency matrices for networks): 91 individual NPP-GPCR networks long-range model
* 07072022_neuropeptide_connectome_long_range_model (FOLDER Adjacency matrices for networks): all grouped NPP-GPCR networks long-range model

