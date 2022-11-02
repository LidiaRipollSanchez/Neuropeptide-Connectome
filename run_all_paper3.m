%Run all at once

%Created by LIDIA RIPOLL SANCHEZ, for CENGEN scRNA-seq
%neuropeptide-receptor network project 12/07/2022
%This script builds the neuropeptide receptor networks, then plots them,
%masks them to make them unweighted and finally groups several networks
%depending on certain conditions given. 
version = 'R2021b'; %MATLAB version when this workspace was saved

%% ****************************************************************************************************************
%TO SAVE DATA 
path = '/scRNAseq_expression_threshold4_EC50_500nM_91/Networks/'; %create folder to save images
path2 = '/scRNAseq_expression_threshold4_EC50_500nM_91/Workspaces/'; %create folder to save images
addpath([ 'Results/scRNAseq_C.elegans' path]) %add the folder to the main path so it is easier to search 
addpath([ 'Results/scRNAseq_C.elegans' path2]) %add the folder to the main path so it is easier to search 
date = '20102022_'; %to save workspace with the date

%IMPORT DATA
expression_byneuron = readtable ('30072020_CENGEN_threshold4_expression_NPP_NPR_MR_LGC_allneurons'); %expression of NPP, GPCR and Mon receptors
connections_NPP = readtable ('20220103_Ligand-receptor interactions large-scale screen_per peptide and receptor gene_500nM_FINAL.csv'); %connections between NPP and GPCR
connections_Mon = readtable ('10052021_Ligand-receptor-interactions_Mon_from_Barry_Bentley.csv'); %connections between Mon and receptors
neuron_ID = readtable('26012022_num_neuronID.txt');
anatomical_class = readtable('072022_anatomical_class.csv'); %import list of neuron cell body and process location (in binary)

%RUN ANALYSIS
% Network_building(connections_NPP, expression_byneuron, path, path2, date, version)
% Network_measures(synA, gapA, Mon_all_grouped_networks, NPP_all_grouped_networks, NPP_all_grouped_networks_mid_range_connections, NPP_all_grouped_networks_short_range_connections, NPP_individual_networks_short_range_connections, path, path2, date, version)
% Network_computeRichClub(NPP_all_grouped_networks_short_range_SL_u, randNPP_NPGPCR_short_rangeu, 'All grouped neuropeptide networks short range rich club', (['Results/scRNAseq_C.elegans' path 'Topological_data_unweighted/Rich_club']), neuron_ID.nodeLabel.', @rich_club_bd, 'deg')
% Network_dimensionality_reduction(NPP_all_grouped_networks, NPP_all_grouped_networks_mid_range_connections, NPP_all_grouped_networks_short_range_connections, indeg_NPP_all_mid_range, ...
%     NPP_individual_networks_mid_range_connections_SL, neuron_ID.nodeLabel, anatomical_class.FinalClassification, pair_names_NPP, path, path2, date, version)
Network_models_comparison(connections_NPP, pair_names_NPP, NPP_all_grouped_networks, NPP_all_grouped_networks_mid_range_connections, ...
    NPP_all_grouped_networks_short_range_connections, NPP_individual_networks_short_range_connections, NPP_individual_networks_nerve_ring_connections, ...
    NPP_individual_networks_mid_range_connections, neuropeptide_networks_all, neuron_ID, class_Head, class_Midbody, class_Tail, path, path2, date, version)
% Individual_networks_analysis(NPP_individual_networks_short_range_connections, NPP_individual_networks_mid_range_connections, 50, path, path2, date)

%FUNCTIONS

function [expression_byneuron_ordered_unweighted_all, neuropeptide_networks_all, monoamine_networks_all, ...
    pair_names_NPP, pair_names_Mon, neuropeptides, receptors_NPP, NPP_all_grouped_names, NPP_all_grouped_networks...
    Mon_all_grouped_names, Mon_all_grouped_networks] = Network_building(connections_NPP, expression_byneuron, path, path2, date, version)

geneID = readtable('20210326_Gene_ID_GPCR_common_name.xlsx'); %NPP and GPCR gene common names in same order as the genes in the expression dataset
expression_byneuron_neurotransmitter = readtable('20210510_Neurotransmitters_Expression_(L.Pereira et al. 2015 Hobert lab).csv', 'ReadRownames', true); %import list of neurotrasmitter per neuron with rownames to sort neurons, includes expression of Mon

%% ****************************************************************************************************************
%ADAPT DATA
%****************************************************************************************************************

%ADAPT EXPRESSION DATA
[expression_byneuron_ordered] = Adapt_expression(expression_byneuron, neuron_ID.('nodeLabel'), geneID); %Adapt the expression values split by neuron type to neuron (from 180 columns to 302 columns)
[expression_byneuron_ordered_unweighted] = Mask_network(expression_byneuron_ordered); %Unweight expression values, when expression value is 1 otherwise is 0

%ADAPT NEUROTRANSMITTER EXPRESSION DATA
expression_byneuron_neurotransmitter = rows2vars(expression_byneuron_neurotransmitter(neuron_ID.nodeLabel, :)); %sort data based on neuron_ID
expression_byneuron_neurotransmitter.Properties.VariableNames([1]) = {'Gene'}; %change the name of the first column to Gene, since after transposing it now has the default name

%CREATE OVERALL EXPRESSION TABLE 
expression_byneuron_ordered_unweighted_all = [expression_byneuron_ordered_unweighted;expression_byneuron_neurotransmitter]; %add the CeNGEN expression unweighted dataset to the monoamine expression dataset

%% ****************************************************************************************************************
%BUILD INDIVIDUAL NETWORKS
%****************************************************************************************************************

[neuropeptide_networks_all,structure_NPP,pair_names_NPP] = Neuropeptide_Networks(connections_NPP, expression_byneuron_ordered_unweighted_all);
[monoamine_networks_all,structure_Mon,pair_names_Mon] = Neuropeptide_Networks(connections_Mon, expression_byneuron_ordered_unweighted_all);


%PLOT NETWORKS 
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Images_neuropeptide_networks_unweighted/Individual_NPP_networks'])  %create folder to save images
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Images_monoamine_networks/Individual_Monoamine_networks']) %create folder to save images
Plot_Networks(neuropeptide_networks_all,pair_names_NPP, ([path 'Images_neuropeptide_networks_unweighted/Individual_NPP_networks'])); %plot and save all adjacency matrices from neuropeptide_networks
Plot_Networks(monoamine_networks_all,pair_names_Mon, ([path 'Images_monoamine_networks/Individual_Monoamine_networks'])); %plot and save all adjacency matrices from monoamine_networks

%% ****************************************************************************************************************
%GROUPED NETWORKS 
%****************************************************************************************************************

%GROUP NETWORKS NEUROPEPTIDE
neuropeptides = matlab.lang.makeValidName(unique(connections_NPP.GeneLigand)); %create list of the neuropeptides uses in the edges to group based on them 
receptors_NPP = matlab.lang.makeValidName(unique(connections_NPP.GeneGPCR, 'stable')); %create list of the receptors used in the edges to group based on them, use stable to prevent unique from reordering alphabetically so that we group them in the order of connections
[NPP_all_grouped_names, NPP_all_grouped_networks] = Group_networks (pair_names_NPP, neuropeptide_networks_all, '\w*_\d*_\w*'); %group all networks
[NPP_grouped_names, NPP_grouped_networks] = Group_networks (pair_names_NPP, neuropeptide_networks_all, neuropeptides); %group based on neuropeptides 
[NPGPCR_grouped_names, NPGPCR_grouped_networks] = Group_networks (pair_names_NPP, neuropeptide_networks_all, receptors_NPP); %group based on neuropeptide receptors 
[NLP_NPP_grouped_names, NLP_NPP_grouped_networks] = Group_networks (pair_names_NPP, neuropeptide_networks_all, '(NLP|SNET|PDF|NTC)_\d*_'); %group based on NLP neuropeptides 
[FLP_NPP_grouped_names, FLP_NPP_grouped_networks] = Group_networks (pair_names_NPP, neuropeptide_networks_all, 'FLP_\d*_'); %group based on FLP neuropeptides 

%PLOT GROUPED NETWORKS
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Images_neuropeptide_networks_unweighted/NPP_all_grouped_networks'])  %create folder to save images
Plot_Networks(NPP_all_grouped_networks, 'NPP_all_grouped_networks', ([path 'Images_neuropeptide_networks_unweighted/NPP_all_grouped_networks']), neuron_ID.nodeLabel); %plot and save all networks grouped
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Images_neuropeptide_networks_unweighted/NPP_grouped_networks'])  %create folder to save images
Plot_Networks(NPP_grouped_networks, NPP_grouped_names, ([path 'Images_neuropeptide_networks_unweighted/NPP_grouped_networks']), neuron_ID.nodeLabel); %plot and save NPP grouped networks
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Images_neuropeptide_networks_unweighted/NPGPCR_grouped_networks'])  %create folder to save images
Plot_Networks(NPGPCR_grouped_networks, NPGPCR_grouped_names, ([path 'Images_neuropeptide_networks_unweighted/NPGPCR_grouped_networks']), neuron_ID.nodeLabel);%plot and save NPGPCR grouped networks
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Images_neuropeptide_networks_unweighted/NLP_NPP_grouped_networks'])  %create folder to save images
Plot_Networks(NLP_NPP_grouped_networks, 'NLP_NPP_grouped_networks', ([path 'Images_neuropeptide_networks_unweighted/NLP_NPP_grouped_networks']), neuron_ID.nodeLabel);%plot and save NPGPCR grouped networks
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Images_neuropeptide_networks_unweighted/FLP_NPP_grouped_networks'])  %create folder to save images
Plot_Networks(FLP_NPP_grouped_networks, 'FLP_NPP_grouped_networks', ([path 'Images_neuropeptide_networks_unweighted/FLP_NPP_grouped_networks']), neuron_ID.nodeLabel);%plot and save NPGPCR grouped networks

%GROUP NETWORKS MONOAMINE
monoamines = matlab.lang.makeValidName(unique(connections_Mon.GeneLigand)); %create list of the neuropeptides uses in the edges to group based on them 
[Mon_all_grouped_names, Mon_all_grouped_networks] = Group_networks (pair_names_Mon, monoamine_networks_all, '\w*_\w*'); %group all networks
[Mon_grouped_names, Mon_grouped_networks] = Group_networks (pair_names_Mon, monoamine_networks_all, monoamines); %group based on individual monoamines  

%PLOT GROUPED NETWORKS
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Images_monoamine_networks/All_monoamine_grouped']) %create folder to save images
Plot_Networks(Mon_all_grouped_networks, 'Monoamine_all_grouped', ([path 'Images_monoamine_networks/All_monoamine_grouped']), neuron_ID.nodeLabel); %plot and save all networks grouped
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Images_monoamine_networks/Monoamine_grouped']) %create folder to save images
Plot_Networks(Mon_grouped_networks, Mon_grouped_names, ([path 'Images_monoamine_networks/Monoamine_grouped']), neuron_ID.nodeLabel); %plot and save all networks grouped

% ****************************************************************************************************************
%SLICE  NETWORKS BASED ON PROCESS LOCATION (SHORT, MID and LONG RANGE
%CONNECTIONS)
%****************************************************************************************************************

%DEFINE MID RANGE SECTIONS
columns_Head = [6:22, 34]; %column numbers that correspond to processes that go through the head (everything up to NR + VNC anterior, includes pharynx)
columns_Midbody = [32:33, 35:40]; %column numbers that correspond to processes that go through the midbody (all sublaterals, VNC, DNC and canal)
columns_Tail = [44:46]; %column numbers that correspond to processes that go through the tail

%BUILD NETWORKS BY NEURON PROCESS POSITION
[NPP_all_grouped_networks_individual_short_range, NPP_all_grouped_networks_individual_mid_range, class_list_All_grouped, class_names_All_grouped, class_Head, class_Midbody, class_Tail] = Obtain_local_networks(NPP_all_grouped_networks, anatomical_class, {columns_Head, columns_Midbody, columns_Tail}); %obtain the short and mid range networks for the all grouped neuropeptide networks
[Neuropeptide_networks_individual_short_range, Neuropeptide_networks_individual_mid_range, class_list_neuropeptide_networks, class_names_neuropeptide_networks] = Obtain_local_networks(neuropeptide_networks_all, anatomical_class, {columns_Head, columns_Midbody, columns_Tail}); %obtain the short and mid range networks for the individual neuropeptide networks

%GROUP NETWORKS
NPP_all_grouped_networks_short_range_connections = Group_local_networks(NPP_all_grouped_networks, {NPP_all_grouped_networks_individual_short_range{2:17}, NPP_all_grouped_networks_individual_short_range{27:35}, NPP_all_grouped_networks_individual_short_range{39:41}}); %group all local networks that have close process bundles
NPP_all_grouped_networks_mid_range_connections = Group_local_networks(NPP_all_grouped_networks, {NPP_all_grouped_networks_individual_mid_range{3}, NPP_all_grouped_networks_individual_short_range{17}, NPP_all_grouped_networks_individual_mid_range{1}, NPP_all_grouped_networks_individual_mid_range{2}}); %group all networks that have short or mid range processes (short_range{15} is the Pharynx, that only has short range processes, the other mid range ones are Head, Tail and Midbody) 

NPP_individual_networks_short_range_connections = cell(size(neuropeptide_networks_all)); %empty cell array to store individual networks matrices
NPP_individual_networks_mid_range_connections = cell(size(neuropeptide_networks_all)); %empty cell array to store individual networks matrices
NPP_individual_networks_nerve_ring_connections = cell(size(neuropeptide_networks_all)); %empty cell array to store individual networks matrices
NPP_individual_networks_pharynx_connections = cell(size(neuropeptide_networks_all)); %empty cell array to store individual networks matrices
 
for i = 1:length(neuropeptide_networks_all)
    NPP_individual_networks_short_range_connections{i,1} = Group_local_networks(neuropeptide_networks_all{i}, {Neuropeptide_networks_individual_short_range{i}{2:17}, Neuropeptide_networks_individual_short_range{i}{27:35}, Neuropeptide_networks_individual_short_range{i}{39:41}}); %select short range networks for all individual neuropeptide networks
    NPP_individual_networks_mid_range_connections{i,1} = Group_local_networks(neuropeptide_networks_all{i}, {Neuropeptide_networks_individual_mid_range{i}{1}, Neuropeptide_networks_individual_short_range{i}{17}, Neuropeptide_networks_individual_mid_range{i}{2}, Neuropeptide_networks_individual_mid_range{i}{3}}); %select mid range networks for all individual neuropeptide networks
    NPP_individual_networks_nerve_ring_connections{i,1} = Group_local_networks(neuropeptide_networks_all{i}, {Neuropeptide_networks_individual_short_range{i}{16}}); %select nerve ring network for all individual neuropeptide networks
    NPP_individual_networks_pharynx_connections{i,1} = Group_local_networks(neuropeptide_networks_all{i}, {Neuropeptide_networks_individual_short_range{i}{17}}); %select pharynx network for all individual neuropeptide networks
end

%PLOT GROUPED NETWORKS
Plot_Networks(NPP_all_grouped_networks_short_range_connections, 'All_NPP_grouped_networks_short_range_connections', ([path 'Images_neuropeptide_networks_unweighted/NPP_all_grouped_networks']), neuron_ID.nodeLabel);%plot and save NPP all grouped networks for neurons that have connections with neurons close to them
Plot_Networks(NPP_all_grouped_networks_mid_range_connections, 'All_NPP_grouped_networks_mid_range_connections', ([path 'Images_neuropeptide_networks_unweighted/NPP_all_grouped_networks']), neuron_ID.nodeLabel);%plot and save NPP all grouped networks for neurons that have connections with neurons close to them


%% ****************************************************************************************************************
%SAVE ALL RESULTS FROM SCRIPT
%****************************************************************************************************************
save((['Results/scRNAseq_C.elegans' path2 date 'Networks_build_EC500nM.mat'])) %save all variables from workspace to results at this point
end
%**************************************************************************
function Network_measures(synA, gapA, Mon_all_grouped_networks, NPP_all_grouped_networks, NPP_all_grouped_networks_mid_range_connections, NPP_all_grouped_networks_short_range_connections, NPP_individual_networks_short_range_connections, path, path2, date, version)

synA = readmatrix('Syn_network.csv'); %We have two datasets of synaptic data (Albertson(AW) and WormWiring (WW) this is the AW with the new classification)
gapA = readmatrix('GJ_network.csv'); %Import the gap junctions table

%BINARISE AND REMOVE SELF LOOPS FROM NETWORK
%networks in binary form and without self loops
synA_SL_u = Mask_network(synA - diag(diag(synA))); %the mask network function binarises the network and the substraction of the diagonal removes self loops
gapA_SL_u = Mask_network(gapA - diag(diag(gapA))); %the mask network function binarises the network and the substraction of the diagonal removes self loops
Mon_grouped_networks_SL_u = Mask_network(Mon_all_grouped_networks - diag(diag(Mon_all_grouped_networks))); %the mask network function binarises the network and the substraction of the diagonal removes self loops
NPP_all_grouped_networks_SL_u = Mask_network (NPP_all_grouped_networks - diag(diag(NPP_all_grouped_networks))); %the mask network function binarises the network and the substraction of the diagonal removes self loops
NPP_all_grouped_networks_mid_range_SL_u = Mask_network (NPP_all_grouped_networks_mid_range_connections - diag(diag(NPP_all_grouped_networks_mid_range_connections))); %the mask network function binarises the network and the substraction of the diagonal removes self loops
NPP_all_grouped_networks_short_range_SL_u = Mask_network (NPP_all_grouped_networks_short_range_connections - diag(diag(NPP_all_grouped_networks_short_range_connections))); %the mask network function binarises the network and the substraction of the diagonal removes self loops
for k = 1:length(NPP_individual_networks_short_range_connections)
    NPP_individual_networks_short_range_connections_SL_u{k,1} = Mask_network(NPP_individual_networks_short_range_connections{k} - diag(diag(NPP_individual_networks_short_range_connections{k}))); %get binarised networks without self loops for each individual receptor_neuropeptide pair
end

%networks in weighted form without self loops 
synA_SL = synA - diag(diag(synA)); %the substraction of the diagonal removes self loops  
gapA_SL = gapA - diag(diag(gapA)); %the substraction of the diagonal removes self loops  
Mon_grouped__SL = Mon_all_grouped_networks- diag(diag(Mon_all_grouped_networks)); %the substraction of the diagonal removes self loops  
NPP_all_grouped_networks_SL = NPP_all_grouped_networks - diag(diag(NPP_all_grouped_networks)); %the substraction of the diagonal removes self loops  
NPP_all_grouped_networks_mid_range_SL = NPP_all_grouped_networks_mid_range_connections - diag(diag(NPP_all_grouped_networks_mid_range_connections)); %the substraction of the diagonal removes self loops  
NPP_all_grouped_networks_short_range_SL = NPP_all_grouped_networks_short_range_connections - diag(diag(NPP_all_grouped_networks_short_range_connections)); %the substraction of the diagonal removes self loops  
for k = 1:length(NPP_individual_networks_short_range_connections)
    NPP_individual_networks_short_range_connections_SL{k,1} = NPP_individual_networks_short_range_connections{k} - diag(diag(NPP_individual_networks_short_range_connections{k})); %get networks without self loops for each individual receptor_neuropeptide pair
end
for k = 1:length(NPP_individual_networks_mid_range_connections)
    NPP_individual_networks_mid_range_connections_SL{k,1} = NPP_individual_networks_mid_range_connections{k} - diag(diag(NPP_individual_networks_mid_range_connections{k})); %get networks without self loops for each individual receptor_neuropeptide pair
end

%CALCULATE RANDOM NETWORKS 
%Unweighted directed and undirected random networks without self loops 
% Number of random networks to generate. Using 100 allows us to get comparisons in percentage, our network being 1 and compared with 100 random ones. 
numNets = 100;
% Number of times to swap each edge to generate the random network (default is 5, 10 is a good number)
numSwaps = 10;

[randSynu, randUnSynu] = RandomNet_building(synA, numNets, numSwaps, 1); %get 100 random networks for each neuropeptide grouped network
[randGapu, randUnGapu] = RandomNet_building(gapA, numNets, numSwaps, 1); %get 100 random networks for each neuropeptide grouped network
[randMonu, randUnMonu] = RandomNet_building(Mon_all_grouped_networks, numNets, numSwaps, 1); %get 100 random networks for each neuropeptide grouped network
[randNPP_NPGPCR_long_rangeu, randUnNPP_NPGPCR_long_rangeu] = RandomNet_building(NPP_all_grouped_networks, numNets, numSwaps, 1); %get 100 random networks for all grouped long range receptor_neuropeptide network
[randNPP_NPGPCR_mid_rangeu, randUnNPP_NPGPCR_mid_rangeu] = RandomNet_building(NPP_all_grouped_networks_mid_range_connections, numNets, numSwaps, 1); %get 100 random networks for all grouped mid range netwokr
[randNPP_NPGPCR_short_rangeu, randUnNPP_NPGPCR_short_rangeu] = RandomNet_building(NPP_all_grouped_networks_short_range_connections, numNets, numSwaps, 1); %get 100 random networks for all grouped short range network 
[randNPP_NPGPCR_all_short_rangeu, randUnNPP_NPGPCR_all_short_rangeu] = RandomNet_building(NPP_individual_networks_short_range_connections_SL, numNets, numSwaps, 1); %get 100 random networks for each individual receptor_neuropeptide network

%****************************************************************************************************************
%DEGREE
%****************************************************************************************************************
%Unweighted: I use the degree_dir.m function from the Brain Connectivity toolbox but I name out and in degree the other way around from what the algorithm states because in my matrix the
%columns are the outdegree (NPP) and the rows the indegree (NPGCPR). For the synaptic and monoamine networks the order is like in the function
%receiver neurons in the rows and releaser/send neurons in the columns. With weighted degree we look at the
%abundance of peptide released or received

[indeg_Syn, outdeg_Syn, degree_Syn] = degrees_dir(synA_SL_u);%use BrainConnectivity function to get directed degrees, indegrees accounts for the synapses the node receives and outdegree for the synapses it forms
Syn_degree_nodes = table(neuron_ID.nodeLabel, transpose(outdeg_Syn), transpose(indeg_Syn), transpose(degree_Syn), 'VariableNames', {'Neuron', 'Outdegree_Syn', 'Indegree_Syn', 'Degree_Syn'}); %transform degree information in table
[degree_GJ] = degrees_und(gapA_SL_u);%use BrainConnectivity function to get undirected degrees,  accounts for the gap junctions in the node 
GJ_degree_nodes = table(neuron_ID.nodeLabel, transpose(degree_GJ), 'VariableNames', {'Neuron', 'Degree_GJ'}); %transform degree information in table
[indeg_Mon, outdeg_Mon, degree_Mon] = degrees_dir(Mon_grouped_networks_SL_u);%use BrainConnectivity function to get directed degrees, indegrees accounts for the monoamines the node receives and outdegree for the monoamines it forms
Mon_degree_nodes = table(neuron_ID.nodeLabel, transpose(outdeg_Mon), transpose(indeg_Mon), transpose(degree_Mon), 'VariableNames', {'Neuron', 'Outdegree_Mon', 'Indegree_Mon', 'Degree_Mon'}); %transform degree information in table

[indeg_NPP_all, outdeg_NPP_all, degree_NPP_all] = degrees_dir(NPP_all_grouped_networks_SL_u); %use BrainConnectivity function to get directed degrees, indegrees accounts for NPP node receives and outdegree for the NPP it releases
NPP_all_degree_nodes = table(neuron_ID.nodeLabel, transpose(outdeg_NPP_all), transpose(indeg_NPP_all), transpose(degree_NPP_all), 'VariableNames', {'Neuron', 'Outdegree_NPP_all', 'Indegree_NPP_all', 'Degree_NPP_all'}); %transform degree information in table
[indeg_NPP_all_mid_range, outdeg_NPP_all_mid_range, degree_NPP_all_mid_range] = degrees_dir(NPP_all_grouped_networks_mid_range_SL_u); %use BrainConnectivity function to get direct degrees, indegrees and outdegrees unweighted
NPP_all_grouped_networks_mid_range_connections_degree_nodes = table(neuron_ID.nodeLabel, transpose(outdeg_NPP_all_mid_range), transpose(indeg_NPP_all_mid_range), transpose(degree_NPP_all_mid_range), 'VariableNames', {'Neuron', 'Outdegree_NPP_all_mid_range', 'Indegree_NPP_all_mid_range', 'Degree_NPP_all_mid_range'}); %transform degree information in table
[indeg_NPP_all_short_range, outdeg_NPP_all_short_range, degree_NPP_all_short_range] = degrees_dir(NPP_all_grouped_networks_short_range_SL_u); %use BrainConnectivity function to get direct degrees, indegrees and outdegrees unweighted
NPP_all_grouped_networks_short_range_connections_degree_nodes = table(neuron_ID.nodeLabel, transpose(outdeg_NPP_all_short_range), transpose(indeg_NPP_all_short_range), transpose(degree_NPP_all_short_range), 'VariableNames', {'Neuron', 'Outdegree_NPP_all_short_range', 'Indegree_NPP_all_short_range', 'Degree_NPP_all_short_range'}); %transform degree information in table

% Connectome_degrees_nodes_unweighted = join(NPP_all_degree_nodes, join(NPP_all_grouped_networks_mid_range_connections_degree_nodes, join(NPP_all_grouped_networks_short_range_connections_degree_nodes, join(FLP_degree_nodes, join(FLP_lessdensereceptors_degree_nodes, join(NLP_degree_nodes, join(Syn_degree_nodes, join(GJ_degree_nodes, Mon_degree_nodes)))))))); %join all tables together into one
Connectome_degrees_nodes_unweighted = join(NPP_all_degree_nodes, join(NPP_all_grouped_networks_mid_range_connections_degree_nodes, join(NPP_all_grouped_networks_short_range_connections_degree_nodes, join(Syn_degree_nodes, join(GJ_degree_nodes, Mon_degree_nodes))))); %join all tables together into one
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Topological_data_unweighted/Topological_data_degree/Excel_tables/']) %create folder to save images
writetable(Connectome_degrees_nodes_unweighted,(['Results/scRNAseq_C.elegans' path 'Topological_data_unweighted/Topological_data_degree/Excel_tables/' date 'scRNAseq_Celegans_connectome_degree_unweighted.csv']),'Delimiter',',','QuoteStrings',true)  %save overall connectome table as a csv file to be able to further analyse in excel

%****************************************************************************************************************
%DENSITY
%****************************************************************************************************************
%Calculates directed density per network. Assumes network has no self loops and doesn't take weight into account. Density is the fraction of present connections to possible connections.
[density_Syn, density_mean_rand_Syn, density_norm_Syn] = normalise_network_measure(synA_SL, randSynu, @density_dir);
[density_GJ, density_mean_rand_GJ, density_norm_GJ] = normalise_network_measure(gapA_SL, randGapu, @density_dir);
[density_Mon, density_mean_rand_Mon, density_norm_Mon] = normalise_network_measure(MonA_SL, randMonu, @density_dir);
[density_NPP_all_grouped_networks, density_mean_rand_All_NPP, density_norm_All_NPP] = normalise_network_measure(NPP_all_grouped_networks, randNPP_NPGPCR_long_rangeu, @density_dir);
[density_NPP_grouped_networks_mid_range, density_mean_rand_NPP_mid_range, density_norm_NPP_mid_range] = normalise_network_measure(NPP_all_grouped_networks_mid_range_SL, randNPP_NPGPCR_mid_rangeu, @density_dir);
[density_NPP_grouped_networks_short_range, density_mean_rand_NPP_short_range, density_norm_NPP_short_range] = normalise_network_measure(NPP_all_grouped_networks_short_range_SL, randNPP_NPGPCR_short_rangeu, @density_dir);

save((['Results/scRNAseq_C.elegans', path2 date 'Networks_measures_EC500nM.mat']))

end
%**************************************************************************
function Network_computeRichClub(A, randA, networkName, savePath, neuronID, equation, type)
    
    if max(max(A)) > 1
        error('! Network must be unweighted !');
    end
    if isDirected(A)
        directed = 1;
    else
        directed = 0;
    end
    
    % Remove self-loops
    A = A - diag(diag(A));

    %**********************************************************************
    %	RICH CLUB COMPUTATION
    %**********************************************************************
	% Calculate rich club coefficient
	% R:	rich-club coefficients for each level
	% Nk:	number of nodes with degree>k
	% Ek:   number of edges remaining in subgraph with degree>k
    if directed
        [phi,Nk,Ek] = equation(A);
    else
        [phi,Nk,Ek] = equation(A);
    end

    %**********************************************************************
    %	NORMALISATION
    %**********************************************************************
    for i=1:length(randA)
		% Calculate rich club metric for each random graph
        if directed
            phi_rand(i,1:length(phi)) = equation(randA{i});
        else
            phi_rand(i,1:length(phi)) = equation(randA{i});
        end
    end
    
    % Get mean and stdev for each value of k over all random networks
	phi_rand_mean = nanmean(phi_rand);
	phi_rand_std = nanstd(phi_rand);
	
	% Normalise phi
	phi_norm = phi ./ phi_rand_mean;
    
    %**********************************************************************
    %	PLOT RICH CLUB CURVE
    %**********************************************************************
	f1 = figure();
    hLine1 = plot(phi, 'b', 'LineWidth', 2);
    hold on;
    hLine2 = plot(phi_rand_mean, 'r', 'LineWidth', 2);
    errorbar(1:5:length(phi_rand_mean), phi_rand_mean(1:5:end), ...
                phi_rand_std(1:5:end), '.r', 'LineWidth', 2);
    hLine3 = plot(phi_norm, 'Color', 'g', 'LineWidth', 2);
    errorbar(1:5:length(phi_rand_mean), phi_norm(1:5:end), ...
                phi_rand_std(1:5:end), '.g', 'LineWidth', 2);
    xlabel('Degree (k)','FontSize',25);
    ylabel('Rich-club coefficient (\Phi)','FontSize',25);
    set(gca,'FontSize',25,'LineWidth',1.25,'XLim', [0 length(phi)]);
    h_legend = legend([hLine1 hLine2 hLine3], ...
            '\Phi(k)_{C. elegans}','\Phi(k)_{random}','\Phi(k)_{norm}',...
            'Location', 'northwest');
    set(h_legend,'FontSize',25, 'Location', 'southeast');
    ylabh = get(gca,'YLabel');
    set(ylabh,'Position',get(ylabh,'Position') - [2 0 0]);

    %**********************************************************************
    %	CALCULATE RICH CLUB REGIMES
    %**********************************************************************
  	% findRegime : local function -- see end of file
    % Check for 1 - 10 sigma
    for i=1:10
        sigma{i} = findRegime(phi_norm, phi_rand_std, i);
    end
    % Get neuron degrees
    if directed
        [~,~,deg] = degrees_dir(A);
        [indeg,~,~] = degrees_dir(A);
        [~,outdeg,~] = degrees_dir(A);
    else
        deg = degrees_und(A);
    end
    
    if isequal(type, 'deg')
        [deg,indx] = sort(deg,'descend');
    elseif isequal(type, 'indeg')
        [indeg,indx] = sort(indeg,'descend');
    else isequal(type, 'outdeg')
        [outdeg,indx] = sort(outdeg,'descend');
    end
    %**********************************************************************
    %	SAVE RESULTS
    %**********************************************************************
    % Create save directory
    mkdir(savePath);
    
    % SAVE FIGURE
    saveas(f1, strcat(savePath,'/Fig_RichClub_',networkName),'png');
    close(f1);
    
    % SAVE STATISTICS
    fid = fopen(strcat(savePath, '/Stats_RichClub_',...
                    networkName,'.csv'),'wt','native');
    fprintf(fid,'Degree,Phi_norm,Phi,Phi_rand,Stdev\n');
    for i=1:length(phi)
        fprintf(fid,'%f,%f,%f,%f,%f\n',...
            i, phi_norm(i), phi(i), phi_rand_mean(i), phi_rand_std(i));
    end
    fclose(fid);
    
    % SAVE NEURONS
    % Get degrees of neurons
    % Output a separate file for each sigma level
    for i=1:length(sigma)
        fid = fopen(strcat(savePath, '/',num2str(i),...
                    '_Sigma_RC_',networkName,'.txt'),'wt','native');
        % Iterate over regimes and get neuron names
        for reg=1:length(sigma{i})
            fprintf(fid,'\n*******************************************\n');
            fprintf(fid,'REGIME %i -- (degree %i to %i)\n', reg, ...
                            sigma{i}{reg}(1), sigma{i}{reg}(2));
            fprintf(fid,'*******************************************\n');
            % Output names of neurons for the last regime
            if isequal(type,'deg')
                if reg==length(sigma{i})
                neurons = find(deg>=sigma{i}{reg}(1));
                neurons = indx(neurons);
                for n=1:length(neurons)
                    fprintf(fid,'%s (%i)\n',...
                    char(neuronID(neurons(n))), int32(deg(n)));
                end
                end
             elseif isequal(type,'indeg')
                if reg==length(sigma{i})
                neurons = find(indeg>=sigma{i}{reg}(1));
                neurons = indx(neurons);
                for n=1:length(neurons)
                    fprintf(fid,'%s (%i)\n',...
                    char(neuronID(neurons(n))), int32(indeg(n)));
                end
                end
              else isequal(type,'outdeg')
                if reg==length(sigma{i})
                neurons = find(outdeg>=sigma{i}{reg}(1));
                neurons = indx(neurons);
                for n=1:length(neurons)
                    fprintf(fid,'%s (%i)\n',...
                    char(neuronID(neurons(n))), int32(outdeg(n)));
                end
                end
            end
        end
        fclose(fid);
    end
    %**********************************************************************
    %	SAVE DATA 
    %**********************************************************************
    % Create save table
    
    Phi_values_table = table([1:length(phi)].', phi.', phi_rand_mean.', phi_norm.', 'VariableNames',{'Degree', 'C.elegans', 'Random', 'Norm'}); %join all tables together into one
    writetable(Phi_values_table,([savePath '/' networkName '.csv']),'Delimiter',',','QuoteStrings',true);  %save overall connectome table as a csv file to be able to further analyse in excel
end

function [R,Nk,Ek] = rich_club_bd_in(CIJ,varargin)
%RICH_CLUB_BD_IN        Rich club coefficients (binary directed graph)
%
%   R = rich_club_bd(CIJ)
%   [R,Nk,Ek] = rich_club_bd(CIJ,klevel)
%
%   The rich club coefficient, R, at level k is the fraction of edges that
%   connect nodes of indegree k or higher out of the maximum number of edges
%   that such nodes might share.
%
%   Input:      CIJ,        connection matrix, binary and directed
%            klevel,        optional input argument. klevel sets the
%                              maximum level at which the rich club
%                              coefficient will be calculated. If klevel is
%                              not included the the maximum level will be
%                              set to the maximum degree of CIJ.
%
%   Output:       R,        vector of rich-club coefficients for levels
%                              1 to klevel.
%                Nk,        number of nodes with indegree>k
%                Ek,        number of edges remaining in subgraph with
%                              indegree>k

N = size(CIJ,1);                    %#ok<NASGU>

% definition of "indegree" as used for RC coefficients
% indegree is taken to be the sum of incoming and outgoing connectons
[indegree,~,~] = degrees_dir(CIJ);

if nargin == 1
    klevel = max(indegree);
elseif nargin == 2
    klevel = varargin{1};
elseif nargin > 2
    error('number of inputs incorrect. Should be [CIJ], or [CIJ, klevel]')
end

R = zeros(1,klevel);
Nk = zeros(1,klevel);
Ek = zeros(1,klevel);
for k = 1:klevel
    SmallNodes=find(indegree<=k);       %get 'small nodes' with indegree <=k
    subCIJ=CIJ;                       %extract subnetwork of nodes >k by removing nodes <=k of CIJ
    subCIJ(SmallNodes,:)=[];          %remove rows
    subCIJ(:,SmallNodes)=[];          %remove columns
    Nk(k)=size(subCIJ,2);             %number of nodes with indegree >k
    Ek(k)=sum(subCIJ(:));             %total number of connections in subgraph
    R(k)=Ek(k)/(Nk(k)*(Nk(k)-1));     %unweighted rich-club coefficient
end
end

function [R,Nk,Ek] = rich_club_bd_out(CIJ,varargin)
%RICH_CLUB_BD_Out        Rich club coefficients (binary directed graph)
%
%   R = rich_club_bd(CIJ)
%   [R,Nk,Ek] = rich_club_bd(CIJ,klevel)
%
%   The rich club coefficient, R, at level k is the fraction of edges that
%   connect nodes of outdegree k or higher out of the maximum number of edges
%   that such nodes might share.
%
%   Input:      CIJ,        connection matrix, binary and directed
%            klevel,        optional input argument. klevel sets the
%                              maximum level at which the rich club
%                              coefficient will be calculated. If klevel is
%                              not included the the maximum level will be
%                              set to the maximum degree of CIJ.
%
%   Output:       R,        vector of rich-club coefficients for levels
%                              1 to klevel.
%                Nk,        number of nodes with outdegree>k
%                Ek,        number of edges remaining in subgraph with
%                              outdegree>k

N = size(CIJ,1);                    %#ok<NASGU>

% definition of "outdegree" as used for RC coefficients
% outdegree is taken to be the sum of incoming and outgoing connectons
[~,outdegree,~] = degrees_dir(CIJ);

if nargin == 1
    klevel = max(outdegree);
elseif nargin == 2
    klevel = varargin{1};
elseif nargin > 2
    error('number of inputs incorrect. Should be [CIJ], or [CIJ, klevel]')
end

R = zeros(1,klevel);
Nk = zeros(1,klevel);
Ek = zeros(1,klevel);
for k = 1:klevel
    SmallNodes=find(outdegree<=k);       %get 'small nodes' with outdegree <=k
    subCIJ=CIJ;                       %extract subnetwork of nodes >k by removing nodes <=k of CIJ
    subCIJ(SmallNodes,:)=[];          %remove rows
    subCIJ(:,SmallNodes)=[];          %remove columns
    Nk(k)=size(subCIJ,2);             %number of nodes with outdegree >k
    Ek(k)=sum(subCIJ(:));             %total number of connections in subgraph
    R(k)=Ek(k)/(Nk(k)*(Nk(k)-1));     %unweighted rich-club coefficient
end
end

function regime = findRegime(phi_norm, phi_rand_std, sigma)

% FINDREGIME returns the start and end degree of each rich club regime
	regime = cell.empty(0,1);
    
    newRegime = 1;
    curRegime = 0;
    
    % Iterate over each degree level k
	for k=1:length(phi_norm)
        % Check for rich club at k
		if phi_norm(k) >= (1 + (sigma * phi_rand_std(k))) 
            % If this is a new regime add start degree to 
            % regime(curRegime)(1) else if it is a continuation, 
            % add as the last degree to (2)
            if newRegime
                curRegime = curRegime + 1;
                regime{curRegime}(1) = k;
                newRegime = 0;
                regime{curRegime}(2) = k;
            else
                regime{curRegime}(2) = k;
            end
        else
            % End of regime. Increment counters / set flags.
            newRegime = 1;
		end
	end
end
%**************************************************************************
function Network_dimensionality_reduction(long_range_network, mid_range_network, short_range_network, indeg_mid, individual_network_mid_SL, neurons, class, names, path, path2, date, version)

    %create directory to save plots 
    mkdir ([ 'Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted'])  %create folder to save images
    %set defaults
    set(0, 'DefaultAxesFontName', 'Arial'); %set default text
    
    %% ****************************************************************************************************************
    %NORMALISE DATA
    %****************************************************************************************************************
    
    %normalise the matrix columns using z-score, it scales (makes standard deviation 1,
    %which is essential for t-sne and not a problem for PCA) and centers 
    %(makes column means 0, which is needed for PCA and not a problem for
    %t-sne). For the incoming connections analysis we transpose first  
    
    NPP_all_grouped_networks_normal_out = zscore(long_range_network,1,1); %get the zscore in the columns dimension, using population standard deviation as we consider the incoming connections of that neuron the population
    NPP_all_grouped_networks_mid_range_normal_out = zscore(mid_range_network,1,1); %get the zscore in the columns dimension, using population standard deviation as we consider the incoming connections of that neuron the population
    NPP_all_grouped_networks_short_range_normal_out = zscore(short_range_network,1,1); %get the zscore in the columns dimension, using population standard deviation as we consider the incoming connections of that neuron the population
    
    NPP_all_grouped_networks_normal_in = zscore(transpose(long_range_network),1,1); %get the zscore in the columns dimension, using population standard deviation as we consider the incoming connections of that neuron the population
    NPP_all_grouped_networks_mid_range_normal_in = zscore(transpose(mid_range_network),1,1); %get the zscore in the columns dimension, using population standard deviation as we consider the incoming connections of that neuron the population
    NPP_all_grouped_networks_short_range_normal_in = zscore(transpose(short_range_network),1,1); %get the zscore in the columns dimension, using population standard deviation as we consider the incoming connections of that neuron the population
    
    % %% ****************************************************************************************************************
    %RUN ANALYSIS
    % %****************************************************************************************************************
    
    %PCA
    [PCA_coeff_NPP_all2, PCA_scores_NPP_all2, PCA_latent_NPP_all2, PCA_tsquared_NPP_all2, PCA_explained_NPP_all2, PCA_mu_NPP_all2] = pca(NPP_all_grouped_networks_normal_in, 'Centered', 'off'); %run PCA equation on the normalised 
    [PCA_coeff_NPP_all_mid2, PCA_scores_NPP_all_mid2, PCA_latent_NPP_all_mid2, PCA_tsquared_NPP_all_mid2, PCA_explained_NPP_all_mid2, PCA_mu_NPP_all_mid2] = pca(NPP_all_grouped_networks_mid_range_normal_in, 'Centered', 'off'); %run PCA equation on the normalised data
    [PCA_coeff_NPP_all_short2, PCA_scores_NPP_all_short2, PCA_latent_NPP_all_short2, PCA_tsquared_NPP_all_short2, PCA_explained_NPP_all_short2, PCA_mu_NPP_all_short2] = pca(NPP_all_grouped_networks_short_range_normal_in, 'Centered', 'off'); %run PCA equation on the normalised data
    
    %t-SNE
    rng('default')% set a seed to 0 so that the plots are comparable 
    tnse_NPP_all_2 = tsne(NPP_all_grouped_networks_normal_in,'Algorithm','exact','Distance','euclidean');%t-sne based on columns 
    rng('default')% set a seed to 0 so that the plots are comparable 
    tnse_NPP_all_mid_2 = tsne(NPP_all_grouped_networks_mid_range_normal_in,'Algorithm','exact','Distance','euclidean');%t-sne based on columns 
    rng('default')% set a seed to 0 so that the plots are comparable 
    tnse_NPP_all_short_2 = tsne(NPP_all_grouped_networks_short_range_normal_in,'Algorithm','exact','Distance','euclidean'); %t-sne based on columns 
    
    %define colour maps based on the mid range network for t-sne (best at
    %splitting groups)
    %get the neurons that are in each PCA and t-sne group 

    [tsne_NPP_all_mid_group_neurons, tsne_NPP_all_mid_index_group] = sortrows(tnse_NPP_all_mid_2, 1); %sort rows of the tsne second variable from smaller to bigger in order to identify modules but keep index in order to identify neurons 
    tsne_NPP_all_mid_group_neurons_table = table (tsne_NPP_all_mid_group_neurons(:,1),tsne_NPP_all_mid_group_neurons(:,2), neurons(tsne_NPP_all_mid_index_group)); %make a table with the name of the neurons next to the t-sne distance values
    for i = 1: length(tsne_NPP_all_mid_group_neurons_table.Var1) %loop over every neuron in the network 
        if tsne_NPP_all_mid_group_neurons_table.Var1(i) < -10 && tsne_NPP_all_mid_group_neurons_table.Var2(i) < -5 %select indices of neurons based on component 1 of t-sne, if larger than 0 module 1
            tsne_NPP_all_mid_group_1(i) = tsne_NPP_all_mid_index_group(i);
        elseif tsne_NPP_all_mid_group_neurons_table.Var1(i) > 16 && tsne_NPP_all_mid_group_neurons_table.Var2(i) < -5 %select indices of neurons if component 2 of t-sne is smaller than 0 and component 1 is smaller than 0 
            tsne_NPP_all_mid_group_2(i) = tsne_NPP_all_mid_index_group(i);
        elseif tsne_NPP_all_mid_group_neurons_table.Var1(i) > 0 && tsne_NPP_all_mid_group_neurons_table.Var1(i) < 16 && tsne_NPP_all_mid_group_neurons_table.Var2(i) < 0 %select indices of neurons based on component 1 of t-sne, if larger than 0 and component 2 < 0 module 2
            tsne_NPP_all_mid_group_3(i) = tsne_NPP_all_mid_index_group(i);
        end
    end

    for i = 1:length(neurons) 
        if ismember(i, tsne_NPP_all_mid_group_1)
            tsne_NPP_all_mid_attributes_group{i} = char('group1');
            tsne_NPP_all_mid_color_group(i,:) = [0.4660, 0.6740, 0.1880];%colormap for the groups
        elseif ismember(i, tsne_NPP_all_mid_group_2)
            tsne_NPP_all_mid_attributes_group{i} = char('group2');
            tsne_NPP_all_mid_color_group(i,:) = [0.9290, 0.6940, 0.1250];%colormap for the groups
        elseif ismember(i,tsne_NPP_all_mid_group_3)
            tsne_NPP_all_mid_attributes_group{i} = char('group3');
            tsne_NPP_all_mid_color_group(i,:) = [0.4940, 0.1840, 0.5560];%colormap for the groups
        else
            tsne_NPP_all_mid_attributes_group{i} = char('none'); 
            tsne_NPP_all_mid_color_group(i,:) = [0.5 0.5 0.5];%colormap for the groups
        end
    end
    
    %Plot
    plot_PCA_texture(PCA_scores_NPP_all2, tsne_NPP_all_mid_color_group, neurons, date, path, 'NPP_all_grouped_long_range_normalised_incoming_connections', 'on', class) %plot PCA 
    plot_PCA_texture(PCA_scores_NPP_all_short2, tsne_NPP_all_mid_color_group, neurons, date, path, 'NPP_all_grouped_short_range_normalised_incoming_connections_2', 'on', class) %plot PCA 
    plot_PCA_texture(PCA_scores_NPP_all_mid2, tsne_NPP_all_mid_color_group, neurons, date, path, 'NPP_all_grouped_mid_range_normalised_incoming_connections', 'on', class) %plot PCA 
    
    plot_tsne_texture(tnse_NPP_all_2, tsne_NPP_all_mid_color_group, neurons, date, path, 'NPP_all_grouped_networks_long_range_incoming_connections', 'on', class) %plot t-sne based on the local function
    plot_tsne_texture(tnse_NPP_all_short_2, tsne_NPP_all_mid_color_group, neurons, date, path, 'NPP_all_grouped_networks_short_range_incoming_connections_2', 'on', class) %plot t-sne based on the local function
    plot_tsne_texture(tnse_NPP_all_mid_2, tsne_NPP_all_mid_color_group, neurons, date, path, 'NPP_all_grouped_networks_mid_range_incoming_connections', 'on', class) %plot t-sne based on the local function
    
    %plot the loadings of PCA to see which neurons affect more the clustering 
    image_biplot = figure('Visible', 'off'); %create figure to save
    biplot(PCA_coeff_NPP_all_mid2(:,1:2), 'Varlabel', neurons) %plot loadings the longer the line the more it affects
    hold on 
    scatter(PCA_coeff_NPP_all_mid2(:,1), PCA_coeff_NPP_all_mid2(:,2), [], tsne_NPP_all_mid_color_group, 'filled') %plot neurons coloured by t-sne group
    title('NPP_all_grouped_mid_range_normalised_incoming_connections_PCA_loadings', 'Interpreter', 'none');
    saveas(image_biplot, ['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date 'NPP_all_grouped_mid_range_normalised_incoming_connections_PCA_loadings_biplot.png']);%save image as pdf

    %reorder and colour spy matrix to highlight the groups
    [tsne_NPP_all_mid_attributes_group_sorted, tsne_NPP_all_mid_range_connections_group_reordered_index] = sort(tsne_NPP_all_mid_attributes_group); %set index to sort matrix based on groups
    tsne_NPP_all_mid_range_connections_group_reordered_2 = mid_range_network(tsne_NPP_all_mid_range_connections_group_reordered_index.', tsne_NPP_all_mid_range_connections_group_reordered_index); %sort matrix based on groups

    %plot cspy of reordered matrix 
    image_cspy = figure ('Visible', 'off', 'Position', [109.8035,94.38999999999999,648.6362499999999,691.935]);
    cspy_Lidia(tsne_NPP_all_mid_range_connections_group_reordered_2, 'markersize', 5, 'ColorMap', 'jet', 'Levels', [5]); %plots a spy plot but colours the edges by regions of weight(that is why we select level 5, so that there are 5 colours, weights are split in 5 categories, we choose 5 because we saw it was a good number for the data) and has colorbar
    xlabel ('Receiving neuron')
    ylabel ('Sensing neuron')
    axis('off')
    colorbar;
    saveas(image_cspy, ['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date 'NPP_all_grouped_networks_mid_range_incoming_connections_reordered_by_group_t-sne_cspy.pdf']);%save image as pdf

    %plot correlations with degree and group  
    image_corr = figure('Visible', 'off');%create figure to plot
    violinplot(indeg_mid, tsne_NPP_all_mid_attributes_group, 'ViolinColor', [0.4940, 0.1840, 0.5560] , 'ShowData', false); %vioplinplot degrees of neurons in each module
    ylabel('Neuropeptide mid range indegree'); %label the y axes
    saveas(image_corr, ['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date 'NPP_all_grouped_networks_mid_range_incoming_connections_degree_correlation_groups_violinplot_purple.pdf']);%save image as pdf

    [violinplot_indegree_p, ~, violinplot_indegree_stats] = kruskalwallis(indeg_mid, tsne_NPP_all_mid_attributes_group); %calculate correlation significance between groups
    violinplot_indegree_pairwise = multcompare(violinplot_indegree_stats, 'CType', 'dunn-sidak'); %correct for multiple comparison
    
    %bar plot of the NPP-GPCR pairs that have more connections on the neurons
    %of each group
    for i = 1:length(individual_network_mid_SL) %loop over every neuron 
        [indeg_NPP_individual_mid_range(i,:), ~,~] = degrees_dir(individual_network_mid_SL{i}); %indegree of mid range connections by NPP_GPCR pairs
    end
    
    indeg_NPP_individual_mid_range_group1 = indeg_NPP_individual_mid_range(:, strcmp(tsne_NPP_all_mid_attributes_group, 'group1')); %select indegree by network for group1 neurons
    indeg_NPP_individual_mid_range_group2 = indeg_NPP_individual_mid_range(:, strcmp(tsne_NPP_all_mid_attributes_group, 'group2')); %select indegree by network for group1 neurons
    indeg_NPP_individual_mid_range_group3 = indeg_NPP_individual_mid_range(:, strcmp(tsne_NPP_all_mid_attributes_group, 'group3')); %select indegree by network for group1 neurons
    
    indeg_NPP_mid_range_group1 = sum(indeg_NPP_individual_mid_range_group1.'); %sum all the incoming connections by NPP-GPCR network
    [indeg_NPP_individual_mid_range_group1_sorted_descend, indeg_NPP_individual_mid_range_group1_sorted_index] = sort(indeg_NPP_mid_range_group1, 'descend'); %sort connections from network that has more to less
    pair_names_NPP_sorted_group1 = categorical(indeg_NPP_individual_mid_range_group1_sorted_index, [1:height(names)], names); %category of sorted NPP in the order stated by index_NPP_to_GPCR after reordering
    
    indeg_NPP_mid_range_group2 = sum(indeg_NPP_individual_mid_range_group2.'); %sum all the incoming connections by NPP-GPCR network
    [indeg_NPP_individual_mid_range_group2_sorted_descend, indeg_NPP_individual_mid_range_group2_sorted_index] = sort(indeg_NPP_mid_range_group2, 'descend'); %sort connections from network that has more to less
    pair_names_NPP_sorted_group2 = categorical(indeg_NPP_individual_mid_range_group2_sorted_index, [1:height(names)], names); %category of sorted NPP in the order stated by index_NPP_to_GPCR after reordering
    
    indeg_NPP_mid_range_group3 = sum(indeg_NPP_individual_mid_range_group3.'); %sum all the incoming connections by NPP-GPCR network
    [indeg_NPP_individual_mid_range_group3_sorted_descend, indeg_NPP_individual_mid_range_group3_sorted_index] = sort(indeg_NPP_mid_range_group3, 'descend'); %sort connections from network that has more to less
    pair_names_NPP_sorted_group3 = categorical(indeg_NPP_individual_mid_range_group3_sorted_index, [1:height(names)], names); %category of sorted NPP in the order stated by index_NPP_to_GPCR after reordering
    
    image_bar = figure('Visible', 'off', 'Position', [862,103,621,859]);%create figure to plot
    t = tiledlayout(3,1);
    t.TileSpacing = 'loose';
    nexttile
    bar(categorical(pair_names_NPP_sorted_group1(1:10), pair_names_NPP_sorted_group1(1:10)), indeg_NPP_individual_mid_range_group1_sorted_descend(1:10), 'FaceColor', [0.4660, 0.6740, 0.1880]) %plot top 20 in bar plot
    set(gca,'TickLabelInterpreter','none', 'TickLength', [0, 0], 'Fontsize', 14)%prevent text interpreter, remove axis ticks and set font size
    nexttile
    bar(categorical(pair_names_NPP_sorted_group2(1:10), pair_names_NPP_sorted_group2(1:10)), indeg_NPP_individual_mid_range_group2_sorted_descend(1:10), 'FaceColor', [0.9290, 0.6940, 0.1250])%plot top 20 in bar plot
    set(gca,'TickLabelInterpreter','none', 'TickLength', [0, 0], 'Fontsize', 14)%prevent text interpreter, remove axis ticks and set font size
    nexttile
    bar(categorical(pair_names_NPP_sorted_group3(1:10), pair_names_NPP_sorted_group3(1:10)), indeg_NPP_individual_mid_range_group3_sorted_descend(1:10), 'FaceColor', [0.4940, 0.1840, 0.5560])%plot top 20 in bar plot
    set(gca,'TickLabelInterpreter','none', 'TickLength', [0, 0], 'Fontsize', 14) %prevent text interpreter, remove axis ticks and set font size
    saveas(image_bar, ['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date 'NPP_all_grouped_networks_mid_range_incoming_connections_by_NPP-GPCR_network.pdf']);%save image as pdf
    
    % %% ****************************************************************************************************************
    %LOCAL NETWORKS
    % %****************************************************************************************************************
    
    function plot_PCA_texture(scores, PCA_color3D, neurons, date, path, name, texture, class)
        if strcmp(texture, 'off') == 1
            image3 = figure('Visible', 'off'); %create figure to plot and select the size
            scatter (scores(:, 1), scores(:, 2), [], PCA_color3D, 'o', 'filled'); %plot t-sne coordinates or points for each connection per neuron in dataset
        %     text(tsne_data(:,1)+0.2, tsne_data(:,2)+0.2, neurons); %plot name of each connection next to each dot
        %     title (name, 'interpreter', 'none'); %adds title to font size 13
            axis equal;
            saveas(image3, ['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date name '_PCA_scatter_2D_no_names.pdf']);%save image as pdf
        elseif strcmp(texture, 'on') == 1
            image3 = figure('Visible', 'off'); %create figure to plot and select the size
            tsne_data_s1 = scores(strcmp(class, 'Pharynx'), 1:2); 
            color_s1 = PCA_color3D(strcmp(class, 'Pharynx'), :);
            s1 = scatter (tsne_data_s1(:, 1), tsne_data_s1(:, 2), [], color_s1, '^', 'filled'); %plot t-sne coordinates or points for each connection per neuron in dataset
            hold on 
            tsne_data_s2 = scores(strcmp(class, 'sensory neuron'), 1:2); 
            color_s2 = PCA_color3D(strcmp(class, 'sensory neuron'), :);
            s2 = scatter (tsne_data_s2(:, 1), tsne_data_s2(:, 2), [], color_s2, 'd', 'filled'); %plot t-sne coordinates or points for each connection per neuron in dataset
            hold on 
            tsne_data_s3 = scores(strcmp(class, 'interneuron'), 1:2); 
            color_s3 = PCA_color3D(strcmp(class, 'interneuron'), :);
            s3 = scatter (tsne_data_s3(:, 1), tsne_data_s3(:, 2), [], color_s3, 's', 'filled'); %plot t-sne coordinates or points for each connection per neuron in dataset
            hold on
            tsne_data_s4 = scores(strcmp(class, 'motor neuron'), 1:2); 
            color_s4 = PCA_color3D(strcmp(class, 'motor neuron'), :);
            s4 = scatter (tsne_data_s4(:, 1), tsne_data_s4(:, 2), [], color_s4, 'o', 'filled'); %plot t-sne coordinates or points for each connection per neuron in dataset
            hold off 
            axis equal;
            saveas(image3, ['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date name '_PCA_scatter_2D_texture_no_names.pdf']);%save image as pdf
        end
    end 
 
    function plot_tsne_texture(tsne_data, color_scheme, neurons, date, path, name, texture, class) 
        if strcmp(texture, 'off') == 1
            image4 = figure('Visible', 'off'); %create figure to plot and select the size
            scatter (tsne_data(:, 1), tsne_data(:, 2), [], color_scheme, 'o', 'filled'); %plot t-sne coordinates or points for each connection per neuron in dataset
        %     text(tsne_data(:,1)+0.2, tsne_data(:,2)+0.2, neurons); %plot name of each connection next to each dot
        %     title (name, 'interpreter', 'none'); %adds title to font size 13
            axis equal;
            saveas(image4, ['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date name '_t-sne_scatter_2D_no_names.pdf']);%save image as png
        elseif strcmp(texture, 'on') == 1
            image4 = figure('Visible', 'off'); %create figure to plot and select the size
            tsne_data_s1 = tsne_data(strcmp(class, 'Pharynx'), :); 
            color_s1 = color_scheme(strcmp(class, 'Pharynx'), :);
            s1 = scatter (tsne_data_s1(:, 1), tsne_data_s1(:, 2), [], color_s1, '^', 'filled'); %plot t-sne coordinates or points for each connection per neuron in dataset
            hold on 
            tsne_data_s2 = tsne_data(strcmp(class, 'sensory neuron'), :); 
            color_s2 = color_scheme(strcmp(class, 'sensory neuron'), :);
            s2 = scatter (tsne_data_s2(:, 1), tsne_data_s2(:, 2), [], color_s2, 'd', 'filled'); %plot t-sne coordinates or points for each connection per neuron in dataset
            hold on 
            tsne_data_s3 = tsne_data(strcmp(class, 'interneuron'), :); 
            color_s3 = color_scheme(strcmp(class, 'interneuron'), :);
            s3 = scatter (tsne_data_s3(:, 1), tsne_data_s3(:, 2), [], color_s3, 's', 'filled'); %plot t-sne coordinates or points for each connection per neuron in dataset
            hold on
            tsne_data_s4 = tsne_data(strcmp(class, 'motor neuron'), :); 
            color_s4 = color_scheme(strcmp(class, 'motor neuron'), :);
            s4 = scatter (tsne_data_s4(:, 1), tsne_data_s4(:, 2), [], color_s4, 'o', 'filled'); %plot t-sne coordinates or points for each connection per neuron in dataset
            hold off 
            axis equal;
            saveas(image4, ['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date name '_t-sne_scatter_2D_texture_no_names.pdf']);%save image as png
        end
    end
    
    function plot_scatter_expression(data_NPP, data_GPCR, exclussion_criteria1, exclussion_criteria2, exclussion_criteria3, yaxesnames1, yaxesnames2, date, path, neurons, name)
        
        [y_NPP, x_NPP ] = find(data_NPP); %get x and y coordinates to plot scatter
        [y_GPCR, x_GPCR ] = find(data_GPCR); %get x and y coordinates to plot scatter
        
        %get colormap for the NPP plot
        for i = 1: length(x_NPP)
            if ismember(x_NPP(i), exclussion_criteria1) 
               NPP_colormap(i,:) = [0.4660, 0.6740, 0.1880];
            elseif ismember(x_NPP(i), exclussion_criteria2)
               NPP_colormap(i,:) = [0.9290, 0.6940, 0.1250];
            elseif ismember(x_NPP(i), exclussion_criteria3)
               NPP_colormap(i,:) = [0.4940, 0.1840, 0.5560];
            else
               NPP_colormap(i,:) = [0.5 0.5 0.5]; 
            end
        end
    
        %get colormap for the GPCR plot
        for i = 1: length(x_GPCR)
            if ismember(x_GPCR(i), exclussion_criteria1) 
               GPCR_colormap(i,:) = [0.4660, 0.6740, 0.1880];
            elseif ismember(x_GPCR(i), exclussion_criteria2) 
               GPCR_colormap(i,:) = [0.9290, 0.6940, 0.1250];
            elseif ismember(x_GPCR(i), exclussion_criteria3)
               GPCR_colormap(i,:) = [0.4940, 0.1840, 0.5560];
            else
               GPCR_colormap(i,:) = [0.5 0.5 0.5];
            end
        end
        
        %plot both datasets as two scatter plots in one figure    
        image_expression = figure('Visible', 'off', 'Resize', 'off', 'Position', [199.625,350,3600,2204]); %create an image to save the figure with the right size and prevent it being resize by the plots
        t = tiledlayout(2,1);
        nexttile 
        scatter(x_NPP, y_NPP, [], NPP_colormap, 'filled'); %plot scatter of short range (scatter instead of spy to choose colour)
        yticks( linspace(1, height(data_NPP), height(data_NPP))); %add all y axis ticks
        yticklabels(yaxesnames1); %labels for y axes
        set(gca, 'TickLength', [0, 0], 'TickLabelInterpreter','none', 'Fontsize', 7); %remove tick marks and prevent unwanted underscores in tick labels
        ylabel('NPP expression') %y axes label
        nexttile
        scatter(x_GPCR, y_GPCR, [], GPCR_colormap, 'filled'); %plot scatter of short range (scatter instead of spy to choose colour)
        yticks( linspace(1, height(data_GPCR), height(data_GPCR))); %add all y axis ticks
        yticklabels(yaxesnames2); %labels for y axes
        set(gca, 'TickLength', [0, 0], 'TickLabelInterpreter','none', 'Fontsize', 7); %remove tick marks and prevent unwanted underscores in tick labels
        ylabel('GPCR expression') %y axes label
        xticks(linspace(1, 302, 302)); %add all x axis ticks
        xticklabels(neurons); %labels for x axes
        xlabel(t, 'Neurons') %x axes label
        xtickangle (90); %change orientation of labels
        set(gca, 'TickLength', [0, 0], 'TickLabelInterpreter','none', 'Fontsize', 7); %remove tick marks and prevent unwanted underscores in tick labels
        saveas(image_expression, ['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date name '_coloured_tsne_group_scatter.png']);%save image as png
    end
    %% ****************************************************************************************************************
    %SAVE ALL RESULTS FROM SCRIPT
    %****************************************************************************************************************
    save((['Results/scRNAseq_C.elegans' path2 date 'Networks_modules_EC500nM.mat'])) %save all variables from workspace to results at this point

end
%**************************************************************************
function Network_models_comparison(connections, pair_names, long_range_network, mid_range_network, short_range_network, short_individual, NR_connections, ...
    mid_individual, long_individual, neurons, class_Head, class_Midbody, class_Tail, path, path2, date, version)

%load data
NR_strata_class = readtable('21032021_Colon-Ramos_nerve_ring_neuron_strata.csv');%import list of nerve ring strata class per neuron
%load([date 'Network_measures_EC500nM.mat'])

%create directory to save plots 
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Comparison_btw_ranges'])  %create folder to save images

%set defaults 
set(0, 'DefaultAxesFontName', 'Arial'); %set default text

%% ****************************************************************************************************************
%CREATE MATRICES WITH CONNECTIONS PER NPP-GPCR PAIR FOR NR
% ****************************************************************************************************************

%Get the adjacency matrix for the number of times each NPP_GPCR pair makes
%a connection for each nerve ring and pharynx neuron and plot it
[NPP_GPCR_connections_NR_Amatrix, NPP_GPCR_NR_Amatrix_rows] = Connections_Amatrix (connections, NR_connections, neurons, pair_names,...
    'Number of connections in nerve ring made by each of the individual NPP-GPCR pairs for each neuron node (91 short range)', [ path 'Comparison_btw_ranges'], date);

%Get the adjancency matrix for the number of times each NPP_GPCR pair makes
%a connection for each nerve ring stratum short range network plot it

NPP_individual_networks_nerve_ring_connections_strata_1 = cell(size(NR_connections)); %empty cell array to store matrices
NPP_individual_networks_nerve_ring_connections_strata_2 = cell(size(NR_connections)); %empty cell array to store matrices
NPP_individual_networks_nerve_ring_connections_strata_3 = cell(size(NR_connections)); %empty cell array to store matrices
NPP_individual_networks_nerve_ring_connections_strata_4 = cell(size(NR_connections)); %empty cell array to store matrices
NPP_individual_networks_nerve_ring_connections_strata_5 = cell(size(NR_connections)); %empty cell array to store matrices

NPP_individual_networks_nerve_ring_connections_strata_1(:,1) = {zeros(302,302)}; %empty variable to create a matrix the size of the whole connectome for each NPP_GPCR network
NPP_individual_networks_nerve_ring_connections_strata_2(:,1) = {zeros(302,302)}; %empty variable to create a matrix the size of the whole connectome for each NPP_GPCR network
NPP_individual_networks_nerve_ring_connections_strata_3(:,1) = {zeros(302,302)}; %empty variable to create a matrix the size of the whole connectome for each NPP_GPCR network
NPP_individual_networks_nerve_ring_connections_strata_4(:,1) = {zeros(302,302)}; %empty variable to create a matrix the size of the whole connectome for each NPP_GPCR network
NPP_individual_networks_nerve_ring_connections_strata_5(:,1) = {zeros(302,302)}; %empty variable to create a matrix the size of the whole connectome for each NPP_GPCR network

for i = 1:length(NR_connections)
    NPP_individual_networks_nerve_ring_connections_strata_1{i,1}(find(NR_strata_class.Strata == 1), find(NR_strata_class.Strata == 1)) = NR_connections{i,1}(find(NR_strata_class.Strata == 1), find(NR_strata_class.Strata == 1)); % network for first strata for each NPP_GPCR pair
    NPP_individual_networks_nerve_ring_connections_strata_2{i,1}(find(NR_strata_class.Strata == 2), find(NR_strata_class.Strata == 2)) = NR_connections{i,1}(find(NR_strata_class.Strata == 2), find(NR_strata_class.Strata == 2)); % network for second strata for each NPP_GPCR pair
    NPP_individual_networks_nerve_ring_connections_strata_3{i,1}(find(NR_strata_class.Strata == 3), find(NR_strata_class.Strata == 3)) = NR_connections{i,1}(find(NR_strata_class.Strata == 3), find(NR_strata_class.Strata == 3)); % network for third strata for each NPP_GPCR pair
    NPP_individual_networks_nerve_ring_connections_strata_4{i,1}(find(NR_strata_class.Strata == 4), find(NR_strata_class.Strata == 4)) = NR_connections{i,1}(find(NR_strata_class.Strata == 4), find(NR_strata_class.Strata == 4)); % network for fourth strata for each NPP_GPCR pair 
    NPP_individual_networks_nerve_ring_connections_strata_5{i,1}(find(NR_strata_class.Strata == 5), find(NR_strata_class.Strata == 5)) = NR_connections{i,1}(find(NR_strata_class.Strata == 5), find(NR_strata_class.Strata == 5)); % network for fifth strata for each NPP_GPCR pair
end

[NPP_in_NR_1_connections_Amatrix, NPP_in_NR_1_connections_Amatrix_rows] = Connections_Amatrix (connections, NPP_individual_networks_nerve_ring_connections_strata_1, neurons, pair_names,...
    'Number of connections in nerve ring strata 1 made by each of the individual NPP-GPCR pairs for each neuron node (91 short range)', [ path 'Comparison_btw_ranges'], date);
[NPP_in_NR_2_connections_Amatrix, NPP_in_NR_2_connections_Amatrix_rows] = Connections_Amatrix (connections, NPP_individual_networks_nerve_ring_connections_strata_2, neurons, pair_names,...
    'Number of connections in nerve ring strata 2 made by each of the individual NPP-GPCR pairs for each neuron node (91 short range)', [ path 'Comparison_btw_ranges'], date);
[NPP_in_NR_3_connections_Amatrix, NPP_in_NR_3_connections_Amatrix_rows] = Connections_Amatrix (connections, NPP_individual_networks_nerve_ring_connections_strata_3, neurons, pair_names,...
    'Number of connections in nerve ring strata 3 made by each of the individual NPP-GPCR pairs for each neuron node (91 short range)', [ path 'Comparison_btw_ranges'], date);
[NPP_in_NR_4_connections_Amatrix, NPP_in_NR_4_connections_Amatrix_rows] = Connections_Amatrix (connections, NPP_individual_networks_nerve_ring_connections_strata_4, neurons, pair_names,...
    'Number of connections in nerve ring strata 4 made by each of the individual NPP-GPCR pairs for each neuron node (91 short range)', [ path 'Comparison_btw_ranges'], date);
[NPP_in_NR_5_connections_Amatrix, NPP_in_NR_5_connections_Amatrix_rows] = Connections_Amatrix (connections, NPP_individual_networks_nerve_ring_connections_strata_5, neurons, pair_names,...
    'Number of connections in nerve ring neurons unassigned to strata and the rest of the nerve ring made by each of the individual NPP-GPCR pairs for each neuron node (91 short range)', [ path 'Comparison_btw_ranges'], date);

%Systematically identify pairs in the NR that only have 1 ligand per
%receptor and prove that neuropeptides can cross strata
%load ('scRNAseq_expression_threshold4_EC50_95/Expression_measures_EC500nM_unweighted.mat', {'GPCR_to_NPP', 'ib', 'count_GPCR'}) %load variables showing number of NPP that bind to each receptor 
%NPP_93_networks_1_ligand = ib(find (accumarray(count_GPCR,1) == 1)); % select rows from NPP_93_networks that have receptors that bind only 1 ligand

[receptors_NPP_string, idx_con, idx_rec] = unique(connections(:, 'GeneGPCR'), 'rows'); % get list of receptors that form a NPP connection and get index to count how many ligands each receptor binds
NPP_individual_networks_1_ligand = sort(idx_con(find(accumarray(idx_rec,1) == 1))); % select rows from NPP_91_networks that have receptors that bind only 1 ligand
receptors_NPP_1_ligand = table2cell(connections(NPP_individual_networks_1_ligand, 'GeneGPCR')); %list of receptors with only 1 ligand
NPP_1_ligand = table2cell(connections(NPP_individual_networks_1_ligand, 'GeneLigand')); %list of NPP that bind to GPCR with only 1 ligand
receptors_NPP_GPCR_pairs_1_ligand = pair_names(NPP_individual_networks_1_ligand, 1); %list of NPP_GPCR pairs with only 1 ligand
[~, idx_con1, idx_rec1] = unique(NPP_1_ligand, 'rows'); % get the neuropeptides that only bind to 1 receptor from the list of NPP that bind receptors that bind only 1 ligand. 
NPP_individual_networks_1_lig_1_rec = NPP_individual_networks_1_ligand(sort(idx_con1(find(accumarray(idx_rec1, 1) == 1)))); % select rows from NPP_91_networks that have receptors that bind only 1 ligand and 1 receptor
NPP_1_lig_1_rec = NPP_1_ligand(sort(idx_con1(find(accumarray(idx_rec1, 1) == 1)))); %get the neuropeptides in which the NPP only has one receptor and the receptor only has 1 ligand
receptors_NPP_GPCR_pairs_1_lig_1_rec = pair_names(NPP_individual_networks_1_lig_1_rec, 1); %list of NPP_GPCR pairs with only 1 ligand
NPP_GPCR_connections_NR_in_btw_strata_Amatrix = NPP_GPCR_connections_NR_Amatrix - (NPP_in_NR_1_connections_Amatrix + NPP_in_NR_2_connections_Amatrix + NPP_in_NR_3_connections_Amatrix + NPP_in_NR_4_connections_Amatrix + NPP_in_NR_5_connections_Amatrix);

for i = 1:length(NPP_individual_networks_1_lig_1_rec) 
    j = NPP_individual_networks_1_lig_1_rec(i);
    [indeg_networks_1_lig_1_rec_short_range(i,:), ~, ~] = degrees_dir(short_individual{j}); %calculate indegree of each of the networks that uses a receptor that only has one ligand
    indeg_networks_1_lig_1_rec_short_range_NR_strata_1(i, :) = NPP_in_NR_1_connections_Amatrix(303:604, j)'; %each row of the indegree matrix is the indegree inside strata 1 of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    indeg_networks_1_lig_1_rec_short_range_NR_strata_2(i, :) = NPP_in_NR_2_connections_Amatrix(303:604, j)'; %each row of the indegree matrix is the indegree inside strata 2 of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    indeg_networks_1_lig_1_rec_short_range_NR_strata_3(i, :) = NPP_in_NR_3_connections_Amatrix(303:604, j)'; %each row of the indegree matrix is the indegree inside strata 3 of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    indeg_networks_1_lig_1_rec_short_range_NR_strata_4(i, :) = NPP_in_NR_4_connections_Amatrix(303:604, j)'; %each row of the indegree matrix is the indegree inside strata 4 of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    indeg_networks_1_lig_1_rec_short_range_NR_strata_5(i, :) = NPP_in_NR_5_connections_Amatrix(303:604, j)'; %each row of the indegree matrix is the indegree in between unassigned neurons and other neurons of the NR of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    indeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata(i, :) = NPP_GPCR_connections_NR_in_btw_strata_Amatrix(303:604, j)'; %each row of the indegree matrix is the indegree in between strata of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    [indeg_networks_1_lig_1_rec_mid_range(i,:), ~, ~] = degrees_dir(mid_individual{j}); %calculate indegree of each of the networks that uses a receptor that only has one ligand
    [indeg_networks_1_lig_1_rec_long_range(i,:), ~, ~] = degrees_dir(long_individual{j}); %calculate indegree of each of the networks that uses a receptor that only has one ligand
end

%do the same for the NPP of those pairs 
for i = 1:length(NPP_individual_networks_1_lig_1_rec) 
    j = NPP_individual_networks_1_lig_1_rec(i);
    [~,outdeg_networks_1_lig_1_rec_short_range(i,:), ~] = degrees_dir(short_individual{j}); %calculate outdegree of each of the networks that uses a receptor that only has one ligand
    outdeg_networks_1_lig_1_rec_short_range_NR_strata_1(i, :) = NPP_in_NR_1_connections_Amatrix(1:302, j)'; %each row of the outdegree matrix is the indegree inside strata 1 of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    outdeg_networks_1_lig_1_rec_short_range_NR_strata_2(i, :) = NPP_in_NR_2_connections_Amatrix(1:302, j)'; %each row of the outdegree matrix is the indegree inside strata 2 of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    outdeg_networks_1_lig_1_rec_short_range_NR_strata_3(i, :) = NPP_in_NR_3_connections_Amatrix(1:302, j)'; %each row of the outdegree matrix is the indegree inside strata 3 of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    outdeg_networks_1_lig_1_rec_short_range_NR_strata_4(i, :) = NPP_in_NR_4_connections_Amatrix(1:302, j)'; %each row of the outdegree matrix is the indegree inside strata 4 of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    outdeg_networks_1_lig_1_rec_short_range_NR_strata_5(i, :) = NPP_in_NR_5_connections_Amatrix(1:302, j)'; %each row of the outdegree matrix is the indegree in between unassigned neurons and other neurons of the NR of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    outdeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata(i, :) = NPP_GPCR_connections_NR_in_btw_strata_Amatrix(1:302, j)'; %each row of the outdegree matrix is the indegree in between strata of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    [~,outdeg_networks_1_lig_1_rec_mid_range(i,:), ~] = degrees_dir(mid_individual{j}); %calculate outdegree of each of the networks that uses a receptor that only has one ligand
    [~,outdeg_networks_1_lig_1_rec_long_range(i,:), ~] = degrees_dir(long_individual{j}); %calculate outdegree of each of the networks that uses a receptor that only has one ligand
end


%Plot with 4 colour, 1 for each range of connections and a
%4th for the inside strata connections

image12 = figure('Visible', 'off', 'Resize', 'off', 'Position', [199.625,302,3600,704]); %create an image to save the figure with the right size and prevent it being resize by the plots
[y_indeg_networks_1_lig_1_rec_long_range_sort, x_indeg_networks_1_lig_1_rec_long_range_sort ] = find(indeg_networks_1_lig_1_rec_long_range); %get x and y coordinates to plot scatter
s1 = scatter(x_indeg_networks_1_lig_1_rec_long_range_sort, y_indeg_networks_1_lig_1_rec_long_range_sort, [], 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'none'); %plot scatter of long range (scatter instead of spy to choose colour)
hold on 
[y_indeg_networks_1_lig_1_rec_mid_range_sort, x_indeg_networks_1_lig_1_rec_mid_range_sort ] = find(indeg_networks_1_lig_1_rec_mid_range); %get x and y coordinates to plot scatter
s2 = scatter(x_indeg_networks_1_lig_1_rec_mid_range_sort, y_indeg_networks_1_lig_1_rec_mid_range_sort, [], 'MarkerFaceColor', 1/255*[215 25 28], 'MarkerEdgeColor', 'none'); %plot scatter of mid range (scatter instead of spy to choose colour)
hold on 
[y_indeg_networks_1_lig_1_rec_short_range_sort, x_indeg_networks_1_lig_1_rec_short_range_sort ] = find(indeg_networks_1_lig_1_rec_short_range); %get x and y coordinates to plot scatter
s3 = scatter(x_indeg_networks_1_lig_1_rec_short_range_sort, y_indeg_networks_1_lig_1_rec_short_range_sort, [], 'MarkerFaceColor', 1/255*[171 217 233], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
hold on 
[y_indeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata_sort, x_indeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata_sort ] = find(indeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata); %get x and y coordinates to plot scatter
s9 = scatter(x_indeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata_sort, y_indeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata_sort, [], 'MarkerFaceColor', 1/255*[253,141,60], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
hold on 
[y_indeg_networks_1_lig_1_rec_short_range_NR_strata_1_sort, x_indeg_networks_1_lig_1_rec_short_range_NR_strata_1_sort ] = find(indeg_networks_1_lig_1_rec_short_range_NR_strata_1); %get x and y coordinates to plot scatter
s4 = scatter(x_indeg_networks_1_lig_1_rec_short_range_NR_strata_1_sort, y_indeg_networks_1_lig_1_rec_short_range_NR_strata_1_sort, [], 'MarkerFaceColor', 1/255*[44 123 182], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
hold on
[y_indeg_networks_1_lig_1_rec_short_range_NR_strata_2_sort, x_indeg_networks_1_lig_1_rec_short_range_NR_strata_2_sort ] = find(indeg_networks_1_lig_1_rec_short_range_NR_strata_2); %get x and y coordinates to plot scatter
s5 = scatter(x_indeg_networks_1_lig_1_rec_short_range_NR_strata_2_sort, y_indeg_networks_1_lig_1_rec_short_range_NR_strata_2_sort, [], 'MarkerFaceColor', 1/255*[44 123 182], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
hold on
[y_indeg_networks_1_lig_1_rec_short_range_NR_strata_3_sort, x_indeg_networks_1_lig_1_rec_short_range_NR_strata_3_sort ] = find(indeg_networks_1_lig_1_rec_short_range_NR_strata_3); %get x and y coordinates to plot scatter
s6 = scatter(x_indeg_networks_1_lig_1_rec_short_range_NR_strata_3_sort, y_indeg_networks_1_lig_1_rec_short_range_NR_strata_3_sort, [], 'MarkerFaceColor', 1/255*[44 123 182], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
hold on
[y_indeg_networks_1_lig_1_rec_short_range_NR_strata_4_sort, x_indeg_networks_1_lig_1_rec_short_range_NR_strata_4_sort ] = find(indeg_networks_1_lig_1_rec_short_range_NR_strata_4); %get x and y coordinates to plot scatter
s7 = scatter(x_indeg_networks_1_lig_1_rec_short_range_NR_strata_4_sort, y_indeg_networks_1_lig_1_rec_short_range_NR_strata_4_sort, [], 'MarkerFaceColor', 1/255*[44 123 182], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
hold on
[y_indeg_networks_1_lig_1_rec_short_range_NR_strata_5_sort, x_indeg_networks_1_lig_1_rec_short_range_NR_strata_5_sort ] = find(indeg_networks_1_lig_1_rec_short_range_NR_strata_5); %get x and y coordinates to plot scatter
s8 = scatter(x_indeg_networks_1_lig_1_rec_short_range_NR_strata_5_sort, y_indeg_networks_1_lig_1_rec_short_range_NR_strata_5_sort, [], 'MarkerFaceColor', 1/255*[44 123 182], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
hold off
title('Type of connections of NPP-GPCR pairs where the GPCR only binds 1 ligand, ligand binds 1 GPCR(indegree)', 'interpreter', 'none'); %label the title of the plot
% yticks(linspace(1, size(indeg_networks_1_ligand_short_range, 1), length(receptors_NPP_1_ligand))); %add all x axis ticks
% yticklabels(receptors_NPP_1_ligand);  %labels for y axes
% ylabel('GPCR') %y axes label
% xticks(linspace(1, 302, 302)); %add all x axis ticks
% xticklabels(neuron_ID.nodeLabel);  %labels for y axes
% xtickangle(90); %change orientation of labels
% xlabel('Neurons') %x axes label
set(gca, 'TickLength', [0, 0], 'TickLabelInterpreter','none', 'Fontsize', 7, 'YColor', 'k', 'XColor', 'k'); %remove tick marks and prevent unwanted underscores in tick labels
box on 
legend([s2 s9 s3 s4], 'Requires mid range', 'Requires in between strata', 'Whithin bundle (outside NR)', 'Whithin stratum', 'FontSize', 40, 'Location', 'bestoutside', 'Orientation', 'horizontal')
saveas(image12, ['Results/scRNAseq_C.elegans/' path 'Comparison_btw_ranges/' date 'All neurons connections for each stratum and each range for each NPP_GPCR pair with 1 ligand 1 GPCR in the neuropeptidergic network (91 short range network indegree new 5 colours).png']);%save image as png

%do the same for the NPP of those pairs 

image13 = figure('Visible', 'off', 'Resize', 'off', 'Position', [199.625,350,3600,704]); %create an image to save the figure with the right size and prevent it being resize by the plots
[y_outdeg_networks_1_lig_1_rec_long_range_sort, x_outdeg_networks_1_lig_1_rec_long_range_sort ] = find(outdeg_networks_1_lig_1_rec_long_range); %get x and y coordinates to plot scatter
s1 = scatter(x_outdeg_networks_1_lig_1_rec_long_range_sort, y_outdeg_networks_1_lig_1_rec_long_range_sort, [], 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'none'); %plot scatter of long range (scatter instead of spy to choose colour)
hold on 
[y_outdeg_networks_1_lig_1_rec_mid_range_sort, x_outdeg_networks_1_lig_1_rec_mid_range_sort ] = find(outdeg_networks_1_lig_1_rec_mid_range); %get x and y coordinates to plot scatter
s2 = scatter(x_outdeg_networks_1_lig_1_rec_mid_range_sort, y_outdeg_networks_1_lig_1_rec_mid_range_sort, [], 'MarkerFaceColor', 1/255*[215 25 28], 'MarkerEdgeColor', 'none'); %plot scatter of mid range (scatter instead of spy to choose colour)
hold on 
[y_outdeg_networks_1_lig_1_rec_short_range_sort, x_outdeg_networks_1_lig_1_rec_short_range_sort ] = find(outdeg_networks_1_lig_1_rec_short_range); %get x and y coordinates to plot scatter
s3 = scatter(x_outdeg_networks_1_lig_1_rec_short_range_sort, y_outdeg_networks_1_lig_1_rec_short_range_sort, [], 'MarkerFaceColor', 1/255*[171 217 233], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
hold on 
[y_outdeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata_sort, x_outdeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata_sort ] = find(outdeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata); %get x and y coordinates to plot scatter
s9 = scatter(x_outdeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata_sort, y_outdeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata_sort, [], 'MarkerFaceColor', 1/255*[253,141,60], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
hold on
[y_outdeg_networks_1_lig_1_rec_short_range_NR_strata_1_sort, x_outdeg_networks_1_lig_1_rec_short_range_NR_strata_1_sort ] = find(outdeg_networks_1_lig_1_rec_short_range_NR_strata_1); %get x and y coordinates to plot scatter
s4 = scatter(x_outdeg_networks_1_lig_1_rec_short_range_NR_strata_1_sort, y_outdeg_networks_1_lig_1_rec_short_range_NR_strata_1_sort, [], 'MarkerFaceColor', 1/255*[44 123 182], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
hold on
[y_outdeg_networks_1_lig_1_rec_short_range_NR_strata_2_sort, x_outdeg_networks_1_lig_1_rec_short_range_NR_strata_2_sort ] = find(outdeg_networks_1_lig_1_rec_short_range_NR_strata_2); %get x and y coordinates to plot scatter
s5 = scatter(x_outdeg_networks_1_lig_1_rec_short_range_NR_strata_2_sort, y_outdeg_networks_1_lig_1_rec_short_range_NR_strata_2_sort, [], 'MarkerFaceColor', 1/255*[44 123 182], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
hold on
[y_outdeg_networks_1_lig_1_rec_short_range_NR_strata_3_sort, x_outdeg_networks_1_lig_1_rec_short_range_NR_strata_3_sort ] = find(outdeg_networks_1_lig_1_rec_short_range_NR_strata_3); %get x and y coordinates to plot scatter
s6 = scatter(x_outdeg_networks_1_lig_1_rec_short_range_NR_strata_3_sort, y_outdeg_networks_1_lig_1_rec_short_range_NR_strata_3_sort, [], 'MarkerFaceColor', 1/255*[44 123 182], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
hold on
[y_outdeg_networks_1_lig_1_rec_short_range_NR_strata_4_sort, x_outdeg_networks_1_lig_1_rec_short_range_NR_strata_4_sort ] = find(outdeg_networks_1_lig_1_rec_short_range_NR_strata_4); %get x and y coordinates to plot scatter
s7 = scatter(x_outdeg_networks_1_lig_1_rec_short_range_NR_strata_4_sort, y_outdeg_networks_1_lig_1_rec_short_range_NR_strata_4_sort, [], 'MarkerFaceColor', 1/255*[44 123 182], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
hold on
[y_outdeg_networks_1_lig_1_rec_short_range_NR_strata_5_sort, x_outdeg_networks_1_lig_1_rec_short_range_NR_strata_5_sort ] = find(outdeg_networks_1_lig_1_rec_short_range_NR_strata_5); %get x and y coordinates to plot scatter
s8 = scatter(x_outdeg_networks_1_lig_1_rec_short_range_NR_strata_5_sort, y_outdeg_networks_1_lig_1_rec_short_range_NR_strata_5_sort, [], 'MarkerFaceColor', 1/255*[44 123 182], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
hold off
title('Type of connections of NPP-GPCR pairs where the GPCR only binds 1 ligand, ligand binds 1 GPCR(outdegree)', 'interpreter', 'none'); %label the title of the plot
% yticks(linspace(1, size(indeg_networks_1_ligand_short_range, 1), length(receptors_NPP_1_ligand))); %add all x axis ticks
% yticklabels(NPP_1_ligand);  %labels for y axes
% ylabel('NPP') %y axes label
% xticks(linspace(1, 302, 302)); %add all x axis ticks
% xticklabels(neuron_ID.nodeLabel);  %labels for y axes
% xtickangle (90); %change orientation of labels
% xlabel('Neurons') %x axes label
box on 
set(gca, 'TickLength', [0, 0], 'TickLabelInterpreter','none', 'Fontsize', 7, 'YColor', 'k', 'XColor', 'k'); %remove tick marks and prevent unwanted underscores in tick labels
legend([s2 s9 s3 s4], 'Requires mid range', 'Requires in between strata', 'Whithin bundle (outside NR)', 'Whithin stratum', 'FontSize', 40, 'Location', 'bestoutside', 'Orientation', 'horizontal')
saveas(image13, ['Results/scRNAseq_C.elegans/' path 'Comparison_btw_ranges/' date 'All neurons connections for each stratum and each range for each NPP_GPCR pair with 1 ligand 1 GPCR in the neuropeptidergic network (91 short range network outdegree new 5 colours).png']);%save image as png



%****************************************************************************************************************
%Scatter plot colour incoming and outgoing connections based on anatomy
%****************************************************************************************************************

%Use the scatter plot to plot all neurons of the 93 networks but colour by
%the anatomical mid range area they belong to (head, tail, midbody)
for i = 1:length(short_individual) %loop over every neuron 
    [indegree_NPP_individual_short_range(i,:), ~,~] = degrees_dir(short_individual{i}); %indegree of short range connections by NPP_GPCR pairs
    [ ~, outdegree_NPP_individual_short_range(i,:), ~] = degrees_dir(short_individual{i}); %outdegree of short range connections by NPP_GPCR pairs
end
[y_indeg_NPP_individual_short_range, x_indeg_NPP_individual_short_range ] = find(indegree_NPP_individual_short_range); %get x and y coordinates to plot scatter
[y_outdeg_NPP_individual_short_range, x_outdeg_NPP_individual_short_range ] = find(outdegree_NPP_individual_short_range); %get x and y coordinates to plot scatter

colormap = lines; %choose colormap

for i = 1: length(x_indeg_NPP_individual_short_range)
    if ismember(x_indeg_NPP_individual_short_range(i), class_Head) 
       colors_neurons_Head_ind(i,:) = colormap(4, :);
    else
       colors_neurons_Head_ind(i,:) = [0.5 0.5 0.5];
    end
    if ismember(x_indeg_NPP_individual_short_range(i), class_Midbody)
       colors_neurons_Midbody_ind(i,:) = colormap(2, :);
    else
       colors_neurons_Midbody_ind(i,:) = [0.5 0.5 0.5];
    end
    if ismember(x_indeg_NPP_individual_short_range(i), class_Tail)
       colors_neurons_Tail_ind(i,:) = colormap(3, :);
    else 
       colors_neurons_Tail_ind(i,:) = [0.5 0.5 0.5];
    end
end
for i = 1: length(x_outdeg_NPP_individual_short_range)
    if ismember(x_outdeg_NPP_individual_short_range(i), class_Head) 
       colors_neurons_Head_out(i,:) = colormap(4, :);
    else
       colors_neurons_Head_out(i,:) = [0.5 0.5 0.5];
    end
    if ismember(x_outdeg_NPP_individual_short_range(i), class_Midbody)
       colors_neurons_Midbody_out(i,:) = colormap(2, :);
    else
       colors_neurons_Midbody_out(i,:) = [0.5 0.5 0.5];
    end
    if ismember(x_outdeg_NPP_individual_short_range(i), class_Tail)
       colors_neurons_Tail_out(i,:) = colormap(3, :);
    else 
       colors_neurons_Tail_out(i,:) = [0.5 0.5 0.5];
    end
end
image15 = figure('Visible', 'off', 'Resize', 'off', 'Position', [199.625,350,3600,2204]); %create an image to save the figure with the right size and prevent it being resize by the plots
t = tiledlayout(3,1);
nexttile 
s1 = scatter(x_indeg_NPP_individual_short_range, y_indeg_NPP_individual_short_range, [], colors_neurons_Head_ind, 'filled'); %plot scatter of short range (scatter instead of spy to choose colour)
yticks( linspace(1, length(pair_names), length(pair_names))); %add all y axis ticks
yticklabels(pair_names); %labels for y axes
set(gca, 'TickLength', [0, 0], 'TickLabelInterpreter','none', 'Fontsize', 7); %remove tick marks and prevent unwanted underscores in tick labels
nexttile
s2 = scatter(x_indeg_NPP_individual_short_range, y_indeg_NPP_individual_short_range, [], colors_neurons_Midbody_ind, 'filled'); %plot scatter of short range (scatter instead of spy to choose colour)
yticks( linspace(1, length(pair_names), length(pair_names))); %add all y axis ticks
yticklabels(pair_names); %labels for y axes
set(gca, 'TickLength', [0, 0], 'TickLabelInterpreter','none', 'Fontsize', 7); %remove tick marks and prevent unwanted underscores in tick labels
nexttile
s3 = scatter(x_indeg_NPP_individual_short_range, y_indeg_NPP_individual_short_range, [], colors_neurons_Tail_ind, 'filled'); %plot scatter of short range (scatter instead of spy to choose colour)
yticks( linspace(1, length(pair_names), length(pair_names))); %add all y axis ticks
yticklabels(pair_names); %labels for y axes
ylabel(t, 'NPP-GPCR expression') %y axes label
% xticks(linspace(1, 302, 302)); %add all x axis ticks
% xticklabels(neurons.nodeLabel); %labels for x axes
% xlabel(t, 'Neurons') %x axes label
% xtickangle (90); %change orientation of labels
set(gca, 'TickLength', [0, 0], 'TickLabelInterpreter','none', 'Fontsize', 7); %remove tick marks and prevent unwanted underscores in tick labels
legend([s1 s2 s3], 'Process bundles in the Head', 'Process bundles in Midbody', 'Process bundles in Tail', 'FontSize', 12, 'Location', 'northoutside', 'Orientation', 'horizontal')
saveas(image15, ['Results/scRNAseq_C.elegans/' path 'Comparison_btw_ranges/' date 'All neurons indegree for all NPP-GPCR pairs coloured by anatomy short range.png']);%save image as png

image16 = figure('Visible', 'off', 'Resize', 'off', 'Position', [199.625,350,3600,2204]); %create an image to save the figure with the right size and prevent it being resize by the plots
t = tiledlayout(3,1);
nexttile 
s1 = scatter(x_outdeg_NPP_individual_short_range, y_outdeg_NPP_individual_short_range, [], colors_neurons_Head_out, 'filled'); %plot scatter of short range (scatter instead of spy to choose colour)
yticks( linspace(1, length(pair_names), length(pair_names))); %add all y axis ticks
yticklabels(pair_names); %labels for y axes
set(gca, 'TickLength', [0, 0], 'TickLabelInterpreter','none', 'Fontsize', 7); %remove tick marks and prevent unwanted underscores in tick labels
nexttile
s2 = scatter(x_outdeg_NPP_individual_short_range, y_outdeg_NPP_individual_short_range, [], colors_neurons_Midbody_out, 'filled'); %plot scatter of short range (scatter instead of spy to choose colour)
yticks( linspace(1, length(pair_names), length(pair_names))); %add all y axis ticks
yticklabels(pair_names); %labels for y axes
set(gca, 'TickLength', [0, 0], 'TickLabelInterpreter','none', 'Fontsize', 7); %remove tick marks and prevent unwanted underscores in tick labels
nexttile
s3 = scatter(x_outdeg_NPP_individual_short_range, y_outdeg_NPP_individual_short_range, [], colors_neurons_Tail_out, 'filled'); %plot scatter of short range (scatter instead of spy to choose colour)
yticks( linspace(1, length(pair_names), length(pair_names))); %add all y axis ticks
yticklabels(pair_names); %labels for y axes
ylabel(t, 'NPP-GPCR expression') %y axes label
% xticks(linspace(1, 302, 302)); %add all x axis ticks
% xticklabels(neurons.nodeLabel); %labels for x axes
% xlabel(t, 'Neurons') %x axes label
% xtickangle (90); %change orientation of labels
set(gca, 'TickLength', [0, 0], 'TickLabelInterpreter','none', 'Fontsize', 7); %remove tick marks and prevent unwanted underscores in tick labels
legend([s1 s2 s3], 'Process bundles in the Head', 'Process bundles in Midbody', 'Process bundles in Tail', 'FontSize', 12, 'Location', 'northoutside', 'Orientation', 'horizontal')
saveas(image16, ['Results/scRNAseq_C.elegans/' path 'Comparison_btw_ranges/' date 'All neurons outdegree for all NPP-GPCR pairs coloured by anatomy short range.png']);%save image as png

%save results
save(['Results/scRNAseq_C.elegans/' path2 date 'Neuropeptide_Networks_Comparison_EC500nM.mat']) %save all variables from workspace to results at this point

end

function [Amatrix,row_names] = Connections_Amatrix (data_connections, data_networks, neuron_names, connection_names, name_plot, folder, date)

%calculate matrix 
vector = zeros(2*length(data_networks{1}), 1); %empty vector to store the data of each iteration
if iscell (data_networks)
    for i = 1:height(data_connections)
        for j = 1:length(data_networks{1})
            vector(j,i) = sum(data_networks{i}(j, :));%add all the edges in a row for the matrix of the NPP_GPCR pair
            vector(j+length(data_networks{1}), i) = sum(data_networks{i}(:,j)); %add all the edges for the same neuron in the other dimension of the matrix of the NPP_GPCR pair
        end
    end
    Amatrix = vector; %end result
    row_names = cat(1, neuron_names.nodeLabel, neuron_names.nodeLabel); %names of matrix rows
else %the input is a matrix not a cell array, output is a column vector
    vector = zeros(2*length(data_networks), 1); %store the data temporary in iteration
        for j = 1:length(data_networks)
            vector(j,1) = sum(data_networks(j, :));%add all the edges in a row for the matrix of the NPP_GPCR pair
            vector(j+length(data_networks), 1) = sum(data_networks(:,j)); %add all the edges for the same neuron in the other dimension of the matrix of the NPP_GPCR pair
        end
     Amatrix = vector; %end result
     row_names = cat(1, neuron_names.nodeLabel, neuron_names.nodeLabel); %names of matrix rows
end 

%plot matrix 
image = figure('Visible', 'off', 'Position', [183,114,1575,904]);
cspy(Amatrix, 'ColorMap', 'jet'); %plots adjacency matrix as a sparse plot
title(name_plot, 'interpreter', 'none'); %adds title
xlabel ('Neuropeptide-receptor pair');
xticks(linspace(1, size(Amatrix, 2), length(connection_names))); %add all x axis ticks
xticklabels(connection_names); %add x axis tick labels with the pair names
xtickangle(90); %rotate axis labels so that they fit
ylabel ('Neuron');
set(gca,'TickLabelInterpreter','none')
drawnow; %draws figure
saveas(image, ['Results/scRNAseq_C.elegans/' folder '/' date name_plot '.png']);%save image as png

end
%**************************************************************************
function Individual_networks_analysis(network1, network2, threshold, path, path2, date)

%% ****************************************************************************************************************
%SET DEFAULTS 
% ****************************************************************************************************************

set(0, 'DefaultAxesFontName', 'Arial'); %set default text

colormap = lines; 

%% ****************************************************************************************************************
%Pervasive plot (assortativity analysis)
%****************************************************************************************************************

%Calculate the number of neurons that express each neuropeptide and each
%receptor in each of the 91 individual neuropeptide receptor networks
for i = 1:length(NPP_individual_networks_short_range_connections)
    NPP_individual_networks_neurons_express_NPP_shortr(i) = nnz(any(network1{i},2)); %number of non zero rows in matrix
    NPP_individual_networks_neurons_express_GPCR_shortr(i) = nnz(any(network1{i},1));%number of non zero columns in matrix
    NPP_individual_networks_neurons_express_NPP_midr(i) = nnz(any(network2{i},2)); %number of non zero rows in matrix
    NPP_individual_networks_neurons_express_GPCR_midr(i) = nnz(any(network2{i},1));%number of non zero columns in matrix
end

%Create a colormap depending if network is broadcaster, integrative,
%pervasive or local depending on the number of neurons that express each
%neuropeptide and each receptor

for i = 1:length(network1)
    if NPP_individual_networks_neurons_express_NPP_shortr(i) <= threshold && NPP_individual_networks_neurons_express_GPCR_shortr(i) <= threshold %local networks
        NPP_individual_networks_short_range_connections_pervasive_noloc{i} = []; %cell array empty for local networks 
        assortativity_type_network_sr{i,1} = 1;
        assortativity_type_network_sr{i,2} = 'local';
        colors_pervasive_shortr{i} = colormap(3,:); %color yellow
        colors_pervasive_shortr_nolocal{i} = []; %empty
    elseif NPP_individual_networks_neurons_express_NPP_shortr(i) >= threshold && NPP_individual_networks_neurons_express_GPCR_shortr(i) >= threshold %pervasive networks
        assortativity_type_network_sr{i,1} = 2;
        assortativity_type_network_sr{i,2} = 'pervasive';
        NPP_individual_networks_short_range_connections_pervasive_noloc{i} = NPP_individual_networks_short_range_connections{i}; % cell array add pervasive networks
        colors_pervasive_shortr{i} = colormap(5,:); %color green
        colors_pervasive_shortr_nolocal{i} = colormap(5,:); %color green
    elseif NPP_individual_networks_neurons_express_NPP_shortr(i) > threshold && NPP_individual_networks_neurons_express_GPCR_shortr(i) < threshold %integrative
        assortativity_type_network_sr{i,1} = 3;
        assortativity_type_network_sr{i,2} = 'integrative';
        NPP_individual_networks_short_range_connections_pervasive_noloc{i} = NPP_individual_networks_short_range_connections{i}; % cell array add integrative networks
        colors_pervasive_shortr{i} = colormap (1, :); %color blue
        colors_pervasive_shortr_nolocal{i} = colormap (1, :); %color blue
    else NPP_individual_networks_neurons_express_NPP_shortr(i) < threshold && NPP_individual_networks_neurons_express_GPCR_shortr(i) > threshold; %broadcaster
        assortativity_type_network_sr{i,1} = 4;
        assortativity_type_network_sr{i,2} = 'broadcaster';
        NPP_individual_networks_short_range_connections_pervasive_noloc{i} = NPP_individual_networks_short_range_connections{i}; % cell array add broadcaster networks
        colors_pervasive_shortr{i} = colormap(7,:); %color red
        colors_pervasive_shortr_nolocal{i} = colormap(7,:); %color red
    end
    if NPP_individual_networks_neurons_express_NPP_midr(i) <= threshold && NPP_individual_networks_neurons_express_GPCR_midr(i) <= threshold %local networks
        assortativity_type_network_mr{i,1} = 1;
        assortativity_type_network_mr{i,2} = 'local';
        colors_pervasive_midr{i} = colormap(3,:); %color yellow
        colors_pervasive_midr_nolocal{i} = [0.5 0.5 0.5]; %grey to be removed afterwards
    elseif NPP_individual_networks_neurons_express_NPP_midr(i) >= threshold && NPP_individual_networks_neurons_express_GPCR_midr(i) >= threshold %pervasive networks
        assortativity_type_network_mr{i,1} = 2;
        assortativity_type_network_mr{i,2} = 'pervasive';
        colors_pervasive_midr{i} = colormap(5,:); %color green
        colors_pervasive_midr_nolocal{i} = [0.5 0.5 0.5]; %grey
    elseif NPP_individual_networks_neurons_express_NPP_midr(i) > threshold && NPP_individual_networks_neurons_express_GPCR_midr(i) < threshold %integrative
        assortativity_type_network_mr{i,1} = 3;
        assortativity_type_network_mr{i,2} = 'integrative';
        colors_pervasive_midr{i} = colormap (1, :); %color blue
        colors_pervasive_midr_nolocal{i} = [0.5 0.5 0.5]; %grey
    else NPP_individual_networks_neurons_express_NPP_midr(i) < threshold && NPP_individual_networks_neurons_express_GPCR_midr(i) > threshold; %broadcaster
        assortativity_type_network_mr{i,1} = 4;
        assortativity_type_network_mr{i,2} = 'broadcaster';
        colors_pervasive_midr{i} = colormap(7,:); %color red
        colors_pervasive_midr_nolocal{i} = [0.5 0.5 0.5]; %grey
    end
end

%create sets of NPP individual networks without local networks to plot
NPP_individual_networks_short_range_connections_pervasive_noloc = network1(~cellfun(@isempty, NPP_individual_networks_short_range_connections_pervasive_noloc)); %remove local networks, that are empty based on the previous selection we did 
NPP_individual_networks_mid_range_connections_pervasive_noloc = network2(~cellfun(@isempty, NPP_individual_networks_short_range_connections_pervasive_noloc)); %remove local networks, that are empty based on the previous selection we did in short range 
colors_pervasive_shortr_nolocal = colors_pervasive_shortr_nolocal(~cellfun(@isempty, colors_pervasive_shortr_nolocal)); %remove local networks empty cells
colors_pervasive_midr_nolocal = colors_pervasive_midr_nolocal(~cellfun(@isempty, colors_pervasive_shortr_nolocal)); %remove local networks empty cells

%plot pervasive analysis scatter plot
image11 = figure('Visible', 'on', 'Resize', 'off', 'Position', [183, 114, 600, 600]); %create an image to save the figure with the right size and prevent it being resize by the plots
scatter(NPP_individual_networks_neurons_express_GPCR_shortr, NPP_individual_networks_neurons_express_NPP_shortr', 'k', 'o', 'filled') %scatter plot, represent points as black filled dots
xlabel('Neurons that express GPCR') %y axes label
ylabel('Neurons that express NPP') %y axes label
xlim([0 250]) %set axes limit to make square plot
ylim([0 250]) %set axes limit to make square plot
set(gca, 'Fontsize', 20)
saveas(image11, ['Results/scRNAseq_C.elegans' path 'Images_neuropeptide_networks_unweighted/Individual_NPP_networks/' date 'NPP_short_range_broadcaster_analysis_plot(colours and labels).pdf']);%save image as png

%Plot the 91 networks coloured by pervasive type for short and mid range 

plot_digraphs (NPP_individual_networks_short_range_connections, 10, 10, colors_pervasive_shortr, path, date, 'NPP_short_range_individual_networks_pervasive_colors')
plot_digraphs (NPP_individual_networks_mid_range_connections, 10, 10, colors_pervasive_midr, path, date, 'NPP_mid_range_individual_networks_pervasive_colors')


%Plot the broadcaster, pervasive and integrative networks for short and mid
%range and overlap making mid range grey 
plot_digraphs (NPP_individual_networks_short_range_connections_pervasive_noloc, 8, 5, colors_pervasive_shortr_nolocal, path, date, 'NPP_short_range_individual_networks_pervasive_colors_nolocal')
plot_digraphs (NPP_individual_networks_mid_range_connections_pervasive_noloc, 8, 5, colors_pervasive_midr_nolocal, path, date, 'NPP_mid_range_individual_networks_pervasive_grey_nolocal')

save((['Results/scRNAseq_C.elegans' path2 date 'Networks_Neuropeptide_individual_analysis_EC500nM.mat']))

end

function plot_digraphs (networks, subplotx, subploty, colors, path, date, name)
    image12 = figure('Visible', 'off', 'Resize', 'off', 'Position', [183,114,1575,1004]); %create an image to save the figure with the right size and prevent it being resize by the plots
    %image12 = figure('Visible', 'off', 'Resize', 'off', 'Position',
    %[183,114,788,1004]); %create an image to save the figure with the
    %right size and prevent it being resize by the plots (size for non
    %local plots
    for i = 1:length(networks) %iterate over the number of networks we have
    subplot(subplotx,subploty,i) %make 100 squares subplot 
    graph = digraph(networks{i}, 'omitselfloops'); %make the directed graph for each network without self loops
    [~,~,degree]= degrees_dir (networks{i}); %get the degree to remove nodes that have degree 0
    zero_nodes = find( degree == 0); %remove nodes that have no edges
    plot(rmnode(graph, zero_nodes), 'Layout', 'force', 'WeightEffect', 'direct', 'UseGravity', true, 'ShowArrows', 'off', 'EdgeColor', 'k', 'NodeColor', colors{i}, 'NodeLabel', {}) %plot using force plot
    end
    saveas(image12, [path 'Images_neuropeptide_networks_unweighted/Individual_NPP_networks/' date name '.png']);%save image as png
end
%**************************************************************************

