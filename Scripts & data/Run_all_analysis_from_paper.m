%Run all analysis for the paper

%Created by LIDIA RIPOLL SANCHEZ, for L.Ripoll-Sanchez et al. 2023, Neuron
% project 21/11/2023
%Builds the neuropeptide network and gets degree and rich club measures for
%it, it also performs the dimensionality reduction analysis to obtain the
%mesoscale structure, compares the different diffusion range models,
%calculates network measures for the individual networks and performs a
%co-expression analysis

version = 'R2023a'; %MATLAB version when this workspace was saved

%% ************************************************************************
%TO SAVE DATA
path = '/scRNAseq_expression_threshold4_EC50_FINAL/Networks/'; %create folder to save images
path2 = '/scRNAseq_expression_threshold4_EC50_FINAL/Workspaces/'; %create folder to save images
addpath([ 'Results/scRNAseq_C.elegans' path]) %add the folder to the main path so it is easier to search 
addpath([ 'Results/scRNAseq_C.elegans' path2]) %add the folder to the main path so it is easier to search 
date = '08062023_'; %to save workspace with the date

%IMPORT DATA

expression_byneuron = readtable ('30072020_CENGEN_threshold4_expression_NPP_NPR_MR_LGC_allneurons'); %expression of NPP, GPCR and Mon receptors
connections_NPP = readtable ('20230602_Ligand_receptors_interactions_all_EC50_explicit_500nM.csv'); %connections between NPP and GPCR
connections_Mon = readtable ('05062023_Ligand-receptor-interactions_Mon_from_Barry_Bentley'); %connections between Mon and receptors
neuron_ID = readtable('26012022_num_neuronID.txt');
anatomical_class = readtable('072022_anatomical_class.csv'); %import list of neuron cell body and process location (in binary)

%RUN ANALYSIS 
Network_building(connections_NPP, connections_Mon, expression_byneuron, neuron_ID.nodeLabel, anatomical_class, path, path2, date, version)
Network_measures(neuron_ID.nodeLabel, path, path2, date, version)
Network_degree(neuron_ID.nodeLabel, path, path2, date, version)
Network_dimensionality_reduction(neuron_ID.nodeLabel, anatomical_class.FinalClassification, path, path2, date, version)
Network_models_comparison(connections_NPP, neuron_ID.nodeLabel, path, path2, date, version)
Individual_networks_analysis(NPP_individual_networks_short_range_connections, NPP_individual_networks_mid_range_connections, 50, path, path2, date)
Co_expression_analysis(connections_NPP, path, path2, date, version)

%% ************************************************************************
%LOCAL FUNCTIONS 
%**************************************************************************
function Network_building(connections1, connections2, expressionneuron, neuron_names, anatomy, path, path2, date, version)
%This function builds the neuropeptide receptor networks, then plots them,
%masks them to make them unweighted and finally groups several networks
%depending on certain conditions given. 
%% ****************************************************************************************************************
%ADAPT AND LOAD DATA
%****************************************************************************************************************
%Load data
geneID = readtable('20230601_Gene_ID_GPCR_common_name.xlsx'); %NPP and GPCR gene common names in same order as the genes in the expression dataset
expression_byneuron_neurotransmitter = readtable('20230605_Neurotransmitters_Expression_(L.Pereira et al. 2015 Hobert lab).csv', 'ReadRownames', true); %import list of neurotrasmitter per neuron with rownames to sort neurons, includes expression of Mon
synA = readmatrix('Syn_network.csv'); %We have two datasets of synaptic data (Albertson(AW) and WormWiring (WW) this is the AW with the new classification)
gapA = readmatrix('GJ_network.csv'); %Import the gap junctions table

%Adapt expression data
[expression_byneuron_ordered] = Adapt_expression(expressionneuron, neuron_names, geneID); %Adapt the expression values split by neuron type to neuron (from 180 columns to 302 columns)
[expression_byneuron_ordered_unweighted] = Mask_network(expression_byneuron_ordered); %Unweight expression values, when expression value is 1 otherwise is 0

%Adapt neurotransmitter data
expression_byneuron_neurotransmitter = rows2vars(expression_byneuron_neurotransmitter(neuron_names, :)); %sort data based on neuron_ID
expression_byneuron_neurotransmitter.Properties.VariableNames([1]) = {'Gene'}; %change the name of the first column to Gene, since after transposing it now has the default name

%Create full expression table
expression_byneuron_ordered_unweighted_all = [expression_byneuron_ordered_unweighted;expression_byneuron_neurotransmitter]; %add the CeNGEN expression unweighted dataset to the monoamine expression dataset

%% ****************************************************************************************************************
%MONOAMINES
%****************************************************************************************************************

%Build individual networks
[monoamine_networks_all,structure_Mon,pair_names_Mon] = Neuropeptide_Networks(connections2, expression_byneuron_ordered_unweighted_all);

%Group networks
monoamines = matlab.lang.makeValidName(unique(connections2.GeneLigand)); %create list of the neuropeptides uses in the edges to group based on them 
[~, Mon_all_grouped_networks] = Group_networks (pair_names_Mon, monoamine_networks_all, '\w*_\w*'); %group all networks
[Mon_grouped_names, Mon_grouped_networks] = Group_networks (pair_names_Mon, monoamine_networks_all, monoamines); %group based on individual monoamines  

%Plot networks 
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Images_monoamine_networks/Individual_Monoamine_networks']) %create folder to save images
names_Mon = regexprep(pair_names_Mon, {'_(?!\D)', '_(?!\d)'},{'-', ' '});
Plot_Networks(monoamine_networks_all, names_Mon, ([path 'Images_monoamine_networks/Individual_Monoamine_networks']), 'Mon_networks'); %plot and save all adjacency matrices from monoamine_networks
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Images_monoamine_networks/All_monoamine_grouped']) %create folder to save images
Plot_Networks(Mon_all_grouped_networks, 'Monoamine_all_grouped', ([path 'Images_monoamine_networks/All_monoamine_grouped']), 'Mon_allgrouped_networks'); %plot and save all networks grouped
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Images_monoamine_networks/Monoamine_grouped']) %create folder to save images
names_Mon_grouped = regexprep(Mon_grouped_names, '_(?!\D)', '-');
Plot_Networks(Mon_grouped_networks, names_Mon_grouped, ([path 'Images_monoamine_networks/Monoamine_grouped']), 'Mon_networks'); %plot and save all networks grouped

%% ****************************************************************************************************************
%NEUROPEPTIDES 
%****************************************************************************************************************
%Build unrestricted individual and grouped networks
[neuropeptide_networks_all,structure_NPP,pair_names_NPP] = Neuropeptide_Networks(connections1, expression_byneuron_ordered_unweighted_all);
[~, NPP_all_grouped_networks] = Group_networks (pair_names_NPP, neuropeptide_networks_all, '\w*_\d*_\w*'); %group all networks

%Anatomically restrict networks (mid, short range)
columns_Head = [6:22, 34]; %column numbers that correspond to processes that go through the head (everything up to NR + VNC anterior, includes pharynx)
columns_Midbody = [32:33, 35:40]; %column numbers that correspond to processes that go through the midbody (all sublaterals, VNC, DNC and canal)
columns_Tail = [44:46]; %column numbers that correspond to processes that go through the tail

[NPP_all_grouped_networks_individual_short_range, NPP_all_grouped_networks_individual_mid_range, class_list_All_grouped, class_names_All_grouped, class_Head, class_Midbody, class_Tail] = Obtain_local_networks(NPP_all_grouped_networks, anatomy, {columns_Head, columns_Midbody, columns_Tail}); %obtain the short and mid range networks for the all grouped neuropeptide networks
[Neuropeptide_networks_individual_short_range, Neuropeptide_networks_individual_mid_range, class_list_neuropeptide_networks, class_names_neuropeptide_networks] = Obtain_local_networks(neuropeptide_networks_all, anatomy, {columns_Head, columns_Midbody, columns_Tail}); %obtain the short and mid range networks for the individual neuropeptide networks

NPP_individual_networks_short_range_connections = cell(size(neuropeptide_networks_all)); %empty cell array to store individual networks matrices
NPP_individual_networks_mid_range_connections = cell(size(neuropeptide_networks_all)); %empty cell array to store individual networks matrices
NPP_individual_networks_nerve_ring_connections = cell(size(neuropeptide_networks_all)); %empty cell array to store individual networks matrices
 
for i = 1:length(neuropeptide_networks_all)
    NPP_individual_networks_short_range_connections{i,1} = Group_local_networks(neuropeptide_networks_all{i}, {Neuropeptide_networks_individual_short_range{i}{2:17}, Neuropeptide_networks_individual_short_range{i}{27:35}, Neuropeptide_networks_individual_short_range{i}{39:41}}); %select short range networks for all individual neuropeptide networks
    NPP_individual_networks_mid_range_connections{i,1} = Group_local_networks(neuropeptide_networks_all{i}, {Neuropeptide_networks_individual_mid_range{i}{1}, Neuropeptide_networks_individual_short_range{i}{17}, Neuropeptide_networks_individual_mid_range{i}{2}, Neuropeptide_networks_individual_mid_range{i}{3}}); %select mid range networks for all individual neuropeptide networks
    NPP_individual_networks_nerve_ring_connections{i,1} = Group_local_networks(neuropeptide_networks_all{i}, {Neuropeptide_networks_individual_short_range{i}{16}}); %select nerve ring network for all individual neuropeptide networks
end

%Plot individual networks 
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Images_neuropeptide_networks_unweighted/Individual_NPP_networks'])  %create folder to save images
names_NPP = regexprep(pair_names_NPP, {'_(?!\D)', '_(?!\d)'},{'-', ' '});
Plot_Networks(neuropeptide_networks_all, names_NPP, ([path 'Images_neuropeptide_networks_unweighted/Individual_NPP_networks']), 'NPP_networks_LR'); %plot and save all adjacency matrices from neuropeptide_networks
Plot_Networks(NPP_individual_networks_short_range_connections, names_NPP, ([path 'Images_neuropeptide_networks_unweighted/Individual_NPP_networks']), 'NPP_networks_SR'); %plot and save all adjacency matrices from neuropeptide_networks short range
Plot_Networks(NPP_individual_networks_mid_range_connections, names_NPP, ([path 'Images_neuropeptide_networks_unweighted/Individual_NPP_networks']), 'NPP_networks_MR'); %plot and save all adjacency matrices from neuropeptide_networks short range

%Group networks NPP
neuropeptides = matlab.lang.makeValidName(unique(connections1.GeneLigand, 'stable')); %create list of the neuropeptides uses in the edges to group based on them 
receptors_NPP = matlab.lang.makeValidName(unique(connections1.GeneGPCR, 'stable')); %create list of the receptors used in the edges to group based on them, use stable to prevent unique from reordering alphabetically so that we group them in the order of connections

NPP_all_grouped_networks_mid_range_connections = Group_local_networks(NPP_all_grouped_networks, {NPP_all_grouped_networks_individual_mid_range{3}, NPP_all_grouped_networks_individual_short_range{17}, NPP_all_grouped_networks_individual_mid_range{1}, NPP_all_grouped_networks_individual_mid_range{2}}); %group all networks that have short or mid range processes (short_range{15} is the Pharynx, that only has short range processes, the other mid range ones are Head, Tail and Midbody) 
NPP_all_grouped_networks_short_range_connections = Group_local_networks(NPP_all_grouped_networks, {NPP_all_grouped_networks_individual_short_range{2:17}, NPP_all_grouped_networks_individual_short_range{27:35}, NPP_all_grouped_networks_individual_short_range{39:41}}); %group all local networks that have close process bundles

[NPP_grouped_names, NPP_grouped_networks_LR] = Group_networks (pair_names_NPP, neuropeptide_networks_all, neuropeptides); %group based on neuropeptides 
[NPGPCR_grouped_names, NPGPCR_grouped_networks_LR] = Group_networks (pair_names_NPP, neuropeptide_networks_all, receptors_NPP); %group based on neuropeptide receptors 
[~, NLP_NPP_grouped_networks_LR] = Group_networks (pair_names_NPP, neuropeptide_networks_all, '(NLP|SNET|PDF|NTC)_\d*_'); %group based on NLP neuropeptides 
[~, FLP_NPP_grouped_networks_LR] = Group_networks (pair_names_NPP, neuropeptide_networks_all, 'FLP_\d*_'); %group based on FLP neuropeptides 
[~, NPP_all_grouped_networks_promiscuous_LR] = Group_networks (pair_names_NPP, neuropeptide_networks_all, '\w*_\d*_(DMSR_7|DMSR_1|FRPR_8)'); %group all networks removing the most dense receptors
[~, NPP_all_grouped_networks_nopromiscuous_LR] = Group_networks (pair_names_NPP, neuropeptide_networks_all, '\w*_\d*_(?!DMSR_7|DMSR_1|FRPR_8)'); %group all networks removing the most dense receptors
[~, NPP_grouped_networks_MR] = Group_networks (pair_names_NPP, NPP_individual_networks_mid_range_connections, neuropeptides); %group based on neuropeptides 
[~, NPGPCR_grouped_networks_MR] = Group_networks (pair_names_NPP, NPP_individual_networks_mid_range_connections, receptors_NPP); %group based on neuropeptide receptors 
[~, NLP_NPP_grouped_networks_MR] = Group_networks (pair_names_NPP, NPP_individual_networks_mid_range_connections, '(NLP|SNET|PDF|NTC)_\d*_'); %group based on NLP neuropeptides 
[~, FLP_NPP_grouped_networks_MR] = Group_networks (pair_names_NPP, NPP_individual_networks_mid_range_connections, 'FLP_\d*_'); %group based on FLP neuropeptides 
[~, NPP_all_grouped_networks_promiscuous_MR] = Group_networks (pair_names_NPP, NPP_individual_networks_mid_range_connections, '\w*_\d*_(DMSR_7|DMSR_1|FRPR_8)'); %group all networks removing the most dense receptors
[~, NPP_all_grouped_networks_nopromiscuous_MR] = Group_networks (pair_names_NPP, NPP_individual_networks_mid_range_connections, '\w*_\d*_(?!DMSR_7|DMSR_1|FRPR_8)'); %group all networks removing the most dense receptors
[~, NPP_grouped_networks_SR] = Group_networks (pair_names_NPP, NPP_individual_networks_short_range_connections, neuropeptides); %group based on neuropeptides 
[~, NPGPCR_grouped_networks_SR] = Group_networks (pair_names_NPP, NPP_individual_networks_short_range_connections, receptors_NPP); %group based on neuropeptide receptors 
[~, NLP_NPP_grouped_networks_SR] = Group_networks (pair_names_NPP, NPP_individual_networks_short_range_connections, '(NLP|SNET|PDF|NTC)_\d*_'); %group based on NLP neuropeptides 
[~, FLP_NPP_grouped_networks_SR] = Group_networks (pair_names_NPP, NPP_individual_networks_short_range_connections, 'FLP_\d*_'); %group based on FLP neuropeptides 
[~, NPP_all_grouped_networks_promiscuous_SR] = Group_networks (pair_names_NPP, NPP_individual_networks_short_range_connections, '\w*_\d*_(DMSR_7|DMSR_1|FRPR_8)'); %group all networks removing the most dense receptors
[~, NPP_all_grouped_networks_nopromiscuous_SR] = Group_networks (pair_names_NPP, NPP_individual_networks_short_range_connections, '\w*_\d*_(?!DMSR_7|DMSR_1|FRPR_8)'); %group all networks removing the most dense receptors

%Plot grouped networks
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Images_neuropeptide_networks_unweighted/NPP_all_grouped_networks'])  %create folder to save images
Plot_Networks(NPP_all_grouped_networks, 'NPP_all_grouped_networks', ([path 'Images_neuropeptide_networks_unweighted/NPP_all_grouped_networks']), 'NPP_allgroupednetworks_LR'); %plot and save all networks grouped
Plot_Networks(NPP_all_grouped_networks_mid_range_connections, 'All_NPP_grouped_networks_mid_range_connections', ([path 'Images_neuropeptide_networks_unweighted/NPP_all_grouped_networks']), 'NPP_allgroupednetworks_MR');%plot and save NPP all grouped networks for neurons that have connections with neurons close to them
Plot_Networks(NPP_all_grouped_networks_short_range_connections, 'All_NPP_grouped_networks_short_range_connections', ([path 'Images_neuropeptide_networks_unweighted/NPP_all_grouped_networks']), 'NPP_allgroupednetworks_SR');%plot and save NPP all grouped networks for neurons that have connections with neurons close to them

mkdir ([ 'Results/scRNAseq_C.elegans' path 'Images_neuropeptide_networks_unweighted/NPP_grouped_networks'])  %create folder to save images
names_NPP_grouped = regexprep(NPP_grouped_names, '_(?!\D)', '-');
Plot_Networks(NPP_grouped_networks_LR, names_NPP_grouped, ([path 'Images_neuropeptide_networks_unweighted/NPP_grouped_networks']), 'NPP_groupnetworks_LR'); %plot and save NPP grouped networks
Plot_Networks(NPP_grouped_networks_MR, names_NPP_grouped, ([path 'Images_neuropeptide_networks_unweighted/NPP_grouped_networks']), 'NPP_groupnetworks_MR'); %plot and save NPP grouped networks
Plot_Networks(NPP_grouped_networks_SR, names_NPP_grouped, ([path 'Images_neuropeptide_networks_unweighted/NPP_grouped_networks']), 'NPP_groupnetworks_SR'); %plot and save NPP grouped networks
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Images_neuropeptide_networks_unweighted/NPGPCR_grouped_networks'])  %create folder to save images
names_NPPGPCR_grouped = regexprep(NPGPCR_grouped_names, '_(?!\D)', '-');
Plot_Networks(NPGPCR_grouped_networks_LR, names_NPPGPCR_grouped, ([path 'Images_neuropeptide_networks_unweighted/NPGPCR_grouped_networks']), 'GPCR_networks_LR');%plot and save NPGPCR grouped networks
Plot_Networks(NPGPCR_grouped_networks_MR, names_NPPGPCR_grouped, ([path 'Images_neuropeptide_networks_unweighted/NPGPCR_grouped_networks']), 'GPCR_networks_MR');%plot and save NPGPCR grouped networks
Plot_Networks(NPGPCR_grouped_networks_SR, names_NPPGPCR_grouped, ([path 'Images_neuropeptide_networks_unweighted/NPGPCR_grouped_networks']), 'GPCR_networks_SR');%plot and save NPGPCR grouped networks
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Images_neuropeptide_networks_unweighted/NLP_NPP_grouped_networks'])  %create folder to save images
Plot_Networks(NLP_NPP_grouped_networks_LR, 'NLP_NPP_grouped_networks', ([path 'Images_neuropeptide_networks_unweighted/NLP_NPP_grouped_networks']), 'NLP_networks_LR');%plot and save NPGPCR grouped networks
Plot_Networks(NLP_NPP_grouped_networks_MR, 'NLP_NPP_grouped_networks', ([path 'Images_neuropeptide_networks_unweighted/NLP_NPP_grouped_networks']), 'NLP_networks_MR');%plot and save NPGPCR grouped networks
Plot_Networks(NLP_NPP_grouped_networks_SR, 'NLP_NPP_grouped_networks', ([path 'Images_neuropeptide_networks_unweighted/NLP_NPP_grouped_networks']), 'NLP_networks_SR');%plot and save NPGPCR grouped networks
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Images_neuropeptide_networks_unweighted/FLP_NPP_grouped_networks'])  %create folder to save images
Plot_Networks(FLP_NPP_grouped_networks_LR, 'FLP_NPP_grouped_networks', ([path 'Images_neuropeptide_networks_unweighted/FLP_NPP_grouped_networks']), 'FLP_networks_LR');%plot and save NPGPCR grouped networks
Plot_Networks(FLP_NPP_grouped_networks_MR, 'FLP_NPP_grouped_networks', ([path 'Images_neuropeptide_networks_unweighted/FLP_NPP_grouped_networks']), 'FLP_networks_MR');%plot and save NPGPCR grouped networks
Plot_Networks(FLP_NPP_grouped_networks_SR, 'FLP_NPP_grouped_networks', ([path 'Images_neuropeptide_networks_unweighted/FLP_NPP_grouped_networks']), 'FLP_networks_SR');%plot and save NPGPCR grouped networks
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Images_neuropeptide_networks_unweighted/NPP_all_noprom_grouped_networks'])  %create folder to save images
Plot_Networks(NPP_all_grouped_networks_nopromiscuous_LR, 'NPP_all_noprom_grouped_networks', ([path 'Images_neuropeptide_networks_unweighted/NPP_all_noprom_grouped_networks']), 'NPP_all_noprom_groupednetworks_LR'); %plot and save all networks grouped
Plot_Networks(NPP_all_grouped_networks_nopromiscuous_MR, 'NPP_all_noprom_grouped_networks', ([path 'Images_neuropeptide_networks_unweighted/NPP_all_noprom_grouped_networks']), 'NPP_all_noprom_groupednetworks_MR'); %plot and save all networks grouped
Plot_Networks(NPP_all_grouped_networks_nopromiscuous_SR, 'NPP_all_noprom_grouped_networks', ([path 'Images_neuropeptide_networks_unweighted/NPP_all_noprom_grouped_networks']), 'NPP_all_noprom_groupednetworks_SR'); %plot and save all networks grouped
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Images_neuropeptide_networks_unweighted/NPP_all_prom_grouped_networks'])  %create folder to save images
Plot_Networks(NPP_all_grouped_networks_promiscuous_LR, 'NPP_all_prom_grouped_networks', ([path 'Images_neuropeptide_networks_unweighted/NPP_all_prom_grouped_networks']), 'NPP_all_prom_groupednetworks_LR'); %plot and save all networks grouped
Plot_Networks(NPP_all_grouped_networks_promiscuous_MR, 'NPP_all_prom_grouped_networks', ([path 'Images_neuropeptide_networks_unweighted/NPP_all_prom_grouped_networks']), 'NPP_all_prom_groupednetworks_MR'); %plot and save all networks grouped
Plot_Networks(NPP_all_grouped_networks_promiscuous_SR, 'NPP_all_prom_grouped_networks', ([path 'Images_neuropeptide_networks_unweighted/NPP_all_prom_grouped_networks']), 'NPP_all_prom_groupednetworks_SR'); %plot and save all networks grouped


%% ****************************************************************************************************************
%SAVE ALL RESULTS FROM SCRIPT
%****************************************************************************************************************
%Export CSV of networks

mkdir ([ 'Results/scRNAseq_C.elegans' path 'CSV_networks']) %create folder to save data
mkdir ([ 'Results/scRNAseq_C.elegans' path 'CSV_networks/Individual_SR']) %create folder to save data
mkdir ([ 'Results/scRNAseq_C.elegans' path 'CSV_networks/Individual_MR']) %create folder to save data
mkdir ([ 'Results/scRNAseq_C.elegans' path 'CSV_networks/Individual_LR']) %create folder to save data
mkdir ([ 'Results/scRNAseq_C.elegans' path 'CSV_networks/NPP_grouped_SR']) %create folder to save data
mkdir ([ 'Results/scRNAseq_C.elegans' path 'CSV_networks/NPP_grouped_MR']) %create folder to save data
mkdir ([ 'Results/scRNAseq_C.elegans' path 'CSV_networks/NPP_grouped_LR']) %create folder to save data

writetable(array2table(Mon_all_grouped_networks, 'RowNames', neuron_names, 'VariableNames', neuron_names), ...
    ['Results/scRNAseq_C.elegans' path 'CSV_networks/' date 'monoamine_connectome.csv'], 'WriteRowNames',true); %save mon aggregate network
writetable(array2table(NPP_all_grouped_networks, 'RowNames', neuron_names, 'VariableNames', neuron_names), ...
    ['Results/scRNAseq_C.elegans' path 'CSV_networks/' date 'neuropeptide_connectome_long_range_model.csv'], 'WriteRowNames',true); %save long-range aggregate network
writetable(array2table(NPP_all_grouped_networks_mid_range_connections, 'RowNames', neuron_names, 'VariableNames', neuron_names), ...
    ['Results/scRNAseq_C.elegans' path 'CSV_networks/' date 'neuropeptide_connectome_mid_range_model.csv'], 'WriteRowNames',true); %save mid-range aggregate network
writetable(array2table(NPP_all_grouped_networks_short_range_connections, 'RowNames', neuron_names, 'VariableNames', neuron_names), ...
    ['Results/scRNAseq_C.elegans' path 'CSV_networks/' date 'neuropeptide_connectome_short_range_model.csv'], 'WriteRowNames',true); %save short-range aggregate network
for i = 1:(length(neuropeptide_networks_all))
writetable(array2table(neuropeptide_networks_all{i}, 'RowNames', neuron_names, 'VariableNames', neuron_names), ...
    strcat('Results/scRNAseq_C.elegans', path, 'CSV_networks/Individual_LR/', date, 'neuropeptide_network', num2str(i, '%03d'), '_long_range_model.csv'), 'WriteRowNames',true); %save short-range aggregate network
writetable(array2table(NPP_individual_networks_mid_range_connections{i}, 'RowNames', neuron_names, 'VariableNames', neuron_names), ...
    strcat('Results/scRNAseq_C.elegans', path, 'CSV_networks/Individual_MR/', date, 'neuropeptide_network', num2str(i, '%03d'), '_mid_range_model.csv'), 'WriteRowNames',true); %save short-range aggregate network
writetable(array2table(NPP_individual_networks_short_range_connections{i}, 'RowNames', neuron_names, 'VariableNames', neuron_names), ...
    strcat('Results/scRNAseq_C.elegans', path, 'CSV_networks/Individual_SR/', date, 'neuropeptide_network', num2str(i, '%03d'), '_short_range_model.csv'), 'WriteRowNames',true); %save short-range aggregate network
end
for i = 1:(length(NPP_grouped_names))
writetable(array2table(NPP_grouped_networks_LR{i}, 'RowNames', neuron_names, 'VariableNames', neuron_names), ...
    strcat('Results/scRNAseq_C.elegans', path, 'CSV_networks/NPP_grouped_LR/', date, 'neuropeptide_grouped_network', num2str(i, '%03d'), '_long_range_model.csv'), 'WriteRowNames',true); %save short-range aggregate network
writetable(array2table(NPP_grouped_networks_MR{i}, 'RowNames', neuron_names, 'VariableNames', neuron_names), ...
    strcat('Results/scRNAseq_C.elegans', path, 'CSV_networks/NPP_grouped_MR/', date, 'neuropeptide_grouped_network', num2str(i, '%03d'), '_mid_range_model.csv'), 'WriteRowNames',true); %save short-range aggregate network
writetable(array2table(NPP_grouped_networks_SR{i}, 'RowNames', neuron_names, 'VariableNames', neuron_names), ...
    strcat('Results/scRNAseq_C.elegans', path, 'CSV_networks/NPP_grouped_SR/', date, 'neuropeptide_grouped_network', num2str(i, '%03d'), '_short_range_model.csv'), 'WriteRowNames',true); %save short-range aggregate network
end

%Save workspace
mkdir(['Results/scRNAseq_C.elegans' path2]); %create folder
save((['Results/scRNAseq_C.elegans' path2 date 'Networks_build_EC500nM.mat'])) %save all variables from workspace to results at this point
end

function Plot_Networks(data, names_list, folder, name)

%Plot adjacency matrices and distribution of weights created by Neuropeptide_Networks in structure format and save
if iscell(data)
    numtiles = 12*8;
    numnet = length(names_list);
    nf = ceil(numnet/numtiles);
    for j = 1: nf
        fig = figure('Visible','off');
        t1 = tiledlayout(12,8,'TileSpacing','compact', 'Padding', 'tight', 'PositionConstraint','innerposition');
        start = (j-1)*(numtiles)+1;
        finish = min(start + numtiles - 1, numnet);
    %Plot the adjacency matrix 
        for i = start:finish
            nexttile(t1)
            spy(data{i}, 1); %plots adjacency matrix as a sparse plot
            title(names_list{i}, 'interpreter', 'none'); %adds title
            subtitle(sprintf('nz= %d', nnz(data{i})))
            xlabel('')
            box on
            set(gca, 'TickLength', [0 0], 'YTickLabel', [], 'XTickLabel', [])
            hold on
            yline(20);
            yline(99);
            yline(181);
            xline(20);
            xline(99);
            xline(181);
            hold off
        end
        t1.XLabel.String = 'Receiving neurons';
        t1.YLabel.String = 'Sending neurons';
        set(findall(gcf,'-property','FontSize'),'FontSize',5);
        print(['Results/scRNAseq_C.elegans/' folder '/Amatrix_' name '_' num2str(j)], '-dpdf', '-fillpage', '-r300');%save image as pdf
    %Plot the distribution of the nonzero values of the adjacency matrix if its interesting
        h = max(cell2mat(cellfun(@max, cellfun(@max, data, 'UniformOutput', false), 'UniformOutput', false))); %get the maximum of all the matrices in the cell array to set colorbar maximum
        if h == 1
            continue
        elseif h > 1
            clf %erase previous figure
            t2 = tiledlayout(12,8,'TileSpacing','Compact', 'Padding', 'tight', 'PositionConstraint','innerposition');
            for i = 1: length(names_list)
                nexttile(t2)
                cspy_Lidia2(data{i}, 'markersize', 1, 'ColorMap', 'jet', 'maxval', h); %plots a spy plot but colours the edges by weight and has colorbar with maxval set for the figure
                colorbar off
                title(names_list{i}, 'interpreter', 'none'); %adds title
                box on
                set(gca, 'TickLength', [0 0], 'YTickLabel', [], 'XTickLabel', [])
                hold on
                yline(20);
                yline(99);
                yline(181);
                xline(20);
                xline(99);
                xline(181);
                hold off
            end
            cb = colorbar;
            cb.Layout.TileSpan = [1 6];
            cb.Layout.Tile = 'east';
            t2.XLabel.String = 'Receiving neurons';
            t2.YLabel.String = 'Sending neurons';
            set(findall(gcf,'-property','FontSize'),'FontSize',5)
            print(['Results/scRNAseq_C.elegans/' folder '/Weights_' name '_' num2str(j)], '-dpdf', '-fillpage', '-r300')
            close(fig) %clean memory
        end
    end
else
  %Plot sparse weighted adjancency matrices
    figure('Visible', 'off');
    spy(data); %plots adjacency matrix as a sparse plot
    title(names_list, 'interpreter', 'none'); %adds title
    xlabel ('Receiving neurons')
    ylabel ('Sender neurons')
    box on
    set(gca, 'TickLength', [0 0], 'YTickLabel', [], 'XTickLabel', [])
    hold on
    yline(20);
    yline(99);
    yline(181);
    xline(20);
    xline(99);
    xline(181);
    hold off
    print(['Results/scRNAseq_C.elegans/' folder '/Amatrix_' name],  '-dpdf', '-fillpage', '-r300');%save image as pdf
  %Plot weighted adjancency matrix in cspy heatmap
    figure ('Visible', 'off', 'Position', [109.8035,94.38999999999999,648.6362499999999,691.935]);
    cspy_Lidia2(data, 'markersize', 5, 'ColorMap', 'jet', 'Levels', [5]); %plots a spy plot but colours the edges by regions of weight(that is why we select level 5, so that there are 5 colours, weights are split in 5 categories, we choose 5 because we saw it was a good number for the data) and has colorbar
    title(names_list, 'interpreter', 'none'); %adds title
    xlabel ('Receiving neurons')
    ylabel ('Sending neurons')
    box on
    set(gca, 'TickLength', [0 0], 'YTickLabel', [], 'XTickLabel', [])
    hold on
    yline(20);
    yline(99);
    yline(181);
    xline(20);
    xline(99);
    xline(181);
    hold off
    colorbar;
    print(['Results/scRNAseq_C.elegans/' folder '/Weights_' name], '-dpdf', '-fillpage', '-r300');
end
end
function [neuropeptide_network, structure, pair_names] = Neuropeptide_Networks (edges, data)
structure = struct;%creates structure to store adjacency matrices
neuropeptide_network = cell(height(edges),1); %creates cell to store adjacency matrices
pair_names = cell(height(edges), 1); %creates cell to store pair names
    for i = 1:height(edges)
        neuropeptide = char(edges.('GeneLigand')(i)); %select one by one each neuropeptide
        receptor = char(edges.('GeneGPCR')(i)); %select the corresponding receptor
        name = append (matlab.lang.makeValidName(neuropeptide), '_', matlab.lang.makeValidName(receptor)); %name for every field in the structure, so every network
        pair_names{i} = name; %saves cell with all the neuropeptide_receptor pairs, index for structure
        neuropeptide_network{i} = buildnetwork (neuropeptide, receptor, data); %create a new dimension of matrix for each iteration
        structure.(name) = buildnetwork (neuropeptide, receptor, data); %create a field in the structure for each iteration
    end
end
function N = neuropeptideexpression (neuropeptide, data)
% Create expression matrix for each neuropeptide and receptor
N = table2array(data(strcmpi(data.(1), neuropeptide), 2:width(data)))';  % get data for neuropeptide value, transform to matrix and transpose to make it be a column, so NPP expression is in the matrix rows
end
function R = receptorexpression (receptor, data)
% Create expression matrix for each neuropeptide and receptor
R = table2array(data(strcmpi(data.(1), receptor), 2:width(data))); %get data for receptor value, transform to matrix and leave as a row, so that receptor expression is in the matrix columns 
end
function B = buildnetwork (neuropeptide, receptor, data)
%Create adjacency matrix for each NPP NPGPCR pair
neExp = neuropeptideexpression(neuropeptide, data); %vector with expression values for neuropeptides
reExp = receptorexpression(receptor, data); %vector with expression values for receptors
if isempty(neExp) == 1 || isempty(reExp) == 1 
    B = 'missing_data'; %if there are missing values for one of the neuropeptide_receptor pairs print missing data
elseif all(ismember(neExp, [0 1])) == 1 || all(ismember(reExp, [0 1])) == 1 
    B = reExp .* neExp; % calculate adjacency matrix if values in the matrix are in binary format
else
    B = reExp .* neExp; % calculate adjacency matrix
    B = full(spfun(@log10, B)); %calculate log10 of the adjacency matrix only in non zero values, and show in full matri
end
end
function [group_names, grouped_network] = Group_networks(network_names, network_used, grouping_based_on)
grouped_network = cell(length(grouping_based_on), 1); %cell array to store the new grouped matrices
if iscell(grouping_based_on)
    for i = 1:length(grouping_based_on)
        grouped = zeros (size(network_used{1}, 1),size(network_used{1}, 2)); %to store the addition of matrices per group
        group = find(~cellfun('isempty', regexpi(network_names, strcat(grouping_based_on{i},'(_|$)')))); %select all the indeces in the networks cell array that match the condition for the grouping
        for j = 1:length(group)
            network = cell2mat((network_used(group(j,:)))); % extract matrix from the cell array based on the index j given by the group indeces
            grouped = grouped + network; %add the matrix extracted in this iteration to the previous ones
        end
    grouped_network{i} = grouped; %cell array to store new matrices
    group_names = grouping_based_on; %cell array to store the names of the new matrices
    end
else
    grouped = zeros (size(network_used{1}, 1),size(network_used{1}, 2)); %to store the addition of matrices per group
    group = find(~cellfun('isempty', regexpi(network_names, grouping_based_on)));%select all the indeces in the networks cell array that match the condition for the grouping
    for j = 1:length(group)
            network = cell2mat((network_used(group(j,:)))); % extract matrix from the cell array based on the index j given by the group indeces
            grouped = grouped + network; %add the matrix extracted in this iteration to the previous ones
    end
    grouped_network = grouped; %cell array to store new matrices
    group_names = regexpi(network_names, grouping_based_on, 'match'); %cell array to store the names of the new matrices
end
end
function [array_networks_short_range, array_networks_mid_range, class_list, class_names, class_Head, class_Midbody, class_Tail] = Obtain_local_networks (network, anatomy, class_mid_range_list)
class_list = cell(width(anatomy) - 5, 1); %empty cell array to save the anatomical neuron numbers
class_names = cell(width(anatomy) - 5, 1); %empty cell array to save the anatomical class name 
if iscell(network) 
    array_networks_short_range = cell(length(network), 1); %empty cell array to save local_networks matrices for each NPP network
    array_networks_mid_range = cell(length(network), 1); %empty cell array to save mid range local network matrices for each NPP network
    local_networks = cell(width(anatomy) - 5, 1); %empty cell array to save local matrices
    mid_range_networks = cell(length(class_mid_range_list), 1); %empty cell array to save the mid range networks
    for j = 1:length(network)
        for i = 6: width(anatomy)%loop over the columns in anatomyc_class file that have process and cell body information, from 6 onwards
          class_list{i-5} = find(anatomy.(i) > 0); %get the index of the neurons that have a value >0 in table, so that have a process in that anatommical location 
          class_names{i-5} = anatomy(:, i).Properties.VariableNames; %store the neuron indeces wiht the anatomical class list 
          local_networks{i-5} = zeros(size(network{j})); %create an empty matrix the size of the total network
          local_networks{i-5}(class_list{i-5}, class_list{i-5}) = network{j}(class_list{i-5}, class_list{i-5}); %%mask everything that is not on the neurons that have processes in the location we are interested in, by only adding the data for that region to the empty matrix 
        end 
    array_networks_short_range{j} = local_networks; %save the version of local_networks for each individual network
        for k = 1:length(class_mid_range_list)%loop over classes of mid range networks
            class_Head = unique(vertcat(class_list{class_mid_range_list{1}-5}));%select the neuron indices for the neurons that have processes in the head of the worm all bundles according to anatomical class, remove repeated indices using unique 
            class_Midbody = unique(vertcat(class_list{class_mid_range_list{2}-5}));%select neuron indices for neurons which processes are in the midbody, remove repeated indices using unique 
            class_Tail = unique(vertcat(class_list{class_mid_range_list{3}-5}));%select the neuron indices for the neurons that have processes in the tail according to anatomical class, remove repeated indices using unique
            class_mid_range = {class_Head, class_Midbody, class_Tail}; % group together the neuron indices of the neurons that have short and mid range connections
            mid_range_networks{k} = zeros(size(network{j})); %create an empty matrix the size of the total network
            mid_range_networks{k}(class_mid_range{k}, class_mid_range{k}) = network{j}(class_mid_range{k}, class_mid_range{k}); %mask everything that is not on the neurons that have processes in the location we are interested in, by only adding the data for that region to the empty matrix
        end
     array_networks_mid_range{j} = mid_range_networks; %save the version of local_networks for each individual network
    end
else
    array_networks_short_range = cell(width(anatomy) - 5, 1); %empty cell array to save local matrices 
    array_networks_mid_range = cell(length(class_mid_range_list), 1); %empty cell array to save the mid range networks
    for i = 6 : width(anatomy) %loop over anatomical class columns 
      class_list{i-5} = find(anatomy.(i) > 0); %get the index of the neurons that have a value >0 in table, so that have a process in that anatommical location 
      class_names{i-5} = anatomy(:, i).Properties.VariableNames; %store the neuron indeces wiht the anatomical class list 
      array_networks_short_range{i-5} = zeros(size(network)); %create an empty matrix the size of the total network
      array_networks_short_range{i-5}(class_list{i-5}, class_list{i-5}) = network(class_list{i-5}, class_list{i-5}); %%mask everything that is not on the neurons that have processes in the location we are interested in, by only adding the data for that region to the empty matrix 
    end
    for k = 1:length(class_mid_range_list) %loop over classes of mid range networks
      class_Head = unique(vertcat(class_list{class_mid_range_list{1}-5}));%select the neuron indices for the neurons that have processes in the head of the worm all bundles according to anatomical class, remove repeated indices using unique 
      class_Midbody = unique(vertcat(class_list{class_mid_range_list{2}-5}));%select neuron indices for neurons which processes are in the midbody, remove repeated indices using unique 
      class_Tail = unique(vertcat(class_list{class_mid_range_list{3}-5}));%select the neuron indices for the neurons that have processes in the tail according to anatomical class, remove repeated indices using unique
      class_mid_range = {class_Head, class_Midbody, class_Tail}; % group together the neuron indices of the neurons that have short and mid range connections
      array_networks_mid_range{k} = zeros(size(network)); %create an empty matrix the size of the total network
      array_networks_mid_range{k}(class_mid_range{k}, class_mid_range{k}) = network(class_mid_range{k}, class_mid_range{k}); %mask everything that is not on the neurons that have processes in the location we are interested in, by only adding the data for that region to the empty matrix
    end
end
end
function grouped_local_network = Group_local_networks(total_network, list_local_networks)
grouped_local_network = zeros(size(total_network)); %create an empty matrix the size of the total network
for i = 1: size(total_network) %loop to select the total matrix of local matrices
    for j = 1:size(total_network)
        for l = 1:length(list_local_networks)
            if grouped_local_network(i,j) == 0 && list_local_networks{1}(i,j) ~= 0 %condition to take the number of pharynx neurons connections into the total matrix if the value in the matrix is 0 (to prevent overlaps between local networks, as some neurons belong to 2 local networks)
            grouped_local_network (i,j) =  list_local_networks{1} (i,j);
            elseif grouped_local_network(i,j) == 0 && list_local_networks{l}(i,j) ~= 0 %condition to take the number of nerve ring neurons connections into the total matrix if the value in the matrix is 0 (to prevent overlaps between local networks, as some neurons belong to 2 local networks)
            grouped_local_network(i,j) =  list_local_networks{l}(i,j);
            elseif grouped_local_network(i,j) ~= 0
            grouped_local_network(i,j) = grouped_local_network(i,j); %condition to take the number of pharynx neurons connections into the total matrix if the value in the matrix is 0 (to prevent overlaps between local networks, as some neurons belong to 2 local networks)
            end
        end
     end
end
end
%**************************************************************************
function Network_degree(neuron_names, path, path2, date, version)
%****************************************************************************************************************
%LOAD DATA, BINARISE AND REMOVE SELF LOOPS FROM NETWORK
%****************************************************************************************************************
%load data 
load([pwd '/Results/scRNAseq_C.elegans' path2 date 'Networks_build_EC500nM.mat'], 'Mon_all_grouped_networks', 'NPP_all_grouped_networks', ...
    'NPP_all_grouped_networks_mid_range_connections', 'NPP_all_grouped_networks_short_range_connections', ...
    'synA', 'gapA')

load([pwd '/Results/scRNAseq_C.elegans' path2 date 'Networks_measures_EC500nM.mat'], 'randNPP_NPGPCR_short_rangeu', 'randNPP_NPGPCR_mid_rangeu', ...
    'randNPP_NPGPCR_long_rangeu')

%networks in binary form and without self loops
synA_SL_u = Mask_network(synA - diag(diag(synA))); %the mask network function binarises the network and the substraction of the diagonal removes self loops
gapA_SL_u = Mask_network(gapA - diag(diag(gapA))); %the mask network function binarises the network and the substraction of the diagonal removes self loops
Mon_grouped_networks_SL_u = Mask_network(Mon_all_grouped_networks - diag(diag(Mon_all_grouped_networks))); %the mask network function binarises the network and the substraction of the diagonal removes self loops
NPP_all_grouped_networks_SL_u = Mask_network (NPP_all_grouped_networks - diag(diag(NPP_all_grouped_networks))); %the mask network function binarises the network and the substraction of the diagonal removes self loops
NPP_all_grouped_networks_mid_range_SL_u = Mask_network (NPP_all_grouped_networks_mid_range_connections - diag(diag(NPP_all_grouped_networks_mid_range_connections))); %the mask network function binarises the network and the substraction of the diagonal removes self loops
NPP_all_grouped_networks_short_range_SL_u = Mask_network (NPP_all_grouped_networks_short_range_connections - diag(diag(NPP_all_grouped_networks_short_range_connections))); %the mask network function binarises the network and the substraction of the diagonal removes self loops

%networks in weighted form without self loops 
synA_SL = synA - diag(diag(synA)); %the substraction of the diagonal removes self loops  
gapA_SL = gapA - diag(diag(gapA)); %the substraction of the diagonal removes self loops  
Mon_grouped__SL = Mon_all_grouped_networks- diag(diag(Mon_all_grouped_networks)); %the substraction of the diagonal removes self loops  
NPP_all_grouped_networks_SL = NPP_all_grouped_networks - diag(diag(NPP_all_grouped_networks)); %the substraction of the diagonal removes self loops  
NPP_all_grouped_networks_mid_range_SL = NPP_all_grouped_networks_mid_range_connections - diag(diag(NPP_all_grouped_networks_mid_range_connections)); %the substraction of the diagonal removes self loops  
NPP_all_grouped_networks_short_range_SL = NPP_all_grouped_networks_short_range_connections - diag(diag(NPP_all_grouped_networks_short_range_connections)); %the substraction of the diagonal removes self loops  

%****************************************************************************************************************
%CALCULATE
%****************************************************************************************************************
%Unweighted: I use the degree_dir.m function from the Brain Connectivity toolbox but I name out and in degree the other way around from what the algorithm states because in my matrix the
%columns are the outdegree (NPP) and the rows the indegree (NPGCPR). For the synaptic and monoamine networks the order is like in the function
%receiver neurons in the rows and releaser/send neurons in the columns. 

[indeg_Syn, outdeg_Syn, degree_Syn] = degrees_dir(synA_SL_u);%use BrainConnectivity function to get directed degrees, indegrees accounts for the synapses the node receives and outdegree for the synapses it forms
Syn_degree_nodes = table(neuron_names, transpose(outdeg_Syn), transpose(indeg_Syn), transpose(degree_Syn), 'VariableNames', {'Neuron', 'Outdegree_Syn', 'Indegree_Syn', 'Degree_Syn'}); %transform degree information in table
[degree_GJ] = degrees_und(gapA_SL_u);%use BrainConnectivity function to get undirected degrees,  accounts for the gap junctions in the node 
GJ_degree_nodes = table(neuron_names, transpose(degree_GJ), 'VariableNames', {'Neuron', 'Degree_GJ'}); %transform degree information in table
[indeg_Mon, outdeg_Mon, degree_Mon] = degrees_dir(Mon_grouped_networks_SL_u);%use BrainConnectivity function to get directed degrees, indegrees accounts for the monoamines the node receives and outdegree for the monoamines it forms
Mon_degree_nodes = table(neuron_names, transpose(outdeg_Mon), transpose(indeg_Mon), transpose(degree_Mon), 'VariableNames', {'Neuron', 'Outdegree_Mon', 'Indegree_Mon', 'Degree_Mon'}); %transform degree information in table

[indeg_NPP_all, outdeg_NPP_all, degree_NPP_all] = degrees_dir(NPP_all_grouped_networks_SL_u); %use BrainConnectivity function to get directed degrees, indegrees accounts for NPP node receives and outdegree for the NPP it releases
NPP_all_degree_nodes = table(neuron_names, transpose(outdeg_NPP_all), transpose(indeg_NPP_all), transpose(degree_NPP_all), 'VariableNames', {'Neuron', 'Outdegree_NPP_all', 'Indegree_NPP_all', 'Degree_NPP_all'}); %transform degree information in table
[indeg_NPP_all_mid_range, outdeg_NPP_all_mid_range, degree_NPP_all_mid_range] = degrees_dir(NPP_all_grouped_networks_mid_range_SL_u); %use BrainConnectivity function to get direct degrees, indegrees and outdegrees unweighted
NPP_all_grouped_networks_mid_range_connections_degree_nodes = table(neuron_names, transpose(outdeg_NPP_all_mid_range), transpose(indeg_NPP_all_mid_range), transpose(degree_NPP_all_mid_range), 'VariableNames', {'Neuron', 'Outdegree_NPP_all_mid_range', 'Indegree_NPP_all_mid_range', 'Degree_NPP_all_mid_range'}); %transform degree information in table
[indeg_NPP_all_short_range, outdeg_NPP_all_short_range, degree_NPP_all_short_range] = degrees_dir(NPP_all_grouped_networks_short_range_SL_u); %use BrainConnectivity function to get direct degrees, indegrees and outdegrees unweighted
NPP_all_grouped_networks_short_range_connections_degree_nodes = table(neuron_names, transpose(outdeg_NPP_all_short_range), transpose(indeg_NPP_all_short_range), transpose(degree_NPP_all_short_range), 'VariableNames', {'Neuron', 'Outdegree_NPP_all_short_range', 'Indegree_NPP_all_short_range', 'Degree_NPP_all_short_range'}); %transform degree information in table

% Connectome_degrees_nodes_unweighted = join(NPP_all_degree_nodes, join(NPP_all_grouped_networks_mid_range_connections_degree_nodes, join(NPP_all_grouped_networks_short_range_connections_degree_nodes, join(FLP_degree_nodes, join(FLP_lessdensereceptors_degree_nodes, join(NLP_degree_nodes, join(Syn_degree_nodes, join(GJ_degree_nodes, Mon_degree_nodes)))))))); %join all tables together into one
Connectome_degrees_nodes_unweighted = join(NPP_all_degree_nodes, join(NPP_all_grouped_networks_mid_range_connections_degree_nodes, join(NPP_all_grouped_networks_short_range_connections_degree_nodes, join(Syn_degree_nodes, join(GJ_degree_nodes, Mon_degree_nodes))))); %join all tables together into one
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Topological_data_unweighted/Topological_data_degree/Excel_tables/']) %create folder to save images
writetable(Connectome_degrees_nodes_unweighted,(['Results/scRNAseq_C.elegans' path 'Topological_data_unweighted/Topological_data_degree/Excel_tables/' date 'scRNAseq_Celegans_connectome_degree_unweighted.csv']),'Delimiter',',','QuoteStrings',true)  %save overall connectome table as a csv file to be able to further analyse in excel

%****************************************************************************************************************
%PLOT
%****************************************************************************************************************

%sort data
[degree_SR_sorted, index1] = sort(degree_NPP_all_short_range, 'descend'); 
[degree_MR_sorted, index2] = sort(degree_NPP_all_mid_range, 'descend'); 
[degree_LR_sorted, index3] = sort(degree_NPP_all, 'descend'); 
[degree_Syn_sorted, index4] = sort(degree_Syn, 'descend'); 
[degree_Mon_sorted, index5] = sort(degree_Mon, 'descend'); 

%get plots
figure('Renderer', 'painters', 'Visible', 'off', 'PaperOrientation', 'portrait')
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'tight');
nexttile;
plot(degree_SR_sorted, 'Color', 1/255*[119 172 48], 'LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor',1/255*[119 172 48], 'MarkerSize', 2);
xlim([0 size(degree_SR_sorted, 2)+5])
ylim([0 600])
hold on
plot(indeg_NPP_all_short_range(index1), 'Color', 1/255*[0 114 189], 'LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor',1/255*[0 114 189], 'MarkerSize', 2);
plot(outdeg_NPP_all_short_range(index1), 'Color', 1/255*[237 177 32], 'LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor',1/255*[237 177 32], 'MarkerSize', 2);
yline(600,'k')
xline(size(degree_SR_sorted, 2)+5,'k')
box off
xlabel('Individual neurons', 'FontSize', 9)
title('Short Range', 'FontSize', 9)
axis square
hold off
set(gca, 'TickDir','out', 'XTick', [], 'XTickLabel', [], 'YColor', 'k', 'XColor', 'k', 'FontName','Arial'); %remove tick marks and prevent unwanted underscores in tick labels
nexttile
plot(degree_MR_sorted, 'Color', 1/255*[119 172 48], 'LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor',1/255*[119 172 48], 'MarkerSize', 2);
xlim([0 size(degree_MR_sorted, 2)+5])
ylim([0 600])
hold on
plot(indeg_NPP_all_mid_range(index2), 'Color', 1/255*[0 114 189], 'LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor',1/255*[0 114 189], 'MarkerSize', 2);
plot(outdeg_NPP_all_mid_range(index2), 'Color', 1/255*[237 177 32], 'LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor',1/255*[237 177 32], 'MarkerSize', 2);
yline(600,'k')
xline(size(degree_MR_sorted, 2)+5,'k')
box off
xlabel('Individual neurons', 'FontSize', 9)
title('Mid Range', 'FontSize', 9)
axis square
hold off
set(gca, 'TickDir','out', 'XTick', [], 'XTickLabel', [], 'YColor', 'k', 'XColor', 'k', 'FontName','Arial'); %remove tick marks and prevent unwanted underscores in tick labels
nexttile
plot(degree_Syn_sorted, 'Color', 1/255*[119 172 48], 'LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor',1/255*[119 172 48], 'MarkerSize', 2);
xlim([0 size(degree_Syn_sorted, 2)+5])
ylim([0 max(degree_Syn_sorted)+10])
hold on
plot(indeg_Syn(index4), 'Color', 1/255*[0 114 189], 'LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor',1/255*[0 114 189], 'MarkerSize', 2);
plot(outdeg_Syn(index4), 'Color', 1/255*[237 177 32], 'LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor',1/255*[237 177 32], 'MarkerSize', 2);
yline(max(degree_Syn_sorted)+10,'k')
xline(size(degree_Syn_sorted, 2)+5,'k')
box off
xlabel('Individual neurons', 'FontSize', 9)
title('Synaptic', 'FontSize', 9)
axis square
hold off
set(gca, 'TickDir','out', 'XTick', [], 'XTickLabel', [], 'YColor', 'k', 'XColor', 'k', 'FontName','Arial'); %remove tick marks and prevent unwanted underscores in tick labels
nexttile
plot(degree_Mon_sorted, 'Color', 1/255*[119 172 48], 'LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor',1/255*[119 172 48], 'MarkerSize', 2);
xlim([0 size(degree_Mon_sorted, 2)+5])
ylim([0 max(degree_Mon_sorted)+20])
hold on
plot(indeg_Mon(index5), 'Color', 1/255*[0 114 189], 'LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor',1/255*[0 114 189], 'MarkerSize', 2);
plot(outdeg_Mon(index5), 'Color', 1/255*[237 177 32], 'LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor',1/255*[237 177 32], 'MarkerSize', 2);
yline(max(degree_Mon_sorted)+20,'k')
xline(size(degree_Mon_sorted, 2)+5,'k')
box off
xlabel('Individual neurons', 'FontSize', 9)
title('Monoamine', 'FontSize', 9)
set(gca, 'TickDir','out', 'XTick', [], 'XTickLabel', [], 'YColor', 'k', 'XColor', 'k', 'FontName','Arial'); %remove tick marks and prevent unwanted underscores in tick labels
axis square
hold off
set(findall(gcf,'-property','FontSize'),'FontSize',9);
print(['Results/scRNAseq_C.elegans/' path 'Topological_data_unweighted/Topological_data_degree/' date 'Degree_plots'], '-dpdf', '-fillpage', '-r300');%save image as pdf

%****************************************************************************************************************
%GET RICH CLUB
%****************************************************************************************************************
[phi_SR, phi_rand_mean_SR, phi_norm_SR] = Network_computeRichClub(NPP_all_grouped_networks_short_range_SL_u, randNPP_NPGPCR_short_rangeu, 'All grouped neuropeptide networks short range rich club', (['Results/scRNAseq_C.elegans' path 'Topological_data_unweighted/Rich_club']), neuron_names.', @rich_club_bd, 'deg');
[phi_MR, phi_rand_mean_MR, phi_norm_MR] = Network_computeRichClub(NPP_all_grouped_networks_mid_range_SL_u, randNPP_NPGPCR_mid_rangeu, 'All grouped neuropeptide networks mid range rich club', (['Results/scRNAseq_C.elegans' path 'Topological_data_unweighted/Rich_club']), neuron_names.', @rich_club_bd, 'deg');
[phi_LR, phi_rand_mean_LR, phi_norm_LR] = Network_computeRichClub(NPP_all_grouped_networks_SL_u, randNPP_NPGPCR_long_rangeu, 'All grouped neuropeptide networks long range rich club', (['Results/scRNAseq_C.elegans' path 'Topological_data_unweighted/Rich_club']), neuron_names.', @rich_club_bd, 'deg');

%****************************************************************************************************************
%PLOT CORRELATION
%****************************************************************************************************************
figure('Renderer', 'painters', 'Visible', 'off', 'PaperOrientation', 'landscape')
scatter(degree_NPP_all_short_range, degree_Syn, [], 'k', 'filled');
hold on
rich = find (degree_NPP_all_short_range>= 203);
scatter(degree_NPP_all_short_range(rich), degree_Syn(rich), [], 1/255*[192 1 0], 'filled');
linear_regression = polyfit(degree_NPP_all_short_range, degree_Syn, 1); %fit polinomial line to data
equation = sprintf('y = %.3f x + %.3f', linear_regression(1), linear_regression(2)); %calculate equation for polinomial fit, transform to string
h = lsline;
set(h(2), 'Color', 'k', 'LineWidth', 3)
set(h(1), 'Color', 'w', 'LineWidth', 0.5)
xlim([0 max(degree_NPP_all_short_range)+20])
ylim([0 max(degree_Syn)+5])
yline(max(degree_Syn)+5,'k')
xline(max(degree_NPP_all_short_range)+20,'k')
hold off
close(figure) %clean memory
print(['Results/scRNAseq_C.elegans/' path 'Topological_data_unweighted/Topological_data_degree/' date 'Correlation_plot'], '-dpdf', '-fillpage', '-r300');%save image as pdf

%****************************************************************************************************************
%SAVE
%****************************************************************************************************************

%Save workspace
mkdir(['Results/scRNAseq_C.elegans' path2]); %create folder
save((['Results/scRNAseq_C.elegans' path2 date 'Networks_degree_EC500nM.mat'])) %save all variables from workspace to results at this point

end

function Network_measures(neuron_names, path, path2, date, version)

load([pwd '/Results/scRNAseq_C.elegans' path2 date 'Networks_build_EC500nM.mat'], 'Mon_all_grouped_networks', 'NPP_all_grouped_networks', ...
    'NPP_all_grouped_networks_mid_range_connections', 'NPP_all_grouped_networks_short_range_connections', ...
    'NPP_individual_networks_mid_range_connections', 'NPP_individual_networks_short_range_connections', 'synA', 'gapA')

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
%DENSITY
%****************************************************************************************************************
%Calculates directed density per network. Assumes network has no self loops and doesn't take weight into account. Density is the fraction of present connections to possible connections.
[density_Syn, density_mean_rand_Syn, density_norm_Syn] = normalise_network_measure(synA_SL, randSynu, @density_dir);
[density_GJ, density_mean_rand_GJ, density_norm_GJ] = normalise_network_measure(gapA_SL, randGapu, @density_dir);
[density_Mon, density_mean_rand_Mon, density_norm_Mon] = normalise_network_measure(Mon_grouped__SL, randMonu, @density_dir);
[density_NPP_all_grouped_networks, density_mean_rand_All_NPP, density_norm_All_NPP] = normalise_network_measure(NPP_all_grouped_networks, randNPP_NPGPCR_long_rangeu, @density_dir);
[density_NPP_grouped_networks_mid_range, density_mean_rand_NPP_mid_range, density_norm_NPP_mid_range] = normalise_network_measure(NPP_all_grouped_networks_mid_range_SL, randNPP_NPGPCR_mid_rangeu, @density_dir);
[density_NPP_grouped_networks_short_range, density_mean_rand_NPP_short_range, density_norm_NPP_short_range] = normalise_network_measure(NPP_all_grouped_networks_short_range_SL, randNPP_NPGPCR_short_rangeu, @density_dir);

save((['Results/scRNAseq_C.elegans', path2 date 'Networks_measures_EC500nM.mat']))

end
function [randnetDir, randnetUndir] = RandomNet_building(network, num_rand_nets, swaps, directed)
%check if the network is a cell, to calculate random network for each element on the cell
if iscell(network) 
    randnetDir = cell(length(network), 1); %create empty cell array to store directed random networks
    randnetUndir = cell(length(network), 1); %create empty cell array to store undirected random networks
    for i = 1:length(network) %loop over every item of the cell array that contains networks
        Ad = network{i} > 0; %binarise network if it is not binarised 
        Ad = Ad - diag(diag(Ad)); %remove self loops. 
        Au = adj2simple(network{i}); %make network binarised and undirected, by removing self loops and symmetrising matrix
        for j = 1:num_rand_nets %for every random network that you want to create
            if directed %if network is directed
                randDir{j} = randmio_dir(Ad, swaps); %calculate random network of the binarised directed network with the number of swaps given
            else %if network input is undirected
                randDir = 0; %the random network directed is 0
            end
            randUndir{j} = randmio_und(Au, swaps); %calculate undirected random network
            randnetDir{i} = randDir; %store randDir for each network in the cell array
            randnetUndir{i} = randUndir; %store randUnDir for each network in the cell array
        end
    end
    
%check if the network is a matrix to calculate random network just for the matrix
else 
    randnetDir = cell(num_rand_nets, 1); %empty cell array to store directed random networks
    randnetUndir = cell(num_rand_nets, 1); %empty cell array to store undirected random networks 
    Ad = network > 0; %binarise network
    Ad = Ad - diag(diag(Ad)); %remove self loops. 
    Au = adj2simple(network); %undirect and binarise network
    
    for i = 1:num_rand_nets %for every random network that you want to create
        if directed %if network is directed
            randnetDir{i} = randmio_dir(Ad, swaps); %calculate random network of the binarised directed network with the number of swaps given
        else %if network input is undirected
            randnetDir = 0; %the random network directed is 0
        end
        randnetUndir{i} = randmio_und(Au, swaps); %calculate undirected random network
    end
end
end

function [phi, phi_rand_mean, phi_norm] = Network_computeRichClub(A, randA, networkName, savePath, neuronID, equation, type)
    
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
    [phi,Nk,Ek] = equation(A);

    %**********************************************************************
    %	NORMALISATION
    %**********************************************************************
    for i=1:length(randA)
		% Calculate rich club metric for each random graph
        phi_rand(i,1:length(phi)) = equation(randA{i});
    end
    
    % Get mean and stdev for each value of k over all random networks
	phi_rand_mean = nanmean(phi_rand);
	phi_rand_std = nanstd(phi_rand);
	
	% Normalise phi
	phi_norm = phi ./ phi_rand_mean;
    
    %**********************************************************************
    %	PLOT RICH CLUB CURVE
    %**********************************************************************
	f1 = figure('Renderer', 'painters', 'Visible', 'on');
    colororder([0 0 0 ; 1/255*[192 1 0]])
    yyaxis left
    patch([203 310 310 203], [min(phi) min(phi) max(phi) max(phi)], 1/255*[211 211 211], 'EdgeColor', 'none');
    hold on
    hLine1 = plot(phi, 'k', 'LineWidth', 2);
    hLine2 = plot(phi_rand_mean, 'Color', 1/255*[127 127 127], 'LineWidth', 2);
    errorbar(1:5:length(phi_rand_mean), phi_rand_mean(1:5:end), ...
                phi_rand_std(1:5:end), 'Color', 1/255*[127 127 127], 'LineWidth', 2);
    xlabel('Neuropeptide degree (k)','FontSize',9);
    ylabel('(\Phi) C.elegans rich-club coeff.','FontSize',9);
    set(gca,'FontSize',9,'LineWidth',1.25,'XLim', [0 length(phi)]);
    yyaxis right
    hLine3 = plot(phi_norm, 'Color', 1/255*[192 1 0], 'LineWidth', 2);
    errorbar(1:5:length(phi_rand_mean), phi_norm(1:5:end), ...
                phi_rand_std(1:5:end), 'Color', 1/255*[192 1 0], 'LineWidth', 2);
    ylabel('(\Phi) Normalised rich-club coeff.','FontSize',9);
    set(gca,'FontSize',9,'LineWidth',1.25,'XLim', [0 length(phi)]);
    yline(1,':k', 'LineWidth', 2)
    hold off
    h_legend = legend([hLine1 hLine2 hLine3], ...
            '\Phi(k)_{C. elegans}','\Phi(k)_{random}','\Phi(k)_{norm}',...
            'Location', 'northeast');
    set(h_legend,'FontSize',7, 'Location', 'southeast');
    ylabh = get(gca,'YLabel');
    set(ylabh,'Position',get(ylabh,'Position') - [2 0 0]);
    axis square
    box on
    set(findall(gcf,'-property','FontName'),'FontName','Arial');

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
    print(f1, strcat(savePath,'/Fig_RichClub_',networkName), '-dpdf', '-fillpage', '-r300');%save image as pdf
    %saveas(f1, strcat(savePath,'/Fig_RichClub_',networkName),'png');
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
function Network_dimensionality_reduction(neurons, class, path, path2, date, version)
% *************************************************************************
    %LOAD DATA
    % *********************************************************************
    load([pwd '/Results/scRNAseq_C.elegans' path2 date 'Networks_build_EC500nM.mat'], 'names_NPP', 'NPP_all_grouped_networks', ...
        'NPP_all_grouped_networks_mid_range_connections', 'NPP_all_grouped_networks_short_range_connections', ...
        'NPP_individual_networks_mid_range_connections');
    load([pwd '/Results/scRNAseq_C.elegans' path2 date 'Networks_degree_EC500nM.mat'], 'indeg_NPP_all_mid_range');

    %create directory to save plots 
    mkdir ([ 'Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted'])  %create folder to save images
    %set defaults
    set(0, 'DefaultAxesFontName', 'Arial'); %set default text
    
    %get mid range networks no SL
    NPP_individual_networks_mid_range_connections_SL = cell(length(NPP_individual_networks_mid_range_connections), 1);
    for k = 1:length(NPP_individual_networks_mid_range_connections)
        NPP_individual_networks_mid_range_connections_SL{k,1} = NPP_individual_networks_mid_range_connections{k} - diag(diag(NPP_individual_networks_mid_range_connections{k})); %get networks without self loops for each individual receptor_neuropeptide pair
    end

    %% ****************************************************************************************************************
    %NORMALISE DATA
    %****************************************************************************************************************
    
    %normalise the matrix columns using z-score, it scales (makes standard deviation 1,
    %which is essential for t-sne and not a problem for PCA) and centers 
    %(makes column means 0, which is needed for PCA and not a problem for
    %t-sne). For the incoming connections analysis we transpose first  
    
    NPP_all_grouped_networks_normal_in = zscore(transpose(NPP_all_grouped_networks),1,1); %get the zscore in the columns dimension, using population standard deviation as we consider the incoming connections of that neuron the population
    NPP_all_grouped_networks_mid_range_normal_in = zscore(transpose(NPP_all_grouped_networks_mid_range_connections),1,1); %get the zscore in the columns dimension, using population standard deviation as we consider the incoming connections of that neuron the population
    NPP_all_grouped_networks_short_range_normal_in = zscore(transpose(NPP_all_grouped_networks_short_range_connections),1,1); %get the zscore in the columns dimension, using population standard deviation as we consider the incoming connections of that neuron the population
    
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
        if tsne_NPP_all_mid_group_neurons_table.Var1(i) < -10 && tsne_NPP_all_mid_group_neurons_table.Var2(i) < 2 %select indices of neurons based on component 1 of t-sne, if larger than 0 module 1
            tsne_NPP_all_mid_group_1(i) = tsne_NPP_all_mid_index_group(i);
        elseif tsne_NPP_all_mid_group_neurons_table.Var1(i) > 15 && tsne_NPP_all_mid_group_neurons_table.Var2(i) < -5 %select indices of neurons if component 2 of t-sne is smaller than 0 and component 1 is smaller than 0 
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
    
    %plot a pie chart to indicate the number of neurons for each neuron
    %type in each group
    g0 = strcmp(tsne_NPP_all_mid_attributes_group,'none');
    g1 = strcmp(tsne_NPP_all_mid_attributes_group,'group1');
    g2 = strcmp(tsne_NPP_all_mid_attributes_group,'group2');
    g3 = strcmp(tsne_NPP_all_mid_attributes_group,'group3');
    
    figure('visible', 'off');
    tiledlayout(2, 2)
    nexttile
    pie(categorical(class(g0)));
    title('Reduced connectivity neurons');
    nexttile
    pie(categorical(class(g1)));
    title('Motorneuron hubs core');
    nexttile
    pie(categorical(class(g2)));
    title('Interneuron hubs core');
    nexttile
    pie(categorical(class(g3)));
    title('Sensory hubs core');
    set(gca, 'FontName', 'Arial', 'FontSize', 9);
    colormap([1/255*[190 186 218]; 1/255*[179 222 105]; 1/255*[128 177 211]; 1/255*[253 180 98]])
    print(['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date 'NPP_all_grouped_mid_range_pie_charts'],  '-dpdf', '-fillpage', '-r300');%save image as pdf
    close(figure);

    %plot the loadings of PCA to see which neurons affect more the clustering 
    figure('Visible', 'off'); %create figure to save
    biplot(PCA_coeff_NPP_all_mid2(:,1:2)) %plot loadings the longer the line the more it affects
    set(gca, 'FontSize', 6);
    hold on 
    scatter(PCA_coeff_NPP_all_mid2(:,1), PCA_coeff_NPP_all_mid2(:,2), [], tsne_NPP_all_mid_color_group, 'filled') %plot neurons coloured by t-sne group
    title('NPP_all_grouped_mid_range_normalised_incoming_connections_PCA_loadings', 'Interpreter', 'none');
    print(['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date 'NPP_all_grouped_mid_range_normalised_incoming_connections_PCA_loadings_biplot'],  '-dpdf', '-fillpage', '-r300');%save image as pdf
    close(figure);

    %reorder and colour spy matrix to highlight the groups
    [tsne_NPP_all_mid_attributes_group_sorted, tsne_NPP_all_mid_range_connections_group_reordered_index] = sort(tsne_NPP_all_mid_attributes_group); %set index to sort matrix based on groups
    tsne_NPP_all_mid_range_connections_group_reordered_2 = NPP_all_grouped_networks_mid_range_connections(tsne_NPP_all_mid_range_connections_group_reordered_index.', tsne_NPP_all_mid_range_connections_group_reordered_index); %sort matrix based on groups

    %plot cspy of reordered matrix 
    figure ('Visible', 'off', 'Position', [109.8035,94.38999999999999,648.6362499999999,691.935]);
    cspy_Lidia(tsne_NPP_all_mid_range_connections_group_reordered_2, 'markersize', 5, 'ColorMap', 'jet', 'Levels', [5]); %plots a spy plot but colours the edges by regions of weight(that is why we select level 5, so that there are 5 colours, weights are split in 5 categories, we choose 5 because we saw it was a good number for the data) and has colorbar
    xlabel ('Receiving neurons')
    ylabel ('Sending neurons')
    box on
    set(gca, 'TickLength', [0 0], 'YTickLabel', [], 'XTickLabel', [])
    hold on
    yline(48);
    yline(100);
    yline(155);
    xline(48);
    xline(100);
    xline(155);
    hold off
    colorbar;
    print(['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date 'NPP_all_grouped_networks_mid_range_incoming_connections_reordered_by_group_t-sne_cspy'], '-dpdf', '-fillpage', '-r300');%save image as pdf
    close(figure);

    %plot correlations with degree and group  
    figure('Visible', 'off');%create figure to plot
    violinplot(indeg_NPP_all_mid_range, tsne_NPP_all_mid_attributes_group, 'ViolinColor', [0.4940, 0.1840, 0.5560] , 'ShowData', false); %vioplinplot degrees of neurons in each module
    ylabel('Neuropeptide mid range indegree'); %label the y axes
    print(['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date 'NPP_all_grouped_networks_mid_range_incoming_connections_indegree_correlation_groups_violinplot_purple'], '-dpdf', '-fillpage', '-r300');%save image as pdf
    close(figure);

    [violinplot_indegree_p, ~, violinplot_indegree_stats] = kruskalwallis(indeg_NPP_all_mid_range, tsne_NPP_all_mid_attributes_group); %calculate correlation significance between groups
    %print(['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date 'NPP_all_grouped_networks_mid_range_violin_corrsignificance'], '-dpdf', '-fillpage', '-r300');%save image as pdf
    violinplot_indegree_pairwise = multcompare(violinplot_indegree_stats, 'CType', 'dunn-sidak'); %correct for multiple comparison
    %print(['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date 'NPP_all_grouped_networks_mid_range_violin_multiple_comparison'], '-dpdf', '-fillpage', '-r300');%save image as pdf
    close(figure);

    figure('Visible', 'off');%create figure to plot
    violinplot(degree_NPP_all_mid_range, tsne_NPP_all_mid_attributes_group, 'ViolinColor', [0.4940, 0.1840, 0.5560] , 'ShowData', false); %vioplinplot degrees of neurons in each module
    ylabel('Neuropeptide mid range degree'); %label the y axes
    print(['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date 'NPP_all_grouped_networks_mid_range_incoming_connections_degree_correlation_groups_violinplot_purple'], '-dpdf', '-fillpage', '-r300');%save image as pdf
    close(figure);

    [violinplot_degree_p, violinplot_degree_table, violinplot_degree_stats] = kruskalwallis(degree_NPP_all_mid_range, tsne_NPP_all_mid_attributes_group); %calculate correlation significance between groups
    %print(['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date 'NPP_all_grouped_networks_mid_range_violin_corrsignificance'], '-dpdf', '-fillpage', '-r300');%save image as pdf
    violinplot_degree_pairwise = multcompare(violinplot_degree_stats, 'CType', 'dunn-sidak'); %correct for multiple comparison
    % print(['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date 'NPP_all_grouped_networks_mid_range_violin_multiple_comparison'], '-dpdf', '-fillpage', '-r300');%save image as pdf
    close(figure);


    %plot correlations with degree and group  
    figure('Visible', 'on');%create figure to plot
    violinplot(outdeg_NPP_all_mid_range, tsne_NPP_all_mid_attributes_group, 'ViolinColor', [0.4940, 0.1840, 0.5560] , 'ShowData', false); %vioplinplot degrees of neurons in each module
    ylabel('Neuropeptide mid range outdegree'); %label the y axes
    print(['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date 'NPP_all_grouped_networks_mid_range_incoming_connections_outdegree_correlation_groups_violinplot_purple'], '-dpdf', '-fillpage', '-r300');%save image as pdf
    close(figure);

    [violinplot_outdegree_p, ~, violinplot_outdegree_stats] = kruskalwallis(outdeg_NPP_all_mid_range, tsne_NPP_all_mid_attributes_group); %calculate correlation significance between groups
    %print(['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date 'NPP_all_grouped_networks_mid_range_violin_corrsignificance'], '-dpdf', '-fillpage', '-r300');%save image as pdf
    violinplot_outdegree_pairwise = multcompare(violinplot_outdegree_stats, 'CType', 'dunn-sidak'); %correct for multiple comparison
    %print(['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date 'NPP_all_grouped_networks_mid_range_violin_multiple_comparison'], '-dpdf', '-fillpage', '-r300');%save image as pdf
    close(figure);

    %bar plot of the NPP-GPCR pairs that have more connections on the neurons
    %of each group
    for i = 1:length(NPP_individual_networks_mid_range_connections_SL) %loop over every neuron 
        [indeg_NPP_individual_mid_range(i,:), ~,~] = degrees_dir(NPP_individual_networks_mid_range_connections_SL{i}); %indegree of mid range connections by NPP_GPCR pairs
    end
    
    indeg_NPP_individual_mid_range_group1 = indeg_NPP_individual_mid_range(:, strcmp(tsne_NPP_all_mid_attributes_group, 'group1')); %select indegree by network for group1 neurons
    indeg_NPP_individual_mid_range_group2 = indeg_NPP_individual_mid_range(:, strcmp(tsne_NPP_all_mid_attributes_group, 'group2')); %select indegree by network for group1 neurons
    indeg_NPP_individual_mid_range_group3 = indeg_NPP_individual_mid_range(:, strcmp(tsne_NPP_all_mid_attributes_group, 'group3')); %select indegree by network for group1 neurons
    
    indeg_NPP_mid_range_group1 = sum(indeg_NPP_individual_mid_range_group1.'); %sum all the incoming connections by NPP-GPCR network
    [indeg_NPP_individual_mid_range_group1_sorted_descend, indeg_NPP_individual_mid_range_group1_sorted_index] = sort(indeg_NPP_mid_range_group1, 'descend'); %sort connections from network that has more to less
    pair_names_NPP_sorted_group1 = categorical(indeg_NPP_individual_mid_range_group1_sorted_index, [1:height(names_NPP)], names_NPP); %category of sorted NPP in the order stated by index_NPP_to_GPCR after reordering
    
    indeg_NPP_mid_range_group2 = sum(indeg_NPP_individual_mid_range_group2.'); %sum all the incoming connections by NPP-GPCR network
    [indeg_NPP_individual_mid_range_group2_sorted_descend, indeg_NPP_individual_mid_range_group2_sorted_index] = sort(indeg_NPP_mid_range_group2, 'descend'); %sort connections from network that has more to less
    pair_names_NPP_sorted_group2 = categorical(indeg_NPP_individual_mid_range_group2_sorted_index, [1:height(names_NPP)], names_NPP); %category of sorted NPP in the order stated by index_NPP_to_GPCR after reordering
    
    indeg_NPP_mid_range_group3 = sum(indeg_NPP_individual_mid_range_group3.'); %sum all the incoming connections by NPP-GPCR network
    [indeg_NPP_individual_mid_range_group3_sorted_descend, indeg_NPP_individual_mid_range_group3_sorted_index] = sort(indeg_NPP_mid_range_group3, 'descend'); %sort connections from network that has more to less
    pair_names_NPP_sorted_group3 = categorical(indeg_NPP_individual_mid_range_group3_sorted_index, [1:height(names_NPP)], names_NPP); %category of sorted NPP in the order stated by index_NPP_to_GPCR after reordering
    
    image1 = figure('Visible', 'off', 'Position', [862,103,621,859]);%create figure to plot
    tiledlayout(3,1,'TileSpacing','loose');
    nexttile
    bar(categorical(pair_names_NPP_sorted_group1(1:10), pair_names_NPP_sorted_group1(1:10)), indeg_NPP_individual_mid_range_group1_sorted_descend(1:10), 'FaceColor', [0.4660, 0.6740, 0.1880]) %plot top 20 in bar plot
    set(gca,'TickLabelInterpreter','none', 'TickLength', [0, 0], 'Fontsize', 14)%prevent text interpreter, remove axis ticks and set font size
    nexttile
    bar(categorical(pair_names_NPP_sorted_group2(1:10), pair_names_NPP_sorted_group2(1:10)), indeg_NPP_individual_mid_range_group2_sorted_descend(1:10), 'FaceColor', [0.9290, 0.6940, 0.1250])%plot top 20 in bar plot
    set(gca,'TickLabelInterpreter','none', 'TickLength', [0, 0], 'Fontsize', 14)%prevent text interpreter, remove axis ticks and set font size
    nexttile
    bar(categorical(pair_names_NPP_sorted_group3(1:10), pair_names_NPP_sorted_group3(1:10)), indeg_NPP_individual_mid_range_group3_sorted_descend(1:10), 'FaceColor', [0.4940, 0.1840, 0.5560])%plot top 20 in bar plot
    set(gca,'TickLabelInterpreter','none', 'TickLength', [0, 0], 'Fontsize', 14) %prevent text interpreter, remove axis ticks and set font size
    print(image1, ['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date 'NPP_all_grouped_networks_mid_range_incoming_connections_by_NPP-GPCR_network'], '-dpdf', '-fillpage', '-r300');%save image as pdf
    
    % %% ****************************************************************************************************************
    %LOCAL FUNCTIONS
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
        figure('Visible', 'off', 'Resize', 'off', 'Position', [199.625,350,3600,2204]); %create an image to save the figure with the right size and prevent it being resize by the plots
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
        xticks(linspace(1, width(data_GPCR), width(data_GPCR))); %add all x axis ticks
        xticklabels(neurons); %labels for x axes
        xlabel(t, 'Neurons') %x axes label
        xtickangle (90); %change orientation of labels
        set(gca, 'TickLength', [0, 0], 'TickLabelInterpreter','none', 'Fontsize', 7); %remove tick marks and prevent unwanted underscores in tick labels
        print(['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date name '_coloured_tsne_group_scatter'], '-dpdf', '-fillpage', '-r300');%save image as png
        close(figure);
    end
    %EXPORT CSV OF NETWORKS
    csvwrite(['Results/scRNAseq_C.elegans' path 'CSV_networks/' date 'neuropeptide_connectome_mid_range_sorted_by_groups'], tsne_NPP_all_mid_range_connections_group_reordered_2); %save mid-range aggregate network reordered
    writetable(table(tnse_NPP_all_mid_2, neurons, indeg_NPP_all_mid_range.', tsne_NPP_all_mid_attributes_group.', PCA_coeff_NPP_all_mid2(:,2), 'VariableNames', {'tSNE_Variables', 'Neuron', 'Indegree', 'Group', 'Biplot'}), ['Results/scRNAseq_C.elegans' path 'Dimensionality_reduction_analysis_unweighted/' date 'Groups_mid_range_network_table.csv']); %save mid-range network group data 

    %% ****************************************************************************************************************
    %SAVE ALL RESULTS FROM SCRIPT
    %****************************************************************************************************************
    save((['Results/scRNAseq_C.elegans' path2 date 'Networks_modules_EC500nM.mat'])) %save all variables from workspace to results at this point

end
%**************************************************************************
%**************************************************************************
function Network_models_comparison(connections, neurons, path, path2, date, version)

% *************************************************************************
%LOAD DATA
% *************************************************************************

%load data
NR_strata_class = readtable('21032021_Colon-Ramos_nerve_ring_neuron_strata.csv');%import list of nerve ring strata class per neuron
load([pwd '/Results/scRNAseq_C.elegans' path2 date 'Networks_build_EC500nM.mat'], 'names_NPP', 'NPP_all_grouped_networks', ...
    'NPP_all_grouped_networks_mid_range_connections', 'NPP_all_grouped_networks_short_range_connections', ...
    'NPP_individual_networks_mid_range_connections', 'NPP_individual_networks_short_range_connections', ...
    'NPP_individual_networks_nerve_ring_connections', 'neuropeptide_networks_all', 'class_Head', 'class_Midbody', 'class_Tail');

%create directory to save plots 
mkdir ([ 'Results/scRNAseq_C.elegans' path 'Comparison_btw_ranges'])  %create folder to save images

%set defaults 
set(0, 'DefaultAxesFontName', 'Arial'); %set default text

%% ************************************************************************
%CREATE MATRICES WITH CONNECTIONS PER NPP-GPCR PAIR FOR NR
% *************************************************************************

%Nerve ring adjacency matrix of NPP/GPCR connections
[NPP_GPCR_connections_NR_Amatrix, NPP_GPCR_NR_Amatrix_rows] = Connections_Amatrix (connections, NPP_individual_networks_nerve_ring_connections, neurons, names_NPP,...
    'Number of connections in nerve ring made by each of the individual NPP-GPCR pairs for each neuron node', [ path 'Comparison_btw_ranges'], date);

%Nerve ring stratum adjacency matrix of NPP/GPCR connections

NPP_individual_networks_nerve_ring_connections_strata_1 = cell(size(NPP_individual_networks_nerve_ring_connections)); %empty cell array to store matrices
NPP_individual_networks_nerve_ring_connections_strata_2 = cell(size(NPP_individual_networks_nerve_ring_connections)); %empty cell array to store matrices
NPP_individual_networks_nerve_ring_connections_strata_3 = cell(size(NPP_individual_networks_nerve_ring_connections)); %empty cell array to store matrices
NPP_individual_networks_nerve_ring_connections_strata_4 = cell(size(NPP_individual_networks_nerve_ring_connections)); %empty cell array to store matrices
NPP_individual_networks_nerve_ring_connections_strata_5 = cell(size(NPP_individual_networks_nerve_ring_connections)); %empty cell array to store matrices

NPP_individual_networks_nerve_ring_connections_strata_1(:,1) = {zeros(length(neurons),length(neurons))}; %empty variable to create a matrix the size of the whole connectome for each NPP_GPCR network
NPP_individual_networks_nerve_ring_connections_strata_2(:,1) = {zeros(length(neurons),length(neurons))}; %empty variable to create a matrix the size of the whole connectome for each NPP_GPCR network
NPP_individual_networks_nerve_ring_connections_strata_3(:,1) = {zeros(length(neurons),length(neurons))}; %empty variable to create a matrix the size of the whole connectome for each NPP_GPCR network
NPP_individual_networks_nerve_ring_connections_strata_4(:,1) = {zeros(length(neurons),length(neurons))}; %empty variable to create a matrix the size of the whole connectome for each NPP_GPCR network
NPP_individual_networks_nerve_ring_connections_strata_5(:,1) = {zeros(length(neurons),length(neurons))}; %empty variable to create a matrix the size of the whole connectome for each NPP_GPCR network

for i = 1:length(NPP_individual_networks_nerve_ring_connections)
    NPP_individual_networks_nerve_ring_connections_strata_1{i,1}(find(NR_strata_class.Strata == 1), find(NR_strata_class.Strata == 1)) = NPP_individual_networks_nerve_ring_connections{i,1}(find(NR_strata_class.Strata == 1), find(NR_strata_class.Strata == 1)); % network for first strata for each NPP_GPCR pair
    NPP_individual_networks_nerve_ring_connections_strata_2{i,1}(find(NR_strata_class.Strata == 2), find(NR_strata_class.Strata == 2)) = NPP_individual_networks_nerve_ring_connections{i,1}(find(NR_strata_class.Strata == 2), find(NR_strata_class.Strata == 2)); % network for second strata for each NPP_GPCR pair
    NPP_individual_networks_nerve_ring_connections_strata_3{i,1}(find(NR_strata_class.Strata == 3), find(NR_strata_class.Strata == 3)) = NPP_individual_networks_nerve_ring_connections{i,1}(find(NR_strata_class.Strata == 3), find(NR_strata_class.Strata == 3)); % network for third strata for each NPP_GPCR pair
    NPP_individual_networks_nerve_ring_connections_strata_4{i,1}(find(NR_strata_class.Strata == 4), find(NR_strata_class.Strata == 4)) = NPP_individual_networks_nerve_ring_connections{i,1}(find(NR_strata_class.Strata == 4), find(NR_strata_class.Strata == 4)); % network for fourth strata for each NPP_GPCR pair 
    NPP_individual_networks_nerve_ring_connections_strata_5{i,1}(find(NR_strata_class.Strata == 5), find(NR_strata_class.Strata == 5)) = NPP_individual_networks_nerve_ring_connections{i,1}(find(NR_strata_class.Strata == 5), find(NR_strata_class.Strata == 5)); % network for fifth strata for each NPP_GPCR pair
end

[NPP_in_NR_1_connections_Amatrix, ~] = Connections_Amatrix (connections, NPP_individual_networks_nerve_ring_connections_strata_1, neurons, names_NPP,...
    'Number of connections in nerve ring strata 1 made by each of the individual NPP-GPCR pairs for each neuron node', [ path 'Comparison_btw_ranges'], date);
[NPP_in_NR_2_connections_Amatrix, ~] = Connections_Amatrix (connections, NPP_individual_networks_nerve_ring_connections_strata_2, neurons, names_NPP,...
    'Number of connections in nerve ring strata 2 made by each of the individual NPP-GPCR pairs for each neuron node', [ path 'Comparison_btw_ranges'], date);
[NPP_in_NR_3_connections_Amatrix, ~] = Connections_Amatrix (connections, NPP_individual_networks_nerve_ring_connections_strata_3, neurons, names_NPP,...
    'Number of connections in nerve ring strata 3 made by each of the individual NPP-GPCR pairs for each neuron node', [ path 'Comparison_btw_ranges'], date);
[NPP_in_NR_4_connections_Amatrix, ~] = Connections_Amatrix (connections, NPP_individual_networks_nerve_ring_connections_strata_4, neurons, names_NPP,...
    'Number of connections in nerve ring strata 4 made by each of the individual NPP-GPCR pairs for each neuron node', [ path 'Comparison_btw_ranges'], date);
[NPP_in_NR_5_connections_Amatrix, ~] = Connections_Amatrix (connections, NPP_individual_networks_nerve_ring_connections_strata_5, neurons, names_NPP,...
    'Number of connections in nerve ring neurons unassigned to strata and the rest of the nerve ring made by each of the individual NPP-GPCR pairs for each neuron node', [ path 'Comparison_btw_ranges'], date);

%Connections of NPP/GPCR pairs with 1 ligand and 1 specific receptor

[receptors_NPP_string, idx_con, idx_rec] = unique(connections(:, 'GeneGPCR'), 'rows'); % get list of receptors that form a NPP connection and get index to count how many ligands each receptor binds
NPP_individual_networks_1_ligand = sort(idx_con(find(accumarray(idx_rec,1) == 1))); % select rows from NPP_91_networks that have receptors that bind only 1 ligand
receptors_NPP_1_ligand = table2cell(connections(NPP_individual_networks_1_ligand, 'GeneGPCR')); %list of receptors with only 1 ligand
NPP_1_ligand = table2cell(connections(NPP_individual_networks_1_ligand, 'GeneLigand')); %list of NPP that bind to GPCR with only 1 ligand
receptors_NPP_GPCR_pairs_1_ligand = names_NPP(NPP_individual_networks_1_ligand, 1); %list of NPP_GPCR pairs with only 1 ligand
[~, idx_con1, idx_rec1] = unique(NPP_1_ligand, 'rows'); % get the neuropeptides that only bind to 1 receptor from the list of NPP that bind receptors that bind only 1 ligand. 
NPP_individual_networks_1_lig_1_rec = NPP_individual_networks_1_ligand(sort(idx_con1(find(accumarray(idx_rec1, 1) == 1)))); % select rows from NPP_91_networks that have receptors that bind only 1 ligand and 1 receptor
NPP_1_lig_1_rec = NPP_1_ligand(sort(idx_con1(find(accumarray(idx_rec1, 1) == 1)))); %get the neuropeptides in which the NPP only has one receptor and the receptor only has 1 ligand
receptors_NPP_GPCR_pairs_1_lig_1_rec = names_NPP(NPP_individual_networks_1_lig_1_rec, 1); %list of NPP_GPCR pairs with only 1 ligand
NPP_GPCR_connections_NR_in_btw_strata_Amatrix = NPP_GPCR_connections_NR_Amatrix - (NPP_in_NR_1_connections_Amatrix + NPP_in_NR_2_connections_Amatrix + NPP_in_NR_3_connections_Amatrix + NPP_in_NR_4_connections_Amatrix + NPP_in_NR_5_connections_Amatrix);

for i = 1:length(NPP_individual_networks_1_lig_1_rec) 
    j = NPP_individual_networks_1_lig_1_rec(i);
    [indeg_networks_1_lig_1_rec_short_range(i,:), ~, ~] = degrees_dir(NPP_individual_networks_short_range_connections{j}); %calculate indegree of each of the networks that uses a receptor that only has one ligand
    indeg_networks_1_lig_1_rec_short_range_NR_strata_1(i, :) = NPP_in_NR_1_connections_Amatrix((length(neurons)+1):(length(neurons)*2), j)'; %each row of the indegree matrix is the indegree inside strata 1 of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    indeg_networks_1_lig_1_rec_short_range_NR_strata_2(i, :) = NPP_in_NR_2_connections_Amatrix((length(neurons)+1):(length(neurons)*2), j)'; %each row of the indegree matrix is the indegree inside strata 2 of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    indeg_networks_1_lig_1_rec_short_range_NR_strata_3(i, :) = NPP_in_NR_3_connections_Amatrix((length(neurons)+1):(length(neurons)*2), j)'; %each row of the indegree matrix is the indegree inside strata 3 of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    indeg_networks_1_lig_1_rec_short_range_NR_strata_4(i, :) = NPP_in_NR_4_connections_Amatrix((length(neurons)+1):(length(neurons)*2), j)'; %each row of the indegree matrix is the indegree inside strata 4 of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    indeg_networks_1_lig_1_rec_short_range_NR_strata_5(i, :) = NPP_in_NR_5_connections_Amatrix((length(neurons)+1):(length(neurons)*2), j)'; %each row of the indegree matrix is the indegree in between unassigned neurons and other neurons of the NR of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    indeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata(i, :) = NPP_GPCR_connections_NR_in_btw_strata_Amatrix(303:604, j)'; %each row of the indegree matrix is the indegree in between strata of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    [indeg_networks_1_lig_1_rec_mid_range(i,:), ~, ~] = degrees_dir(NPP_individual_networks_mid_range_connections{j}); %calculate indegree of each of the networks that uses a receptor that only has one ligand
    [indeg_networks_1_lig_1_rec_long_range(i,:), ~, ~] = degrees_dir(neuropeptide_networks_all{j}); %calculate indegree of each of the networks that uses a receptor that only has one ligand
end


for i = 1:length(NPP_individual_networks_1_lig_1_rec) 
    j = NPP_individual_networks_1_lig_1_rec(i);
    [~,outdeg_networks_1_lig_1_rec_short_range(i,:), ~] = degrees_dir(NPP_individual_networks_short_range_connections{j}); %calculate outdegree of each of the networks that uses a receptor that only has one ligand
    outdeg_networks_1_lig_1_rec_short_range_NR_strata_1(i, :) = NPP_in_NR_1_connections_Amatrix(1:length(neurons), j)'; %each row of the outdegree matrix is the indegree inside strata 1 of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    outdeg_networks_1_lig_1_rec_short_range_NR_strata_2(i, :) = NPP_in_NR_2_connections_Amatrix(1:length(neurons), j)'; %each row of the outdegree matrix is the indegree inside strata 2 of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    outdeg_networks_1_lig_1_rec_short_range_NR_strata_3(i, :) = NPP_in_NR_3_connections_Amatrix(1:length(neurons), j)'; %each row of the outdegree matrix is the indegree inside strata 3 of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    outdeg_networks_1_lig_1_rec_short_range_NR_strata_4(i, :) = NPP_in_NR_4_connections_Amatrix(1:length(neurons), j)'; %each row of the outdegree matrix is the indegree inside strata 4 of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    outdeg_networks_1_lig_1_rec_short_range_NR_strata_5(i, :) = NPP_in_NR_5_connections_Amatrix(1:length(neurons), j)'; %each row of the outdegree matrix is the indegree in between unassigned neurons and other neurons of the NR of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    outdeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata(i, :) = NPP_GPCR_connections_NR_in_btw_strata_Amatrix(1:length(neurons), j)'; %each row of the outdegree matrix is the indegree in between strata of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    [~,outdeg_networks_1_lig_1_rec_mid_range(i,:), ~] = degrees_dir(NPP_individual_networks_mid_range_connections{j}); %calculate outdegree of each of the networks that uses a receptor that only has one ligand
    [~,outdeg_networks_1_lig_1_rec_long_range(i,:), ~] = degrees_dir(neuropeptide_networks_all{j}); %calculate outdegree of each of the networks that uses a receptor that only has one ligand
end


%Plot overlap matrix of receptors that bind 1 ligand and their cognate NPP
%that only bind that receptor

figure('Renderer', 'painters', 'Visible', 'off', 'PaperOrientation', 'landscape')
tiledlayout(2, 1, 'TileSpacing','compact', 'Padding', 'tight');
nexttile
[y_indeg_networks_1_lig_1_rec_long_range_sort, x_indeg_networks_1_lig_1_rec_long_range_sort ] = find(indeg_networks_1_lig_1_rec_long_range); %get x and y coordinates to plot scatter
scatter(x_indeg_networks_1_lig_1_rec_long_range_sort, y_indeg_networks_1_lig_1_rec_long_range_sort, 15, 'MarkerFaceColor', 1/255*[211 211 211], 'MarkerEdgeColor', 'none'); %plot scatter of long range (scatter instead of spy to choose colour)
hold on 
[y_indeg_networks_1_lig_1_rec_mid_range_sort, x_indeg_networks_1_lig_1_rec_mid_range_sort ] = find(indeg_networks_1_lig_1_rec_mid_range); %get x and y coordinates to plot scatter
scatter(x_indeg_networks_1_lig_1_rec_mid_range_sort, y_indeg_networks_1_lig_1_rec_mid_range_sort, 15, 'MarkerFaceColor', 1/255*[215 25 28], 'MarkerEdgeColor', 'none'); %plot scatter of mid range (scatter instead of spy to choose colour)
[y_indeg_networks_1_lig_1_rec_short_range_sort, x_indeg_networks_1_lig_1_rec_short_range_sort ] = find(indeg_networks_1_lig_1_rec_short_range); %get x and y coordinates to plot scatter
scatter(x_indeg_networks_1_lig_1_rec_short_range_sort, y_indeg_networks_1_lig_1_rec_short_range_sort, 15, 'MarkerFaceColor', 1/255*[171 217 233], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
[y_indeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata_sort, x_indeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata_sort ] = find(indeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata); %get x and y coordinates to plot scatter
scatter(x_indeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata_sort, y_indeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata_sort, 15, 'MarkerFaceColor', 1/255*[171 221 164], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
[y_indeg_networks_1_lig_1_rec_short_range_NR_strata_1_sort, x_indeg_networks_1_lig_1_rec_short_range_NR_strata_1_sort ] = find(indeg_networks_1_lig_1_rec_short_range_NR_strata_1); %get x and y coordinates to plot scatter
scatter(x_indeg_networks_1_lig_1_rec_short_range_NR_strata_1_sort, y_indeg_networks_1_lig_1_rec_short_range_NR_strata_1_sort, 15, 'MarkerFaceColor', 1/255*[43 131 186], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
[y_indeg_networks_1_lig_1_rec_short_range_NR_strata_2_sort, x_indeg_networks_1_lig_1_rec_short_range_NR_strata_2_sort ] = find(indeg_networks_1_lig_1_rec_short_range_NR_strata_2); %get x and y coordinates to plot scatter
scatter(x_indeg_networks_1_lig_1_rec_short_range_NR_strata_2_sort, y_indeg_networks_1_lig_1_rec_short_range_NR_strata_2_sort, 15, 'MarkerFaceColor', 1/255*[43 131 186], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
[y_indeg_networks_1_lig_1_rec_short_range_NR_strata_3_sort, x_indeg_networks_1_lig_1_rec_short_range_NR_strata_3_sort ] = find(indeg_networks_1_lig_1_rec_short_range_NR_strata_3); %get x and y coordinates to plot scatter
scatter(x_indeg_networks_1_lig_1_rec_short_range_NR_strata_3_sort, y_indeg_networks_1_lig_1_rec_short_range_NR_strata_3_sort, 15, 'MarkerFaceColor', 1/255*[43 131 186], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
[y_indeg_networks_1_lig_1_rec_short_range_NR_strata_4_sort, x_indeg_networks_1_lig_1_rec_short_range_NR_strata_4_sort ] = find(indeg_networks_1_lig_1_rec_short_range_NR_strata_4); %get x and y coordinates to plot scatter
scatter(x_indeg_networks_1_lig_1_rec_short_range_NR_strata_4_sort, y_indeg_networks_1_lig_1_rec_short_range_NR_strata_4_sort, 15, 'MarkerFaceColor', 1/255*[43 131 186], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
[y_indeg_networks_1_lig_1_rec_short_range_NR_strata_5_sort, x_indeg_networks_1_lig_1_rec_short_range_NR_strata_5_sort ] = find(indeg_networks_1_lig_1_rec_short_range_NR_strata_5); %get x and y coordinates to plot scatter
scatter(x_indeg_networks_1_lig_1_rec_short_range_NR_strata_5_sort, y_indeg_networks_1_lig_1_rec_short_range_NR_strata_5_sort, 15, 'MarkerFaceColor', 1/255*[43 131 186], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
xlim([0 size(indeg_networks_1_lig_1_rec_long_range, 2)+1])
ylim([0 size(indeg_networks_1_lig_1_rec_long_range, 1)+1])
hold off
title('Diffusion distance comparison based on GPCR expression', 'interpreter', 'none', 'FontSize', 5); %label the title of the plot
set(gca, 'TickLength', [0, 0], 'YTickLabel', [], 'XTickLabel', [], 'YColor', 'k', 'XColor', 'k'); %remove tick marks and prevent unwanted underscores in tick labels
box on 
nexttile
[y_outdeg_networks_1_lig_1_rec_long_range_sort, x_outdeg_networks_1_lig_1_rec_long_range_sort ] = find(outdeg_networks_1_lig_1_rec_long_range); %get x and y coordinates to plot scatter
scatter(x_outdeg_networks_1_lig_1_rec_long_range_sort, y_outdeg_networks_1_lig_1_rec_long_range_sort, 15, 'MarkerFaceColor', 1/255*[211 211 211], 'MarkerEdgeColor', 'none'); %plot scatter of long range (scatter instead of spy to choose colour)
hold on 
[y_outdeg_networks_1_lig_1_rec_mid_range_sort, x_outdeg_networks_1_lig_1_rec_mid_range_sort ] = find(outdeg_networks_1_lig_1_rec_mid_range); %get x and y coordinates to plot scatter
s2 = scatter(x_outdeg_networks_1_lig_1_rec_mid_range_sort, y_outdeg_networks_1_lig_1_rec_mid_range_sort, 15, 'MarkerFaceColor', 1/255*[215 25 28], 'MarkerEdgeColor', 'none'); %plot scatter of mid range (scatter instead of spy to choose colour)
[y_outdeg_networks_1_lig_1_rec_short_range_sort, x_outdeg_networks_1_lig_1_rec_short_range_sort ] = find(outdeg_networks_1_lig_1_rec_short_range); %get x and y coordinates to plot scatter
s3 = scatter(x_outdeg_networks_1_lig_1_rec_short_range_sort, y_outdeg_networks_1_lig_1_rec_short_range_sort, 15, 'MarkerFaceColor', 1/255*[171 217 233], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
[y_outdeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata_sort, x_outdeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata_sort ] = find(outdeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata); %get x and y coordinates to plot scatter
s9 = scatter(x_outdeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata_sort, y_outdeg_networks_1_lig_1_rec_short_range_NR_in_btw_strata_sort, 15, 'MarkerFaceColor', 1/255*[171 221 164], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
[y_outdeg_networks_1_lig_1_rec_short_range_NR_strata_1_sort, x_outdeg_networks_1_lig_1_rec_short_range_NR_strata_1_sort ] = find(outdeg_networks_1_lig_1_rec_short_range_NR_strata_1); %get x and y coordinates to plot scatter
s4 = scatter(x_outdeg_networks_1_lig_1_rec_short_range_NR_strata_1_sort, y_outdeg_networks_1_lig_1_rec_short_range_NR_strata_1_sort, 15, 'MarkerFaceColor', 1/255*[43 131 186], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
[y_outdeg_networks_1_lig_1_rec_short_range_NR_strata_2_sort, x_outdeg_networks_1_lig_1_rec_short_range_NR_strata_2_sort ] = find(outdeg_networks_1_lig_1_rec_short_range_NR_strata_2); %get x and y coordinates to plot scatter
scatter(x_outdeg_networks_1_lig_1_rec_short_range_NR_strata_2_sort, y_outdeg_networks_1_lig_1_rec_short_range_NR_strata_2_sort, 15, 'MarkerFaceColor', 1/255*[43 131 186], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
[y_outdeg_networks_1_lig_1_rec_short_range_NR_strata_3_sort, x_outdeg_networks_1_lig_1_rec_short_range_NR_strata_3_sort ] = find(outdeg_networks_1_lig_1_rec_short_range_NR_strata_3); %get x and y coordinates to plot scatter
scatter(x_outdeg_networks_1_lig_1_rec_short_range_NR_strata_3_sort, y_outdeg_networks_1_lig_1_rec_short_range_NR_strata_3_sort, 15, 'MarkerFaceColor', 1/255*[43 131 186], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
[y_outdeg_networks_1_lig_1_rec_short_range_NR_strata_4_sort, x_outdeg_networks_1_lig_1_rec_short_range_NR_strata_4_sort ] = find(outdeg_networks_1_lig_1_rec_short_range_NR_strata_4); %get x and y coordinates to plot scatter
scatter(x_outdeg_networks_1_lig_1_rec_short_range_NR_strata_4_sort, y_outdeg_networks_1_lig_1_rec_short_range_NR_strata_4_sort, 15, 'MarkerFaceColor', 1/255*[43 131 186], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
[y_outdeg_networks_1_lig_1_rec_short_range_NR_strata_5_sort, x_outdeg_networks_1_lig_1_rec_short_range_NR_strata_5_sort ] = find(outdeg_networks_1_lig_1_rec_short_range_NR_strata_5); %get x and y coordinates to plot scatter
scatter(x_outdeg_networks_1_lig_1_rec_short_range_NR_strata_5_sort, y_outdeg_networks_1_lig_1_rec_short_range_NR_strata_5_sort, 15, 'MarkerFaceColor', 1/255*[43 131 186], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
xlim([0 size(indeg_networks_1_lig_1_rec_long_range, 2)+1])
ylim([0 size(indeg_networks_1_lig_1_rec_long_range, 1)+1])
hold off
title('Diffusion distance comparison based on NPP expression', 'interpreter', 'none', 'FontSize', 5); %label the title of the plot
box on 
set(gca, 'TickLength', [0, 0], 'YTickLabel', [], 'XTickLabel', [], 'YColor', 'k', 'XColor', 'k'); %remove tick marks and prevent unwanted underscores in tick labels
legend([s2 s9 s3 s4], 'Requires mid range', 'Requires in between strata', 'Whithin bundle (outside NR)', 'Whithin stratum', 'FontSize', 5, 'Location', 'southoutside', 'Orientation', 'horizontal')
print(['Results/scRNAseq_C.elegans/' path 'Comparison_btw_ranges/' 'Scatter_model_comparison_exp'], '-dpdf', '-bestfit', '-r300');%save image as pdf

%Plot overlap matrix of receptors that bind 1 ligand 
for i = 1:length(NPP_individual_networks_1_ligand) 
    j = NPP_individual_networks_1_ligand(i);
    [indeg_networks_1_lig_short_range(i,:), ~, ~] = degrees_dir(NPP_individual_networks_short_range_connections{j}); %calculate indegree of each of the networks that uses a receptor that only has one ligand
    indeg_networks_1_lig_short_range_NR_strata_1(i, :) = NPP_in_NR_1_connections_Amatrix(303:604, j)'; %each row of the indegree matrix is the indegree inside strata 1 of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    indeg_networks_1_lig_short_range_NR_strata_2(i, :) = NPP_in_NR_2_connections_Amatrix(303:604, j)'; %each row of the indegree matrix is the indegree inside strata 2 of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    indeg_networks_1_lig_short_range_NR_strata_3(i, :) = NPP_in_NR_3_connections_Amatrix(303:604, j)'; %each row of the indegree matrix is the indegree inside strata 3 of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    indeg_networks_1_lig_short_range_NR_strata_4(i, :) = NPP_in_NR_4_connections_Amatrix(303:604, j)'; %each row of the indegree matrix is the indegree inside strata 4 of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    indeg_networks_1_lig_short_range_NR_strata_5(i, :) = NPP_in_NR_5_connections_Amatrix(303:604, j)'; %each row of the indegree matrix is the indegree in between unassigned neurons and other neurons of the NR of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    indeg_networks_1_lig_short_range_NR_in_btw_strata(i, :) = NPP_GPCR_connections_NR_in_btw_strata_Amatrix(303:604, j)'; %each row of the indegree matrix is the indegree in between strata of the nerve ring for the specific network of one of the receptors that only has 1 ligand
    [indeg_networks_1_lig_mid_range(i,:), ~, ~] = degrees_dir(NPP_individual_networks_mid_range_connections{j}); %calculate indegree of each of the networks that uses a receptor that only has one ligand
    [indeg_networks_1_lig_long_range(i,:), ~, ~] = degrees_dir(neuropeptide_networks_all{j}); %calculate indegree of each of the networks that uses a receptor that only has one ligand
end

figure('Renderer', 'painters', 'Visible', 'off', 'PaperOrientation', 'landscape')
tiledlayout(2, 1, 'TileSpacing','compact', 'Padding', 'tight');
nexttile
[y_indeg_networks_1_lig_long_range_sort, x_indeg_networks_1_lig_long_range_sort ] = find(indeg_networks_1_lig_long_range); %get x and y coordinates to plot scatter
scatter(x_indeg_networks_1_lig_long_range_sort, y_indeg_networks_1_lig_long_range_sort, 15, 'MarkerFaceColor', 1/255*[211 211 211], 'MarkerEdgeColor', 'none'); %plot scatter of long range (scatter instead of spy to choose colour)
hold on 
[y_indeg_networks_1_lig_mid_range_sort, x_indeg_networks_1_lig_mid_range_sort ] = find(indeg_networks_1_lig_mid_range); %get x and y coordinates to plot scatter
s2 = scatter(x_indeg_networks_1_lig_mid_range_sort, y_indeg_networks_1_lig_mid_range_sort, 15, 'MarkerFaceColor', 1/255*[215 25 28], 'MarkerEdgeColor', 'none'); %plot scatter of mid range (scatter instead of spy to choose colour)
[y_indeg_networks_1_lig_short_range_sort, x_indeg_networks_1_lig_short_range_sort ] = find(indeg_networks_1_lig_short_range); %get x and y coordinates to plot scatter
s3 = scatter(x_indeg_networks_1_lig_short_range_sort, y_indeg_networks_1_lig_short_range_sort, 15, 'MarkerFaceColor', 1/255*[171 217 233], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
[y_indeg_networks_1_lig_short_range_NR_in_btw_strata_sort, x_indeg_networks_1_lig_short_range_NR_in_btw_strata_sort ] = find(indeg_networks_1_lig_short_range_NR_in_btw_strata); %get x and y coordinates to plot scatter
s9 = scatter(x_indeg_networks_1_lig_short_range_NR_in_btw_strata_sort, y_indeg_networks_1_lig_short_range_NR_in_btw_strata_sort, 15, 'MarkerFaceColor', 1/255*[255 243 13], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
[y_indeg_networks_1_lig_short_range_NR_strata_1_sort, x_indeg_networks_1_lig_short_range_NR_strata_1_sort ] = find(indeg_networks_1_lig_short_range_NR_strata_1); %get x and y coordinates to plot scatter
s4 = scatter(x_indeg_networks_1_lig_short_range_NR_strata_1_sort, y_indeg_networks_1_lig_short_range_NR_strata_1_sort, 15, 'MarkerFaceColor', 1/255*[43 131 186], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
[y_indeg_networks_1_lig_short_range_NR_strata_2_sort, x_indeg_networks_1_lig_short_range_NR_strata_2_sort ] = find(indeg_networks_1_lig_short_range_NR_strata_2); %get x and y coordinates to plot scatter
scatter(x_indeg_networks_1_lig_short_range_NR_strata_2_sort, y_indeg_networks_1_lig_short_range_NR_strata_2_sort, 15, 'MarkerFaceColor', 1/255*[43 131 186], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
[y_indeg_networks_1_lig_short_range_NR_strata_3_sort, x_indeg_networks_1_lig_short_range_NR_strata_3_sort ] = find(indeg_networks_1_lig_short_range_NR_strata_3); %get x and y coordinates to plot scatter
scatter(x_indeg_networks_1_lig_short_range_NR_strata_3_sort, y_indeg_networks_1_lig_short_range_NR_strata_3_sort, 15, 'MarkerFaceColor', 1/255*[43 131 186], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
[y_indeg_networks_1_lig_short_range_NR_strata_4_sort, x_indeg_networks_1_lig_short_range_NR_strata_4_sort ] = find(indeg_networks_1_lig_short_range_NR_strata_4); %get x and y coordinates to plot scatter
scatter(x_indeg_networks_1_lig_short_range_NR_strata_4_sort, y_indeg_networks_1_lig_short_range_NR_strata_4_sort, 15, 'MarkerFaceColor', 1/255*[43 131 186], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
[y_indeg_networks_1_lig_short_range_NR_strata_5_sort, x_indeg_networks_1_lig_short_range_NR_strata_5_sort ] = find(indeg_networks_1_lig_short_range_NR_strata_5); %get x and y coordinates to plot scatter
scatter(x_indeg_networks_1_lig_short_range_NR_strata_5_sort, y_indeg_networks_1_lig_short_range_NR_strata_5_sort, 15, 'MarkerFaceColor', 1/255*[43 131 186], 'MarkerEdgeColor', 'none'); %plot scatter of short range (scatter instead of spy to choose colour)
xlim([0 size(indeg_networks_1_lig_1_rec_long_range, 2)+1])
ylim([0 size(indeg_networks_1_lig_1_rec_long_range, 1)+1])
hold off
title('Diffusion distance comparison based on GPCR expression 1 ligand', 'interpreter', 'none', 'FontSize', 5); %label the title of the plot
box on 
set(gca, 'TickLength', [0, 0], 'YTickLabel', [], 'XTickLabel', [], 'YColor', 'k', 'XColor', 'k'); %remove tick marks and prevent unwanted underscores in tick labels
legend([s2 s9 s3 s4], 'Requires mid range', 'Requires in between strata', 'Whithin bundle (outside NR)', 'Whithin stratum', 'FontSize', 5, 'Location', 'southoutside', 'Orientation', 'horizontal')
print(['Results/scRNAseq_C.elegans/' path 'Comparison_btw_ranges/' 'Scatter_model_comparison_GPCR_1lig'], '-dpdf', '-bestfit', '-r300');%save image as pdf

%**************************************************************************
%SCATTER PLOT NPP/GPCR CONNECTIONS BASED ON ANATOMY
%**************************************************************************

%Plot all NPP/GPCR but colour by anatomical mid range area 
for i = 1:length(NPP_individual_networks_short_range_connections) %loop over every neuron 
    [indegree_NPP_individual_short_range(i,:), ~,~] = degrees_dir(NPP_individual_networks_short_range_connections{i}); %indegree of short range connections by NPP_GPCR pairs
    [ ~, outdegree_NPP_individual_short_range(i,:), ~] = degrees_dir(NPP_individual_networks_short_range_connections{i}); %outdegree of short range connections by NPP_GPCR pairs
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

image15 = figure('Visible', 'off', 'Resize', 'off', 'Position', [199.625,350,3600,4408]); %create an image to save the figure with the right size and prevent it being resize by the plots
t = tiledlayout(6, 1,'TileSpacing','compact', 'Padding', 'tight');
nexttile 
s1 = scatter(x_indeg_NPP_individual_short_range, y_indeg_NPP_individual_short_range, 3, colors_neurons_Head_ind, 'filled'); %plot scatter of short range (scatter instead of spy to choose colour)
yticks( linspace(1, length(names_NPP), length(names_NPP))); %add all y axis ticks
yticklabels(names_NPP); %labels for y axes
set(gca, 'TickLength', [0, 0], 'TickLabelInterpreter','none', 'XTickLabel', [], 'YColor', 'k', 'XColor', 'k', 'FontSize', 1); %remove tick marks and prevent unwanted underscores in tick labels
title('GPCR expression by neuron with anatomical resolution', FontSize=5)
xlim([0 size(indegree_NPP_individual_short_range, 2)+1])
ylim([0 size(indegree_NPP_individual_short_range, 1)+1])
xline(20);
xline(99);
xline(181);
box on 
nexttile
s2 = scatter(x_indeg_NPP_individual_short_range, y_indeg_NPP_individual_short_range, 3, colors_neurons_Midbody_ind, 'filled'); %plot scatter of short range (scatter instead of spy to choose colour)
yticks( linspace(1, length(names_NPP), length(names_NPP))); %add all y axis ticks
yticklabels(names_NPP); %labels for y axes
set(gca, 'TickLength', [0, 0], 'TickLabelInterpreter','none', 'XTickLabel', [], 'YColor', 'k', 'XColor', 'k', 'FontSize', 1); %remove tick marks and prevent unwanted underscores in tick labels
xlim([0 size(indegree_NPP_individual_short_range, 2)+1])
ylim([0 size(indegree_NPP_individual_short_range, 1)+1])
xline(20);
xline(99);
xline(181);
box on 
nexttile
s3 = scatter(x_indeg_NPP_individual_short_range, y_indeg_NPP_individual_short_range, 3, colors_neurons_Tail_ind, 'filled'); %plot scatter of short range (scatter instead of spy to choose colour)
yticks( linspace(1, length(names_NPP), length(names_NPP))); %add all y axis ticks
yticklabels(names_NPP); %labels for y axes
set(gca, 'TickLength', [0, 0], 'TickLabelInterpreter','none', 'XTickLabel', [], 'YColor', 'k', 'XColor', 'k', 'FontSize', 1); %remove tick marks and prevent unwanted underscores in tick labels
xlim([0 size(indegree_NPP_individual_short_range, 2)+1])
ylim([0 size(indegree_NPP_individual_short_range, 1)+1])
xline(20);
xline(99);
xline(181);
box on 
nexttile 
s1 = scatter(x_outdeg_NPP_individual_short_range, y_outdeg_NPP_individual_short_range, 3, colors_neurons_Head_out, 'filled'); %plot scatter of short range (scatter instead of spy to choose colour)
yticks( linspace(1, length(names_NPP), length(names_NPP))); %add all y axis ticks
yticklabels(names_NPP); %labels for y axes
set(gca, 'TickLength', [0, 0], 'TickLabelInterpreter','none', 'XTickLabel', [], 'YColor', 'k', 'XColor', 'k', 'FontSize', 1); %remove tick marks and prevent unwanted underscores in tick labels
title('NPP expression by neuron with anatomical resolution', FontSize=5)
xlim([0 size(outdegree_NPP_individual_short_range, 2)+1])
ylim([0 size(outdegree_NPP_individual_short_range, 1)+1])
xline(20);
xline(99);
xline(181);
box on 
nexttile
s2 = scatter(x_outdeg_NPP_individual_short_range, y_outdeg_NPP_individual_short_range, 3, colors_neurons_Midbody_out, 'filled'); %plot scatter of short range (scatter instead of spy to choose colour)
yticks( linspace(1, length(names_NPP), length(names_NPP))); %add all y axis ticks
yticklabels(names_NPP); %labels for y axes
set(gca, 'TickLength', [0, 0], 'TickLabelInterpreter','none', 'XTickLabel', [], 'YColor', 'k', 'XColor', 'k', 'FontSize', 1); %remove tick marks and prevent unwanted underscores in tick labels
xlim([0 size(outdegree_NPP_individual_short_range, 2)+1])
ylim([0 size(outdegree_NPP_individual_short_range, 1)+1])
xline(20);
xline(99);
xline(181);
box on 
nexttile
s3 = scatter(x_outdeg_NPP_individual_short_range, y_outdeg_NPP_individual_short_range, 3, colors_neurons_Tail_out, 'filled'); %plot scatter of short range (scatter instead of spy to choose colour)
yticks( linspace(1, length(names_NPP), length(names_NPP))); %add all y axis ticks
yticklabels(names_NPP); %labels for y axes
ylabel(t, 'NPP-GPCR expression') %y axes label
set(gca, 'TickLength', [0, 0], 'TickLabelInterpreter','none', 'XTickLabel', [], 'YColor', 'k', 'XColor', 'k', 'FontSize', 1); %remove tick marks and prevent unwanted underscores in tick labels
xlim([0 size(outdegree_NPP_individual_short_range, 2)+1])
ylim([0 size(outdegree_NPP_individual_short_range, 1)+1])
xline(20);
xline(99);
xline(181);
box on 
leg = legend([s1 s2 s3], 'Process bundles in the Head', 'Process bundles in Midbody', 'Process bundles in Tail', 'FontSize', 5, 'Orientation', 'horizontal');
leg.Layout.Tile = 'south';
print(['Results/scRNAseq_C.elegans/' path 'Comparison_btw_ranges/' 'Scatter_exp_NPP_GPCR_anatomy'], '-dpdf', '-fillpage', '-r300');%save image as pdf

%**************************************************************************
%SAVE
%**************************************************************************
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
    row_names = cat(1, neuron_names, neuron_names); %names of matrix rows
else %the input is a matrix not a cell array, output is a column vector
    vector = zeros(2*length(data_networks), 1); %store the data temporary in iteration
        for j = 1:length(data_networks)
            vector(j,1) = sum(data_networks(j, :));%add all the edges in a row for the matrix of the NPP_GPCR pair
            vector(j+length(data_networks), 1) = sum(data_networks(:,j)); %add all the edges for the same neuron in the other dimension of the matrix of the NPP_GPCR pair
        end
     Amatrix = vector; %end result
     row_names = cat(1, neuron_names, neuron_names); %names of matrix rows
end 

end
%**************************************************************************
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
function Co_expression_analysis(connections_NPP, path, path2, date, version)

%Load data
load([pwd '/Results/scRNAseq_C.elegans' path2 date 'Networks_build_EC500nM.mat'], 'expression_byneuron_ordered', 'geneID');

expression_neuron_classes = readtable ('29072020_threshold4_CENGEN_data_TPM_expression_final_release'); %import data in neuron class format
expression_neuron_classes = join(geneID, expression_neuron_classes);  %add the gene name that corresponds to each WBidentifier
expression_neuron_classes = expression_neuron_classes(:,[1 3:129]); %remove the WBidentifier column

neuropeptides_unsorted = unique(connections_NPP.GeneLigand, 'stable'); %get neuropeptide names
receptors_NPP_unsorted = unique(connections_NPP.GeneGPCR, 'stable'); %get neuropeptide names

mkdir ([ 'Results/scRNAseq_C.elegans' path 'Co-expression_analysis']) %create folder to save images

%Build co-expression map and correct
[Fisher_logical_values_005_byclass, Fisher_pvalues_005_byclass, Fisher_stats_005_byclass, Fisher_logical_values_001_byclass, ...
    Fisher_pvalues_001_byclass, Fisher_stats_001_byclass, Fisher_Contingency_table_NPP_GPCR_byclass, Fisher_pFDR_005_byclass, ...
    Fisher_pFDR_001_byclass, Fisher_logical_values_005_corrected_byclass, Fisher_logical_values_001_corrected_byclass] = Fishertestwithtable(expression_neuron_classes, ...
    neuropeptides_unsorted, receptors_NPP_unsorted, 'NPP_to_GPCR_all_genes', 'NPP expression', 'GPCR expression'); %use the Fishertestwithtable function 
[Fisher_logical_values_005_byneuron, Fisher_pvalues_005_byneuron, Fisher_stats_005_byneuron, Fisher_logical_values_001_byneuron,...
    Fisher_pvalues_001_byneuron, Fisher_stats_001_byneuron, Fisher_Contingency_table_NPP_GPCR_byneuron, Fisher_pFDR_005_byneuron,...
    Fisher_pFDR_001_byneuron, Fisher_logical_values_005_corrected_byneuron, Fisher_logical_values_001_corrected_byneuron] = Fishertestwithtable(expression_byneuron_ordered,...
    neuropeptides_unsorted, receptors_NPP_unsorted, 'NPP_to_GPCR_all_genes_302_neurons', 'NPP expression', 'GPCR expression'); %use the Fishertestwithtable function 

%calculate the FDR for the byclass co-expression that overlaps with the by
%neuron FDR corrected co-expression instances
for i = 1:height(Fisher_logical_values_005_byclass) %loop over every co-expression instance
    for j = 1:length(Fisher_logical_values_005_byclass)%loop over every co-expression instance
    if Fisher_logical_values_005_byclass(i,j) > 0 && Fisher_logical_values_005_corrected_byneuron(i,j) > 0 %select only the ones that overlap btw both distributions
        Fisher_pvalues_005_byclass_overlap(i,j) = Fisher_pvalues_005_byclass(i,j); %give those the p-value of by class
    else 
        Fisher_pvalues_005_byclass_overlap(i,j) = 1; %others 1 because they are not significant
    end
    end
end
Fisher_pvalues_005_byclass_overlap_r = reshape(Fisher_pvalues_005_byclass_overlap.',1,[]); % transform matrix into vector to run FDR on all values
Fisher_pvalues_005_byclass_overlap_rr = Fisher_pvalues_005_byclass_overlap_r; %create copy of matrix
Fisher_pvalues_005_byclass_overlap_rr(Fisher_pvalues_005_byclass_overlap_rr==1) = []; %remove instances where the values do not overlap so that they do not spoil FDR making significant p-value smaller
[Fisher_pvalues_005_byclass_overlap_corr] = mafdr(Fisher_pvalues_005_byclass_overlap_rr); %calculate FDR over only overlap values
Fisher_pvalues_005_byclass_overlap_e = ones(1, height(neuropeptides_unsorted)*height(receptors_NPP_unsorted)); % empty matrix to give the vector again the origianl size
Fisher_pvalues_005_byclass_overlap_e(Fisher_pvalues_005_byclass_overlap_r~=1) = Fisher_pvalues_005_byclass_overlap_corr;
Fisher_logical_values_005_byclass_overlap_FDR005 = reshape(Fisher_pvalues_005_byclass_overlap_e <= 0.05, [height(receptors_NPP_unsorted), height(neuropeptides_unsorted)]).'; %reshape and threshold
Fisher_pvalues_005_byclass_overlap_FDR005 = reshape(Fisher_pvalues_005_byclass_overlap_e, [height(receptors_NPP_unsorted), height(neuropeptides_unsorted)]).'; %reshape

%Build networks based on deorphanisation data using the co-expression
%significant instances as nodes
plot_coexpression_data(Fisher_pvalues_005_byclass_overlap_FDR005, Fisher_logical_values_005_byclass_overlap_FDR005, ...
    neuropeptides_unsorted, receptors_NPP_unsorted, [], [], [path 'Co-expression_analysis'], 'NPP-GPCR Comparison');

[Fisher_005_byclasssoverlap_corrected_network, Fisher_005_byclassoverlap_corrected_network_unsort, Fisher_005_byclasssoverlap_corrected_network_x, Fisher_005_byclasssoverlap_corrected_network_y, Fisher_005_byclasssoverlap_corrected_networks_name_GPCR_t, Fisher_005_byclassoverlap_corrected_networks_name_NPP_t] = build_coexpression_network(Fisher_logical_values_005_byclass_overlap_FDR005, neuropeptides_unsorted, receptors_NPP_unsorted, connections_NPP); %run for all coexpression FDR corrected by class
plot_coexpression_network (Fisher_005_byclasssoverlap_corrected_network, [date, 'Co-occurrence network for 5% significance by class overlap'], [path 'Co-expression_analysis'], Fisher_005_byclasssoverlap_corrected_network_x, 0); %plot network

%For comparison between NPP vs NPP and GPCR vs GPCR
[Fisher_logical_values_005_byneuron_NPP2, Fisher_pvalues_005_byneuron_NPP2, Fisher_stats_005_byneuron_NPP2, Fisher_logical_values_001_byneuron_NPP2,...
    Fisher_pvalues_001_byneuron_NPP2, Fisher_stats_001_byneuron_NPP2, Fisher_Contingency_table_NPP_GPCR_byneuron_NPP2, Fisher_pFDR_005_byneuron_NPP2,...
    Fisher_pFDR_001_byneuron_NPP2, Fisher_logical_values_005_corrected_byneuron_NPP2, Fisher_logical_values_001_corrected_byneuron_NPP2] = Fishertestwithtable(expression_byneuron_ordered,...
    neuropeptides_unsorted, neuropeptides_unsorted, 'NPP_to_NPP_all_genes_302_neurons', 'NPP expression', 'NPP expression'); %use the Fishertestwithtable function 
[Fisher_logical_values_005_byneuron_GPCR2, Fisher_pvalues_005_byneuron_GPCR2, Fisher_stats_005_byneuron_GPCR2, Fisher_logical_values_001_byneuron_GPCR2,...
    Fisher_pvalues_001_byneuron_GPCR2, Fisher_stats_001_byneuron_GPCR2, Fisher_Contingency_table_NPP_GPCR_byneuron_GPCR2, Fisher_pFDR_005_byneuron_GPCR2,...
    Fisher_pFDR_001_byneuron_GPCR2, Fisher_logical_values_005_corrected_byneuron_GPCR2, Fisher_logical_values_001_corrected_byneuron_GPCR2] = Fishertestwithtable(expression_byneuron_ordered,...
    receptors_NPP_unsorted, receptors_NPP_unsorted, 'GPCR_to_GPCR_all_genes_302_neurons', 'GPCR expression', 'GPCR expression'); %use the Fishertestwithtable function 

plot_coexpression_data(Fisher_pFDR_005_byneuron_NPP2, Fisher_logical_values_005_corrected_byneuron_NPP2, ...
    neuropeptides_unsorted, receptors_NPP_unsorted, connections_NPP.GeneGPCR, connections_NPP.GeneLigand, [path 'Co-expression_analysis'], 'Neuropeptides');

plot_coexpression_data(Fisher_pFDR_005_byneuron_GPCR2, Fisher_logical_values_005_corrected_byneuron_GPCR2, ...
    receptors_NPP_unsorted, neuropeptides_unsorted, connections_NPP.GeneLigand, connections_NPP.GeneGPCR, [path 'Co-expression_analysis'], 'Receptors');

%Save data
save((['Results/scRNAseq_C.elegans' path2 date 'Networks_Neuropeptide_coexpression_T4_EC500nM.mat']))

end

function [network_sorted, network_unsorted, network_x, network_y, data_network_GPCR_t, data_network_NPP_t] = build_coexpression_network(data, neuropeptides, receptors, connection)
   %build the list of nodes for the x and y dimensions of the network (y NPP
   %and x GPCR in the co-occuring pair) 
    for i = 1:height(data)
        for j= 1:width(data)
            if data(i,j) == 1 
                data_network_NPP_y{j,i} = neuropeptides{i}; %get neuropeptide name 
                data_network_NPP_x{j,i} = receptors{j}; %take receptor names that bind to the NPP
            end
        end
        data_network_NPP_t = table(data_network_NPP_y(~cellfun('isempty', data_network_NPP_y)),...
            data_network_NPP_x(~cellfun('isempty', data_network_NPP_x)), 'VariableNames', {'NPP', 'GPCR'}); %remove empty cells
        data_network_GPCR_y = data_network_NPP_y.';
        data_network_GPCR_x = data_network_NPP_x.';
        data_network_GPCR_t = table(data_network_GPCR_y(~cellfun('isempty', data_network_GPCR_y)),...
            data_network_GPCR_x(~cellfun('isempty', data_network_GPCR_x)), 'VariableNames', {'NPP', 'GPCR'}); %remove empty cells
    end
    
    %build network
    for i = 1: height(data_network_GPCR_t)
        data_network_pairing{i} = connection.GeneLigand(find(matches(connection.GeneGPCR, data_network_GPCR_t.GPCR(i)))); %find the NPP that binds the GPCR in the first node
        network_sorted(:,i) = double(matches(data_network_NPP_t.NPP, data_network_pairing{i})); %logic array of NPP in vector of co-occurrence
        network_unsorted(:,i) = double(matches(data_network_GPCR_t.NPP, data_network_pairing{i})); %logic array of NPP in vector of co-occurrence
        network_y = append (matlab.lang.makeValidName(data_network_NPP_t.NPP), '_', matlab.lang.makeValidName(data_network_NPP_t.GPCR)); %name for every node in the network in y axes, so every co-expression
        network_x = append (matlab.lang.makeValidName(data_network_GPCR_t.NPP), '_', matlab.lang.makeValidName(data_network_GPCR_t.GPCR)); %name for every node in the network in x axes, so every co-expression
    end
end
function plot_coexpression_network (data, names_list, folder, neuron_names, type)
  if type == 1 %plot both plots
    %Plot unweighted adjacency matrix in spy plot
    image1 = figure('Visible', 'off');
    spy(data); %plots adjacency matrix as a sparse plot
    title(names_list, 'interpreter', 'none'); %adds title
    xlabel ('GPCR expression')
    ylabel ('NPP expression')
    drawnow; %draws figure
    saveas(image1, ['Results/scRNAseq_C.elegans/' folder '/Amatrix' names_list '.pdf']);%save image as png
    %Plot the network nodes with the directed edges connecting them 
      if size(data, 1) == size(data, 2) %to make sure digraphs are only made of full networks not sections of them
        image3 = figure('Visible', 'off');
        Digraph_plot = digraph(transpose(data), 'omitselfloops'); %plot digraph with the transposed matrix (row NPP, column receptor), omit selfloops in nodes that do not make biological sense and appear because neuron releases the NPP for a receptor it expresses
        Digraph_plot.Nodes.Name = neuron_names; %make neuron names as node names
        plot(Digraph_plot, 'NodeLabel', Digraph_plot.Nodes.Name, 'NodeFontSize', 5); %plot digraph with the node labels as neuron names
        saveas(image3, ['Results/scRNAseq_C.elegans/' folder '/network' names_list '.pdf']);%save image as png
      end
  elseif type == 0 %plot only spy plot 
    %Plot unweighted adjacency matrix in spy plot
    image1 = figure('Visible', 'off');
    spy(data); %plots adjacency matrix as a sparse plot
    title(names_list, 'interpreter', 'none'); %adds title
    xlabel ('GPCR expression')
    ylabel ('NPP expression')
    drawnow; %draws figure
    saveas(image1, ['Results/scRNAseq_C.elegans/' folder '/' names_list '.pdf']);%save image as png
  end
  close(figure);
end
function [logical_values_005, pvalues_005, stats_005, logical_values_001, pvalues_001, stats_001, Contingency_table_NPP_GPCR, pFDR_values_005_reshape, pFDR_values_001_reshape, logical_values_005_corrected, logical_values_001_corrected] = Fishertestwithtable (data, NPP_names, GPCR_names, name, xaxes, yaxes)
    for i = 1:height(NPP_names)
        for j = 1:height(GPCR_names)
            neuropeptide = NPP_names(i); %select one by one each neuropeptide
            GPCR = GPCR_names(j); %select one by one each receptor
            Neurons_with_NPP_Fisher = find(table2array(data(strcmpi(data.(1), neuropeptide), 2:end))); %get the id of neurons that express that NPP (non 0 elements in expression vector)
            Neurons_without_NPP_Fisher = find(~table2array(data(strcmpi(data.(1), neuropeptide), 2:end))); %get the id of neurons that do not express that NPP ( number of 0 elements in expression vector)
            Neurons_with_NPPandGPCR_Fisher = nnz(table2array(data(strcmpi(data.(1), GPCR), (Neurons_with_NPP_Fisher+1)))); %get the number of neurons that express that NPP & GPCR counting non 0 elements in GPCR expressiong vector for neurons that express NPP 
            Neurons_with_NPP_noGPCR_Fisher = nnz(~table2array(data(strcmpi(data.(1), GPCR), (Neurons_with_NPP_Fisher+1)))); %get the number of neurons that express NPP but do not express that GPCR counting 0 elements in GPCR expression vector for neurons that express NPP
            Neurons_with_GPCR_noNPP_Fisher = nnz(table2array(data(strcmpi(data.(1), GPCR), (Neurons_without_NPP_Fisher+1)))); %get the number of neurons that do not express that NPP but do express that GPCR counting non 0 elements in GPCR expressiong vector for neurons that do not express NPP 
            Neurons_without_NPPandGPCR_Fisher = nnz(~table2array(data(strcmpi(data.(1), GPCR),(Neurons_without_NPP_Fisher+1)))); %get the number of neurons that do not express NPP & GPCR counting 0 elements in GPCR expression vector for neurons that do not express NPP
            Contingency_table_NPP_GPCR{i,j,:} = [Neurons_with_NPPandGPCR_Fisher, Neurons_with_GPCR_noNPP_Fisher; Neurons_with_NPP_noGPCR_Fisher, Neurons_without_NPPandGPCR_Fisher]; %plot contingency table 
            [logical_values_005(i,j), pvalues_005(i,j), stats_005{i,j}] = fishertest([Neurons_with_NPPandGPCR_Fisher, Neurons_with_GPCR_noNPP_Fisher; Neurons_with_NPP_noGPCR_Fisher, Neurons_without_NPPandGPCR_Fisher], 'Alpha', 0.05); % run Fisher's exact test and yield result, p-value and stats, significance 0.05.
            [logical_values_001(i,j), pvalues_001(i,j), stats_001{i,j}] = fishertest([Neurons_with_NPPandGPCR_Fisher, Neurons_with_GPCR_noNPP_Fisher; Neurons_with_NPP_noGPCR_Fisher, Neurons_without_NPPandGPCR_Fisher], 'Alpha', 0.01); % run Fisher's exact test and yield result, p-value and stats, significance 0.01.
        end
    end
    pFDR_values_005(1, :) = mafdr(reshape(pvalues_005.',1,[])); %calculate the positive false discovery rate for each p-value, make matrix into vector
    pFDR_values_001(1, :) = mafdr(reshape(pvalues_001.',1,[])); %calculate the positive false discovery rate for each p-value
    pFDR_values_005_reshape = reshape(pFDR_values_005, [height(GPCR_names), height(NPP_names)]).'; %reshape  into matrix keep all q-values
    pFDR_values_001_reshape = reshape(pFDR_values_001, [height(GPCR_names), height(NPP_names)]).'; %reshape  into matrix keep all q-values
    logical_values_005_corrected = reshape((pFDR_values_005 <= 0.05), [height(GPCR_names), height(NPP_names)]).'; %corrected to 0.05 false discovery rate and reshape  into matrix
    logical_values_001_corrected = reshape((pFDR_values_001 <= 0.05), [height(GPCR_names), height(NPP_names)]).'; %corrected to 0.05 false discovery rate and reshape into matrix
end
function plot_coexpression_data(pFDR_values, logical_values, labely, labelx, connections1, connections2, folder, name)

data = pFDR_values.*logical_values;
if isequal(length(data), height(data))
    data1 = log(data);
    data1(isinf(data1))= 0;
    data1 = tril(data1, -1);
    data1(logical(eye(size(data1)))) = 1;
    data2 = tril(logical_values, -1);
    data3 = cell(size(data2));
    for i = 1:(size(data2,2))
        for j = 1:(size(data2,1))
            if data2(i,j)== 1
                data3{i,j} = {labely{i}, labely{j}};
            end
        end
    end
    
    NPP_coexp_names = vertcat(data3(find(~cellfun(@isempty, data3))));
    
    connections_split_byrec = cell(length(labelx), 1);
    for l = 1:length(labelx)
        t = regexpi(connections1, labelx{l});
        connections_split_byrec{l, 1} = connections2(find(~cellfun(@isempty, t)));
    end
    
    NPP_coexp_coGPCR = false(nnz(data2), length(connections_split_byrec), 2);
    for k = 1:2
        for h = 1:nnz(data2)
            for m = 1:length(connections_split_byrec)
            NPP_coexp_coGPCR(h, m, k) = any(strcmpi(NPP_coexp_names{h}{k}, connections_split_byrec{m}));
            end
        end
    end
    
    % Calculate the number of occurrences where both k=1 and k=2 are true (logical AND)
    NPP_coexp_coGPCR_n = sum((NPP_coexp_coGPCR(:,:,1) + NPP_coexp_coGPCR(:, :, 2)) == 2, 2) + 1; %add 1 to make sure all coexp events are represented (in legend substract 1)
    [xG, yG] = find(data2 == 1);
    data4 = accumarray([xG(:),yG(:)],NPP_coexp_coGPCR_n(:),[size(data2)]); %represent the sizes in the matrix format to use them as size for the plot
    data4(logical(eye(size(data4)))) = 1; %make the diagonal 1 to represent all the self coexpression
    bubbleSizes = abs(data4(:)); % Adjust the scaling factor as needed
    bubbleSizes(bubbleSizes == 0) = nan;

    ticks = {labely; labely};
else
    data1 = log(data);
    data1(isinf(data1))= 0;
    bubbleSizes = abs(data1(:)); % Adjust the scaling factor as needed
    bubbleSizes(bubbleSizes == 0) = nan;
    ticks = {labely; labelx};
end

[x, y] = meshgrid(1:size(data1,2), 1:size(data1,1));
x(data1 == 0 ) = nan;
y(data1 == 0 ) = nan;
bubbleColors = data1(:);
bubbleColors(bubbleColors == 0) = nan;

%Plot bubble chart
figure('Visible', 'off')
bubblechart(x(:), y(:), bubbleSizes, bubbleColors)
bubblesize([5 15])
bubblelegend('Cooperativity','Location','northeastoutside')

% Cosmetics
colormap('winter')
grid on
box off
set(gca,'xtick', 1:size(data1,2), ...
    'ytick', 1:size(data1,1), ...
    'YDir', 'Reverse'); % typically corr matrices use flipped y axes
yticklabels(ticks{1})
xticklabels(ticks{2})
xlabel('GPCR receptors')
ylabel('Neuropeptides')
cb = colorbar; 
title(cb, 'Log corrected p-values ')
set(cb, 'Location', 'eastoutside', 'Position', [0.8 0.5 0.06 0.18])
set(findall(gcf,'-property','FontSize'),'FontSize',8)

print(['Results/scRNAseq_C.elegans/' folder '/Coexpression_' name], '-dpdf', '-bestfit', '-r300');

end