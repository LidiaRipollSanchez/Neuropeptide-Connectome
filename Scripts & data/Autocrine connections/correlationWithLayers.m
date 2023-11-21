%% Script to analyze the correlation with NP - synaptic - gap and MA degree
% Jan Watteyne - jan.watteyne@kuleuven.be and/or jwatteyne@gmail.com
% 20230607

close all
clear all

%% Import data
%-------------------------------------------
% import general variables from network build
%-------------------------------------------
% import gene ID
load ('05062023_Networks_build_EC500nM.mat','geneID');
% import neuron_names
load ('05062023_Networks_build_EC500nM.mat','neuron_names');
% import anatomical_class
load ('05062023_Networks_build_EC500nM.mat','anatomy');

% import synaptic adjacency matrix(load from .mat file Barry Bentley)
load ('10022022_Network_measures_EC500nM.mat','synA');

% import gap junction adjacency matrix(load from .mat file Lidia)
load ('10022022_Network_measures_EC500nM.mat','gapA');
% 
% % import monoamine adjacency matrix(load from .mat file Lidia)
load ('05062023_Networks_build_EC500nM.mat','monoamine_networks_all');

% import Neuropeptide adjacency matrix - short - mid - long range(load from .mat file Lidia)
load ('05062023_Networks_build_EC500nM.mat','NPP_all_grouped_networks');
load ('05062023_Networks_build_EC500nM.mat','NPP_all_grouped_networks_short_range_connections');
load ('05062023_Networks_build_EC500nM.mat','NPP_all_grouped_networks_mid_range_connections');

%----------------------------------
% import CeNGEN expression matrices
%----------------------------------
expression_threshold4 = readtable ('CeNGEN_threshold4.csv'); 

%-------------------------------------------------------------------------------------------
% import NPP-GPCR connections with their respective  EC50 value (Beets et al., 2022 bioRxiv)
%-------------------------------------------------------------------------------------------
connections_EC50 = readtable('Data S6_WBID.xlsx');

% using geneID list, convert WB_ID to common gene name
for row = 1:size(connections_EC50, 1)
    try
    connections_EC50.GPCR_geneName(row) = geneID.Gene(find(strcmp(connections_EC50.GPCRWBID(row),geneID.Var1)));
    catch
        connections_EC50.GPCR_geneName(row) = {'dmsr-5'};
    end
end

for row = 1:size(connections_EC50, 1)
    try
    connections_EC50.NPP_geneName(row) = geneID.Gene(find(strcmp(connections_EC50.NPPWBID(row),geneID.Var1)));
    catch end
end

% subset original table to only retain connections falling below said EC50
% threshold
connections_EC50_500nM = connections_EC50(connections_EC50.EC50_M_ <=  500e-09, :); 

% each of these tables will have multiple entries for specific NPP-GPCR
% pairs (because multiple isoforms/peptides can be produced from one gene).
% Therefore, only retain unique combinations

connections_EC50_500nM = table(connections_EC50_500nM.GPCR_geneName, connections_EC50_500nM.NPP_geneName);
connections_EC50_500nM = unique(connections_EC50_500nM,'rows'); 

%% Convert expression data to binary expression matrix using functions Lidia Ripoll-Sanchez
[expression_threshold4_ordered] = Adapt_expression(expression_threshold4, neuron_names, geneID); %Adapt the expression values split by neuron type to neuron (from 180 columns to 302 columns) % J: column number doesn't change
[expression_threshold4_ordered_unweighted] = Mask_network(expression_threshold4_ordered); %Unweight expression values, when expression value is 1 otherwise is 0

%% calculate number of self loops for each neuron
% original code Lidia 
for i = 1:height(connections_EC50_500nM) % iterate over each NP - GPCR pair
    NPP = table2array(expression_threshold4_ordered_unweighted(strcmpi(expression_threshold4_ordered_unweighted.(1), char(connections_EC50_500nM.Var2(i))), 2:303)); %get the vector for each neuropeptide
    GPCR = table2array(expression_threshold4_ordered_unweighted(strcmpi(expression_threshold4_ordered_unweighted.(1), char(connections_EC50_500nM.Var1(i))), 2:303)); %get the vector for each receptor
    pair = NPP + GPCR;
    index = find(pair >= 2); % get index of self loop neurons
    self_loops_per_neuron{i} = neuron_names(index); % for each NP - GPCR couple, display all the neurons in which a self loop is present
    self_loops_per_pair{i} = sum(pair(:) >= 2); %count total number of neurons that show a self loop for that particular NP - GPCR couple
end
clear NPP
clear GPCR
clear pair
clear index
clear i

%% Compute the number of different self loops each neuron has
neuronList = [];
for i = 1:length(self_loops_per_neuron) % concatenate all the neuron names in one big cell
    neuronList = [neuronList, self_loops_per_neuron{1, i}'];
end

% count how many times each neuron name is found in this list + bundle all
% data in a table
[NeuronIDLoops,~,ic] = unique(neuronList); 
neuronList_counts = accumarray(ic,1);

% also add all the neurons which don't have any self loops 
NeuronsWithoutSelfLoops = setdiff(neuron_names, NeuronIDLoops); %retrieve names of neurons without self loop
zeroCounts_otherNeurons = zeros(length(NeuronsWithoutSelfLoops), 1);

%tableNumberLoopsForNeuron = table(NeuronIDLoops', neuronList_counts);
tableNumberLoopsForNeuron = table([NeuronIDLoops, NeuronsWithoutSelfLoops']', [neuronList_counts; zeroCounts_otherNeurons]);
tableNumberLoopsForNeuron.Properties.VariableNames = {'NeuronID' 'HowManyLoops'};

tableNumberLoopsForNeuron = sortrows(tableNumberLoopsForNeuron, 2); % order the table in ascending order

% figure('Renderer', 'painters', 'Position', [10 10 2500 600])
% plot(tableNumberLoopsForNeuron.HowManyLoops)
% xlim([1 height(tableNumberLoopsForNeuron)])
% ylabel('Number of self loops for each neuron')
% xticks(1:height(tableNumberLoopsForNeuron))
% xticklabels(string(tableNumberLoopsForNeuron.NeuronID)')
% ax=gca;
% ax.XAxis.FontSize = 6;
% xtickangle(90)
% set(gcf,'color','w');
% %export_fig 'NumberLoopsForEachNeuron.pdf' -painters

%% Correlation number of self loops with neuropeptidergic degree
% retrieves adjacency matrices: x-axes is neuropeptide expression, y-axes
% GPCR expression
AdjMatrix = NPP_all_grouped_networks_short_range_connections; % select which type of connection you want (short, mid-range, all)

%All_grouped_networks_short_range_connections;
BinaryAdjMatrix = AdjMatrix; % binarize the adjacency matrix
BinaryAdjMatrix(BinaryAdjMatrix~= 0) = 1;

%spy(BinaryAdjMatrix), xlabel('neuropeptide expression'), ylabel('GPCR expression')
NPInDegree = sum(BinaryAdjMatrix, 2);
NPOutDegree = sum(BinaryAdjMatrix, 1);
NPTotalDegree = NPInDegree + NPOutDegree';

% add NP in/out/total degree to table
tableNumberLoopsForNeuron.NPInDegree = nan(1, size(tableNumberLoopsForNeuron, 1))';
tableNumberLoopsForNeuron.NPOutDegree = nan(1, size(tableNumberLoopsForNeuron, 1))';
tableNumberLoopsForNeuron.NPTotalDegree = nan(1, size(tableNumberLoopsForNeuron, 1))';
for row = 1:size(tableNumberLoopsForNeuron, 1) % iterate over each neuron within the table
    try
    neuron = string(tableNumberLoopsForNeuron.NeuronID(row));
    
    % find neuron class for that particular neuron, and add it to the
    % variable tableNumberLoopsForNeuron
    index = find(strcmp(string(neuron_names), neuron));
    
    tableNumberLoopsForNeuron.NPInDegree(row) = NPInDegree(index);
    tableNumberLoopsForNeuron.NPOutDegree(row) = NPOutDegree(index);
    tableNumberLoopsForNeuron.NPTotalDegree(row) = NPTotalDegree(index);
    catch end
end

figure('Renderer', 'painters', 'Position', [10 10 900 600])

subplot(231)
plot(tableNumberLoopsForNeuron.HowManyLoops, tableNumberLoopsForNeuron.NPTotalDegree, 'b.')
ylim([0 600])
xlabel('number of self loops')
ylabel('NP Total degree')
lsline % plot regression line
mdl = fitlm(tableNumberLoopsForNeuron.HowManyLoops,tableNumberLoopsForNeuron.NPTotalDegree); % fit linear regression model
Rsquared = mdl.Rsquared.Ordinary;
RsquaredAdj = mdl.Rsquared.Adjusted;
set(gcf,'color','w');

str = ['Adj. R-squared = ',num2str(RsquaredAdj)];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 8, 'verticalalignment', 'top', 'horizontalalignment', 'left');

subplot(234)
% retrieve NPTotalDegree for neurons with and without self loops
NPDegreeSelfLoops = tableNumberLoopsForNeuron.NPTotalDegree(tableNumberLoopsForNeuron.HowManyLoops ~= 0);
NPDegreeNoSelfLoops = tableNumberLoopsForNeuron.NPTotalDegree(tableNumberLoopsForNeuron.HowManyLoops == 0);

% format data for plotting
numrows = max([length(NPDegreeSelfLoops), length(NPDegreeNoSelfLoops)]);
DataMatrix = nan(numrows, 2);
% insert data in DataMatrix
DataMatrix(1:length(NPDegreeNoSelfLoops), 1) = NPDegreeNoSelfLoops';
DataMatrix(1:length(NPDegreeSelfLoops), 2) = NPDegreeSelfLoops';

h = notBoxPlot(DataMatrix,[], [], 'sdline'); %external function
h(1).data.MarkerSize = 4;
h(2).data.MarkerSize = 4;
set(gca,'XTickLabel',{'no self loops', 'self loops'}, 'XTickLabelRotation', 45);
ylabel('NP total degree')
ylim([0 600])

subplot(232)
plot(tableNumberLoopsForNeuron.HowManyLoops, tableNumberLoopsForNeuron.NPInDegree, 'b.')
ylim([0 600])
xlabel('number of self loops')
ylabel('NP In degree')
lsline % plot regression line
mdl = fitlm(tableNumberLoopsForNeuron.HowManyLoops,tableNumberLoopsForNeuron.NPInDegree); % fit linear regression model
Rsquared = mdl.Rsquared.Ordinary;
RsquaredAdj = mdl.Rsquared.Adjusted;
str = ['Adj. R-squared = ',num2str(RsquaredAdj)];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 8, 'verticalalignment', 'top', 'horizontalalignment', 'left');

subplot(235)
% retrieve NPInDegree for neurons with and without self loops
NPDegreeSelfLoops = tableNumberLoopsForNeuron.NPInDegree(tableNumberLoopsForNeuron.HowManyLoops ~= 0);
NPDegreeNoSelfLoops = tableNumberLoopsForNeuron.NPInDegree(tableNumberLoopsForNeuron.HowManyLoops == 0);

% format data for plotting
numrows = max([length(NPDegreeSelfLoops), length(NPDegreeNoSelfLoops)]);
DataMatrix = nan(numrows, 2);
% insert data in DataMatrix
DataMatrix(1:length(NPDegreeNoSelfLoops), 1) = NPDegreeNoSelfLoops';
DataMatrix(1:length(NPDegreeSelfLoops), 2) = NPDegreeSelfLoops';

h = notBoxPlot(DataMatrix,[], [], 'sdline'); %external function
h(1).data.MarkerSize = 4;
h(2).data.MarkerSize = 4;
set(gca,'XTickLabel',{'no self loops', 'self loops'}, 'XTickLabelRotation', 45);
ylabel('NP In degree')
ylim([0 600])

subplot(233)
plot(tableNumberLoopsForNeuron.HowManyLoops, tableNumberLoopsForNeuron.NPOutDegree, 'b.')
ylim([0 600])
xlabel('number of self loops')
ylabel('NP Out degree')
lsline % plot regression line
mdl = fitlm(tableNumberLoopsForNeuron.HowManyLoops,tableNumberLoopsForNeuron.NPOutDegree); % fit linear regression model
Rsquared = mdl.Rsquared.Ordinary;
RsquaredAdj = mdl.Rsquared.Adjusted;

str = ['Adj. R-squared = ',num2str(RsquaredAdj)];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 8, 'verticalalignment', 'top', 'horizontalalignment', 'left');

subplot(236)
% retrieve NPOutDegree for neurons with and without self loops
NPDegreeSelfLoops = tableNumberLoopsForNeuron.NPOutDegree(tableNumberLoopsForNeuron.HowManyLoops ~= 0);
NPDegreeNoSelfLoops = tableNumberLoopsForNeuron.NPOutDegree(tableNumberLoopsForNeuron.HowManyLoops == 0);

% format data for plotting
numrows = max([length(NPDegreeSelfLoops), length(NPDegreeNoSelfLoops)]);
DataMatrix = nan(numrows, 2);
% insert data in DataMatrix
DataMatrix(1:length(NPDegreeNoSelfLoops), 1) = NPDegreeNoSelfLoops';
DataMatrix(1:length(NPDegreeSelfLoops), 2) = NPDegreeSelfLoops';

h = notBoxPlot(DataMatrix,[], [], 'sdline'); %external function
h(1).data.MarkerSize = 4;
h(2).data.MarkerSize = 4;
set(gca,'XTickLabel',{'no self loops', 'self loops'}, 'XTickLabelRotation', 45);
ylabel('NP Out degree')
ylim([0 600])

%export_fig 'CorrelationWithNPDegree_All.pdf' -painters

%% make a composite figure for the correlation with NP degree:
figure('Renderer', 'painters', 'Position', [10 10 600 600])

% NP degree OUT
plot(tableNumberLoopsForNeuron.HowManyLoops, tableNumberLoopsForNeuron.NPOutDegree, '<', 'Color', [0.3882 0.3882 0.3882], 'MarkerFaceColor',[0.3882 0.3882 0.3882]), hold on
xlim([0 15])
xticks(0:5:15)
xticklabels('')
ylim([0 600])
yticks(0:200:600)
yticklabels('')
h = lsline; % plot regression line
mdl = fitlm(tableNumberLoopsForNeuron.HowManyLoops,tableNumberLoopsForNeuron.NPInDegree); % fit linear regression model
Rsquared = mdl.Rsquared.Ordinary;
RsquaredAdj = mdl.Rsquared.Adjusted;

% NP degree IN
plot(tableNumberLoopsForNeuron.HowManyLoops, tableNumberLoopsForNeuron.NPInDegree, '>', 'Color', [0.7412 0.7412 0.7412], 'MarkerFaceColor', [0.7412 0.7412 0.7412])
h = lsline; % plot regression line
mdl = fitlm(tableNumberLoopsForNeuron.HowManyLoops,tableNumberLoopsForNeuron.NPInDegree); % fit linear regression model
Rsquared = mdl.Rsquared.Ordinary;
RsquaredAdj = mdl.Rsquared.Adjusted;

% all NP degree
plot(tableNumberLoopsForNeuron.HowManyLoops, tableNumberLoopsForNeuron.NPTotalDegree, 'ko', 'MarkerFaceColor', 'k')
h = lsline; % plot regression line
mdl = fitlm(tableNumberLoopsForNeuron.HowManyLoops,tableNumberLoopsForNeuron.NPTotalDegree); % fit linear regression model
Rsquared = mdl.Rsquared.Ordinary;
RsquaredAdj = mdl.Rsquared.Adjusted;

h(3).LineWidth = 2;
h(3).LineStyle = '--';
h(2).LineWidth = 2;
h(2).LineStyle = '-.';
h(1).LineWidth = 2;

set(gcf,'color','w');
%export_fig 'MainFig_CorrelationWithNPDegree.pdf' -painters
%% Correlation number of self loops with wired degree

% assume adjacency matrices are binary (the number of different inputs, not
% total number of input)
BinarySynA = synA;
BinarySynA = double(BinarySynA ~=0); %spy(BinarySynA)

BinaryGapA = gapA;
BinaryGapA = double(BinaryGapA ~=0); %spy(BinaryGapA)

synInDegree = sum(BinarySynA, 2);
synOutDegree = sum(BinarySynA, 1);
synTotalDegree = synInDegree + synOutDegree';

gapTotalDegree = sum(BinaryGapA, 2) + sum(BinaryGapA, 1)'; %gap junctions are not directed

% add wired degrees to table
tableNumberLoopsForNeuron.synInDegree = nan(1, size(tableNumberLoopsForNeuron, 1))';
tableNumberLoopsForNeuron.synOutDegree = nan(1, size(tableNumberLoopsForNeuron, 1))';
tableNumberLoopsForNeuron.synTotalDegree = nan(1, size(tableNumberLoopsForNeuron, 1))';
tableNumberLoopsForNeuron.gapTotalDegree = nan(1, size(tableNumberLoopsForNeuron, 1))';
for row = 1:size(tableNumberLoopsForNeuron, 1) % iterate over each neuron within the table
    try
    neuron = string(tableNumberLoopsForNeuron.NeuronID(row));
    
    % find neuron class for that particular neuron, and add it to the
    % variable tableNumberLoopsForNeuron
    index = find(strcmp(string(neuron_names), neuron));
    
    tableNumberLoopsForNeuron.synInDegree(row) = synInDegree(index);
    tableNumberLoopsForNeuron.synOutDegree(row) = synOutDegree(index);
    tableNumberLoopsForNeuron.synTotalDegree(row) = synTotalDegree(index);
    tableNumberLoopsForNeuron.gapTotalDegree(row) = gapTotalDegree(index);
    catch end
end

figure('Renderer', 'painters', 'Position', [10 10 1200 600])

subplot(241)
plot(tableNumberLoopsForNeuron.HowManyLoops, tableNumberLoopsForNeuron.synTotalDegree, 'b.')
xlabel('number of self loops')
ylabel('Synaptic Total degree')
lsline % plot regression line
mdl = fitlm(tableNumberLoopsForNeuron.HowManyLoops,tableNumberLoopsForNeuron.synTotalDegree); % fit linear regression model
%Rsquared = mdl.Rsquared.Ordinary;
RsquaredAdj = mdl.Rsquared.Adjusted;
xlim([0 15])
ylim([0 100])

str = ['Adj. R-squared = ',num2str(RsquaredAdj)];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 8, 'verticalalignment', 'top', 'horizontalalignment', 'left');

subplot(245)
% retrieve SynTotalDegree for neurons with and without self loops
synTotalDegreeSelfLoops = tableNumberLoopsForNeuron.synTotalDegree(tableNumberLoopsForNeuron.HowManyLoops ~= 0);
synTotalDegreeNoSelfLoops = tableNumberLoopsForNeuron.synTotalDegree(tableNumberLoopsForNeuron.HowManyLoops == 0);

% format data for plotting
numrows = max([length(synTotalDegreeSelfLoops), length(synTotalDegreeNoSelfLoops)]);
DataMatrix = nan(numrows, 2);
% insert data in DataMatrix
DataMatrix(1:length(synTotalDegreeNoSelfLoops), 1) = synTotalDegreeNoSelfLoops';
DataMatrix(1:length(synTotalDegreeSelfLoops), 2) = synTotalDegreeSelfLoops';

h = notBoxPlot(DataMatrix,[], [], 'sdline'); %external function
h(1).data.MarkerSize = 4;
h(2).data.MarkerSize = 4;
set(gca,'XTickLabel',{'no self loops', 'self loops'}, 'XTickLabelRotation', 45);
ylabel('Synaptic Total degree')
ylim([0 100])

subplot(242)
plot(tableNumberLoopsForNeuron.HowManyLoops, tableNumberLoopsForNeuron.synInDegree, 'b.')
xlabel('number of self loops')
ylabel('Synaptic In degree')
lsline % plot regression line
mdl = fitlm(tableNumberLoopsForNeuron.HowManyLoops,tableNumberLoopsForNeuron.synInDegree); % fit linear regression model
%Rsquared = mdl.Rsquared.Ordinary;
RsquaredAdj = mdl.Rsquared.Adjusted;
xlim([0 15])
ylim([0 100])

str = ['Adj. R-squared = ',num2str(RsquaredAdj)];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 8, 'verticalalignment', 'top', 'horizontalalignment', 'left');

subplot(246)
% retrieve SynTotalDegree for neurons with and without self loops
synInDegreeSelfLoops = tableNumberLoopsForNeuron.synInDegree(tableNumberLoopsForNeuron.HowManyLoops ~= 0);
synInDegreeNoSelfLoops = tableNumberLoopsForNeuron.synInDegree(tableNumberLoopsForNeuron.HowManyLoops == 0);

% format data for plotting
numrows = max([length(synInDegreeSelfLoops), length(synInDegreeNoSelfLoops)]);
DataMatrix = nan(numrows, 2);
% insert data in DataMatrix
DataMatrix(1:length(synInDegreeNoSelfLoops), 1) = synInDegreeNoSelfLoops';
DataMatrix(1:length(synInDegreeSelfLoops), 2) = synInDegreeSelfLoops';

h = notBoxPlot(DataMatrix,[], [], 'sdline'); %external function
h(1).data.MarkerSize = 4;
h(2).data.MarkerSize = 4;
set(gca,'XTickLabel',{'no self loops', 'self loops'}, 'XTickLabelRotation', 45);
ylabel('Synaptic In Degree')
ylim([0 100])

subplot(243)
plot(tableNumberLoopsForNeuron.HowManyLoops, tableNumberLoopsForNeuron.synOutDegree, 'b.')
xlabel('number of self loops')
ylabel('Synaptic Out degree')
lsline % plot regression line
mdl = fitlm(tableNumberLoopsForNeuron.HowManyLoops,tableNumberLoopsForNeuron.synOutDegree); % fit linear regression model
Rsquared = mdl.Rsquared.Ordinary;
RsquaredAdj = mdl.Rsquared.Adjusted;
xlim([0 15])
ylim([0 100])

str = ['Adj. R-squared = ',num2str(RsquaredAdj)];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 8, 'verticalalignment', 'top', 'horizontalalignment', 'left');

subplot(247)
% retrieve SynTotalDegree for neurons with and without self loops
synOutDegreeSelfLoops = tableNumberLoopsForNeuron.synOutDegree(tableNumberLoopsForNeuron.HowManyLoops ~= 0);
synOutDegreeNoSelfLoops = tableNumberLoopsForNeuron.synOutDegree(tableNumberLoopsForNeuron.HowManyLoops == 0);

% format data for plotting
numrows = max([length(synOutDegreeSelfLoops), length(synOutDegreeNoSelfLoops)]);
DataMatrix = nan(numrows, 2);
% insert data in DataMatrix
DataMatrix(1:length(synOutDegreeNoSelfLoops), 1) = synOutDegreeNoSelfLoops';
DataMatrix(1:length(synOutDegreeSelfLoops), 2) = synOutDegreeSelfLoops';

h = notBoxPlot(DataMatrix,[], [], 'sdline'); %external function
h(1).data.MarkerSize = 4;
h(2).data.MarkerSize = 4;
set(gca,'XTickLabel',{'no self loops', 'self loops'}, 'XTickLabelRotation', 45);
ylabel('Synaptic Out Degree')
ylim([0 100])

subplot(244)
plot(tableNumberLoopsForNeuron.HowManyLoops, tableNumberLoopsForNeuron.gapTotalDegree, 'b.')
xlabel('number of self loops')
ylabel('Gap Total degree')
lsline % plot regression line
mdl = fitlm(tableNumberLoopsForNeuron.HowManyLoops,tableNumberLoopsForNeuron.gapTotalDegree); % fit linear regression model
Rsquared = mdl.Rsquared.Ordinary;
RsquaredAdj = mdl.Rsquared.Adjusted;
xlim([0 15])
ylim([0 100])

str = ['Adj. R-squared = ',num2str(RsquaredAdj)];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 8, 'verticalalignment', 'top', 'horizontalalignment', 'left');

subplot(248)
% retrieve SynTotalDegree for neurons with and without self loops
gapTotalDegreeSelfLoops = tableNumberLoopsForNeuron.gapTotalDegree(tableNumberLoopsForNeuron.HowManyLoops ~= 0);
gapTotalDegreeNoSelfLoops = tableNumberLoopsForNeuron.gapTotalDegree(tableNumberLoopsForNeuron.HowManyLoops == 0);

% format data for plotting
numrows = max([length(gapTotalDegreeSelfLoops), length(gapTotalDegreeNoSelfLoops)]);
DataMatrix = nan(numrows, 2);
% insert data in DataMatrix
DataMatrix(1:length(gapTotalDegreeNoSelfLoops), 1) = gapTotalDegreeNoSelfLoops';
DataMatrix(1:length(gapTotalDegreeSelfLoops), 2) = gapTotalDegreeSelfLoops';

h = notBoxPlot(DataMatrix,[], [], 'sdline'); %external function
h(1).data.MarkerSize = 4;
h(2).data.MarkerSize = 4;
set(gca,'XTickLabel',{'no self loops', 'self loops'}, 'XTickLabelRotation', 45);
ylabel('Gap Total Degree')
ylim([0 100])

set(gcf,'color','w');
%export_fig 'CorrelationWithWiredDegree.pdf' -painters

%% make a composite figure for the correlation with wired degree (synapse):
figure('Renderer', 'painters', 'Position', [10 10 600 600])

% NP degree OUT
plot(tableNumberLoopsForNeuron.HowManyLoops, tableNumberLoopsForNeuron.synOutDegree, '<', 'Color', [0.3882 0.3882 0.3882], 'MarkerFaceColor',[0.3882 0.3882 0.3882]), hold on
xlim([0 15])
xticks(0:5:15)
xticklabels('')
ylim([0 100])
yticks(0:25:100)
yticklabels('')
h = lsline; % plot regression line
mdl = fitlm(tableNumberLoopsForNeuron.HowManyLoops,tableNumberLoopsForNeuron.synOutDegree); % fit linear regression model
Rsquared = mdl.Rsquared.Ordinary;
RsquaredAdj = mdl.Rsquared.Adjusted;

% NP degree IN
plot(tableNumberLoopsForNeuron.HowManyLoops, tableNumberLoopsForNeuron.synInDegree, '>', 'Color', [0.7412 0.7412 0.7412], 'MarkerFaceColor', [0.7412 0.7412 0.7412])
h = lsline; % plot regression line
mdl = fitlm(tableNumberLoopsForNeuron.HowManyLoops,tableNumberLoopsForNeuron.synInDegree); % fit linear regression model
Rsquared = mdl.Rsquared.Ordinary;
RsquaredAdj = mdl.Rsquared.Adjusted;

% all NP degree
plot(tableNumberLoopsForNeuron.HowManyLoops, tableNumberLoopsForNeuron.synTotalDegree, 'ko', 'MarkerFaceColor', 'k')
h = lsline; % plot regression line
mdl = fitlm(tableNumberLoopsForNeuron.HowManyLoops,tableNumberLoopsForNeuron.synTotalDegree); % fit linear regression model
Rsquared = mdl.Rsquared.Ordinary;
RsquaredAdj = mdl.Rsquared.Adjusted;

h(3).LineWidth = 2;
h(3).LineStyle = '--';
h(2).LineWidth = 2;
h(2).LineStyle = '-.';
h(1).LineWidth = 2;

set(gcf,'color','w');
%export_fig 'MainFig_CorrelationWithSynapticDegree.pdf' -painters

%% make a composite figure for the correlation with wired degree (gap):
figure('Renderer', 'painters', 'Position', [10 10 600 600])

% Gap degree 
plot(tableNumberLoopsForNeuron.HowManyLoops, tableNumberLoopsForNeuron.gapTotalDegree, 'o', 'Color', 'k', 'MarkerFaceColor','k'), hold on
xlim([0 15])
xticks(0:5:15)
xticklabels('')
ylim([0 100])
yticks(0:25:100)
yticklabels('')
h = lsline; % plot regression line
mdl = fitlm(tableNumberLoopsForNeuron.HowManyLoops,tableNumberLoopsForNeuron.synOutDegree); % fit linear regression model
Rsquared = mdl.Rsquared.Ordinary;
RsquaredAdj = mdl.Rsquared.Adjusted;

h(1).LineWidth = 2;

set(gcf,'color','w');
%export_fig 'SupplFig_CorrelationWithGapDegree.pdf' -painters
%% Correlation number of self loops with monoamines

AllMonoamineA = zeros(size(monoamine_networks_all{1}));
for ii = 1:size(monoamine_networks_all, 1)
    AllMonoamineA = AllMonoamineA + monoamine_networks_all{ii};
end
BinaryAllMonoamineA = AllMonoamineA;
BinaryAllMonoamineA = double(BinaryAllMonoamineA ~=0); %spy(BinaryAllMonoamineA)

MAInDegree = sum(BinaryAllMonoamineA, 1);
MAOutDegree = sum(BinaryAllMonoamineA, 2);
MATotalDegree = MAInDegree + MAOutDegree';

% add MonoAmine degrees to table
tableNumberLoopsForNeuron.MAInDegree = nan(1, size(tableNumberLoopsForNeuron, 1))';
tableNumberLoopsForNeuron.MAOutDegree = nan(1, size(tableNumberLoopsForNeuron, 1))';
tableNumberLoopsForNeuron.MATotalDegree = nan(1, size(tableNumberLoopsForNeuron, 1))';
for row = 1:size(tableNumberLoopsForNeuron, 1) % iterate over each neuron within the table
    try
    neuron = string(tableNumberLoopsForNeuron.NeuronID(row));
    
    % find neuron class for that particular neuron, and add it to the
    % variable tableNumberLoopsForNeuron
    index = find(strcmp(string(neuron_names), neuron));
    
    tableNumberLoopsForNeuron.MAInDegree(row) = MAInDegree(index);
    tableNumberLoopsForNeuron.MAOutDegree(row) = MAOutDegree(index);
    tableNumberLoopsForNeuron.MATotalDegree(row) = MATotalDegree(index);
    catch end
end

figure('Renderer', 'painters', 'Position', [10 10 900 600])
subplot(231)
plot(tableNumberLoopsForNeuron.HowManyLoops, tableNumberLoopsForNeuron.MATotalDegree, 'b.')
ylim([0 200])
xlabel('number of self loops')
ylabel('MA Total degree')
lsline % plot regression line
mdl = fitlm(tableNumberLoopsForNeuron.HowManyLoops,tableNumberLoopsForNeuron.MATotalDegree); % fit linear regression model
Rsquared = mdl.Rsquared.Ordinary;
RsquaredAdj = mdl.Rsquared.Adjusted;

str = ['Adj. R-squared = ',num2str(RsquaredAdj)];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 8, 'verticalalignment', 'top', 'horizontalalignment', 'left');

subplot(234)
MATotalDegreeSelfLoops = tableNumberLoopsForNeuron.MATotalDegree(tableNumberLoopsForNeuron.HowManyLoops ~= 0);
MATotalDegreeNoSelfLoops = tableNumberLoopsForNeuron.MATotalDegree(tableNumberLoopsForNeuron.HowManyLoops == 0);

% format data for plotting
numrows = max([length(MATotalDegreeSelfLoops), length(MATotalDegreeNoSelfLoops)]);
DataMatrix = nan(numrows, 2);
% insert data in DataMatrix
DataMatrix(1:length(MATotalDegreeNoSelfLoops), 1) = MATotalDegreeNoSelfLoops';
DataMatrix(1:length(MATotalDegreeSelfLoops), 2) = MATotalDegreeSelfLoops';

h = notBoxPlot(DataMatrix,[], [], 'sdline'); %external function
h(1).data.MarkerSize = 4;
h(2).data.MarkerSize = 4;
set(gca,'XTickLabel',{'no self loops', 'self loops'}, 'XTickLabelRotation', 45);
ylabel('MA Total Degree')
ylim([0 200])

subplot(232)
plot(tableNumberLoopsForNeuron.HowManyLoops, tableNumberLoopsForNeuron.MAInDegree, 'b.')
ylim([0 200])
xlabel('number of self loops')
ylabel('MA In degree')
lsline % plot regression line
mdl = fitlm(tableNumberLoopsForNeuron.HowManyLoops,tableNumberLoopsForNeuron.MAInDegree); % fit linear regression model
Rsquared = mdl.Rsquared.Ordinary;
RsquaredAdj = mdl.Rsquared.Adjusted;
str = ['Adj. R-squared = ',num2str(RsquaredAdj)];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 8, 'verticalalignment', 'top', 'horizontalalignment', 'left');

subplot(235)
MAInDegreeSelfLoops = tableNumberLoopsForNeuron.MAInDegree(tableNumberLoopsForNeuron.HowManyLoops ~= 0);
MAInDegreeNoSelfLoops = tableNumberLoopsForNeuron.MAInDegree(tableNumberLoopsForNeuron.HowManyLoops == 0);

% format data for plotting
numrows = max([length(MAInDegreeSelfLoops), length(MAInDegreeNoSelfLoops)]);
DataMatrix = nan(numrows, 2);
% insert data in DataMatrix
DataMatrix(1:length(MAInDegreeNoSelfLoops), 1) = MAInDegreeNoSelfLoops';
DataMatrix(1:length(MAInDegreeSelfLoops), 2) = MAInDegreeSelfLoops';

h = notBoxPlot(DataMatrix,[], [], 'sdline'); %external function
h(1).data.MarkerSize = 4;
h(2).data.MarkerSize = 4;
set(gca,'XTickLabel',{'no self loops', 'self loops'}, 'XTickLabelRotation', 45);
ylabel('MA In Degree')
ylim([0 200])

subplot(233)
plot(tableNumberLoopsForNeuron.HowManyLoops, tableNumberLoopsForNeuron.MAOutDegree, 'b.')
ylim([0 200])
xlabel('number of self loops')
ylabel('MA Out degree')
lsline % plot regression line
mdl = fitlm(tableNumberLoopsForNeuron.HowManyLoops,tableNumberLoopsForNeuron.MAOutDegree); % fit linear regression model
Rsquared = mdl.Rsquared.Ordinary;
RsquaredAdj = mdl.Rsquared.Adjusted;

str = ['Adj. R-squared = ',num2str(RsquaredAdj)];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 8, 'verticalalignment', 'top', 'horizontalalignment', 'left');

subplot(236)
MAOutDegreeSelfLoops = tableNumberLoopsForNeuron.MAOutDegree(tableNumberLoopsForNeuron.HowManyLoops ~= 0);
MAOutDegreeNoSelfLoops = tableNumberLoopsForNeuron.MAOutDegree(tableNumberLoopsForNeuron.HowManyLoops == 0);

% format data for plotting
numrows = max([length(MAOutDegreeSelfLoops), length(MAOutDegreeNoSelfLoops)]);
DataMatrix = nan(numrows, 2);
% insert data in DataMatrix
DataMatrix(1:length(MAOutDegreeNoSelfLoops), 1) = MAOutDegreeNoSelfLoops';
DataMatrix(1:length(MAOutDegreeSelfLoops), 2) = MAOutDegreeSelfLoops';

h = notBoxPlot(DataMatrix,[], [], 'sdline'); %external function
h(1).data.MarkerSize = 4;
h(2).data.MarkerSize = 4;
set(gca,'XTickLabel',{'no self loops', 'self loops'}, 'XTickLabelRotation', 45);
ylabel('MA Out Degree')
ylim([0 200])

set(gcf,'color','w');
%export_fig 'CorrelationWithMonoAmineDegree.pdf' -painters
%% make a composite figure for the correlation with Monoamine degree:
figure('Renderer', 'painters', 'Position', [10 10 600 600])

% NP degree OUT
plot(tableNumberLoopsForNeuron.HowManyLoops, tableNumberLoopsForNeuron.MAOutDegree, '<', 'Color', [0.3882 0.3882 0.3882], 'MarkerFaceColor',[0.3882 0.3882 0.3882]), hold on
xlim([0 15])
xticks(0:5:15)
xticklabels('')
ylim([0 300])
yticks(0:100:200)
yticklabels('')
h = lsline; % plot regression line
mdl = fitlm(tableNumberLoopsForNeuron.HowManyLoops,tableNumberLoopsForNeuron.MAOutDegree); % fit linear regression model
Rsquared = mdl.Rsquared.Ordinary;
RsquaredAdj = mdl.Rsquared.Adjusted;

% NP degree IN
plot(tableNumberLoopsForNeuron.HowManyLoops, tableNumberLoopsForNeuron.MAInDegree, '>', 'Color', [0.7412 0.7412 0.7412], 'MarkerFaceColor', [0.7412 0.7412 0.7412])
h = lsline; % plot regression line
mdl = fitlm(tableNumberLoopsForNeuron.HowManyLoops,tableNumberLoopsForNeuron.MAInDegree); % fit linear regression model
Rsquared = mdl.Rsquared.Ordinary;
RsquaredAdj = mdl.Rsquared.Adjusted;

% all NP degree
plot(tableNumberLoopsForNeuron.HowManyLoops, tableNumberLoopsForNeuron.MATotalDegree, 'ko', 'MarkerFaceColor', 'k')
h = lsline; % plot regression line
mdl = fitlm(tableNumberLoopsForNeuron.HowManyLoops,tableNumberLoopsForNeuron.MATotalDegree); % fit linear regression model
Rsquared = mdl.Rsquared.Ordinary;
RsquaredAdj = mdl.Rsquared.Adjusted;

h(3).LineWidth = 2;
h(3).LineStyle = '--';
h(2).LineWidth = 2;
h(2).LineStyle = '-.';
h(1).LineWidth = 2;

set(gcf,'color','w');
export_fig 'SupplFig_CorrelationWithMADegree.pdf' -painters