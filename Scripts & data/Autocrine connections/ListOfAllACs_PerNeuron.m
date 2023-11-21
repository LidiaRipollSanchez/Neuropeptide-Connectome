%% Make list Autocrine connections for chosen neurons
% Jan Watteyne
% 20230608

clear all
close all
workspace
%% Import data
% import neuron_names
load ('05062023_Networks_build_EC500nM.mat','neuron_names');
% import gene ID
load ('05062023_Networks_build_EC500nM.mat','geneID');
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

%% Compute matrices
% compute binary matrices for either NPP or GPCR that show whether a single
% gene (y-axis) is expressed in a particular neuron (x-axis). Do this for
% different CeNGEN and EC50 thresholds

numRows = size(connections_EC50_500nM, 1);
numCols = size(expression_threshold4_ordered_unweighted, 2)-1;

% initiate empty matrices to fill
expressionNPP = nan(numRows, numCols);
expressionGPCR = expressionNPP;

%-------------------------------
% CeNGEN threshold 4 - EC50 500nM
%-------------------------------
for row = 1:numRows
    nameNPP = lower(connections_EC50_500nM.Var2{row}); % get name of NPP, and put it to lower case
    [~,Index] = intersect(string(expression_threshold4_ordered_unweighted.Gene), nameNPP);
    expressionNPP(row, :) = table2array(expression_threshold4_ordered_unweighted(Index,2:end));
 
    nameGPCR = connections_EC50_500nM.Var1{row}; % get name of GPCR
    [~,Index] = intersect(string(expression_threshold4_ordered_unweighted.Gene), nameGPCR);
    expressionGPCR(row, :) = table2array(expression_threshold4_ordered_unweighted(Index,2:end));
end
 
clear numCols
clear numRows
clear nameGPCR
clear nameNPP
clear row
clear Index


CoExpressionMatrix = expressionGPCR + expressionNPP;
CoExpressionMatrix(CoExpressionMatrix ~= 2) = 0;
CoExpressionMatrix(CoExpressionMatrix == 2) = 1;

numberAutocrineConnections = sum(CoExpressionMatrix,1)';
numberAutocrineConnections(numberAutocrineConnections == 0) = NaN;

NeuronName = neuron_names;
NeuronAC = [NeuronName, num2cell(numberAutocrineConnections)];

%% Make big cell with all AC and the host neurons

Cell_AC_Neuron = [];
for i = 1:length(NeuronName)
   Name = connections_EC50_500nM{i,:};
   NeuronCell = NeuronName{i};
   Cell_AC_Neuron = [Cell_AC_Neuron; repmat(Name, size(NeuronCell, 1), 1), NeuronCell];
end
