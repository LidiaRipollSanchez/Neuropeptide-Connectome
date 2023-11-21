%% Plot co-expression NPP - GPCRs
% Jan Watteyne - jan.watteyne@kuleuven.be or jwatteyne@gmail.com

% 20230607

clear all
close all
workspace

%% Import data
%-------------------------------------------
% import general variables from network build
%-------------------------------------------
% import gene ID
load ('05062023_Networks_build_EC500nM.mat','geneID');
% import neuron_ID
load ('05062023_Networks_build_EC500nM.mat','neuron_names');
% import anatomical_class
load ('05062023_Networks_build_EC500nM.mat','anatomy');

%----------------------------------
% import CeNGEN expression matrices
%----------------------------------
expression_threshold3 = readtable ('CeNGEN_threshold3.csv'); 
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
connections_EC50_1uM = connections_EC50(connections_EC50.EC50_M_ <=  1.0e-06, :); 
connections_EC50_500nM = connections_EC50(connections_EC50.EC50_M_ <=  500e-09, :); 
connections_EC50_100nM = connections_EC50(connections_EC50.EC50_M_ <=  100e-09, :); 
connections_EC50_50nM = connections_EC50(connections_EC50.EC50_M_ <=  50e-09, :); 
connections_EC50_10nM = connections_EC50(connections_EC50.EC50_M_ <=  10e-09, :); 

% each of these tables will have multiple entries for specific NPP-GPCR
% pairs (because multiple isoforms/peptides can be produced from one gene).
% Therefore, only retain unique combinations
connections_EC50_1uM = table(connections_EC50_1uM.GPCR_geneName, connections_EC50_1uM.NPP_geneName);
connections_EC50_1uM = unique(connections_EC50_1uM,'rows'); 

connections_EC50_500nM = table(connections_EC50_500nM.GPCR_geneName, connections_EC50_500nM.NPP_geneName);
connections_EC50_500nM = unique(connections_EC50_500nM,'rows'); 

connections_EC50_100nM = table(connections_EC50_100nM.GPCR_geneName, connections_EC50_100nM.NPP_geneName);
connections_EC50_100nM = unique(connections_EC50_100nM,'rows'); 

connections_EC50_50nM = table(connections_EC50_50nM.GPCR_geneName, connections_EC50_50nM.NPP_geneName);
connections_EC50_50nM = unique(connections_EC50_50nM,'rows'); 

connections_EC50_10nM = table(connections_EC50_10nM.GPCR_geneName, connections_EC50_10nM.NPP_geneName);
connections_EC50_10nM = unique(connections_EC50_10nM,'rows'); 
%% Convert expression data to binary expression matrix using functions Lidia Ripoll-Sanchez
[expression_threshold3_ordered] = Adapt_expression(expression_threshold3, neuron_names, geneID); %Adapt the expression values split by neuron type to neuron (from 180 columns to 302 columns) % J: column number doesn't change
[expression_threshold4_ordered] = Adapt_expression(expression_threshold4, neuron_names, geneID); %Adapt the expression values split by neuron type to neuron (from 180 columns to 302 columns) % J: column number doesn't change

[expression_threshold3_ordered_unweighted] = Mask_network(expression_threshold3_ordered); %Unweight expression values, when expression value is 1 otherwise is 0
[expression_threshold4_ordered_unweighted] = Mask_network(expression_threshold4_ordered); %Unweight expression values, when expression value is 1 otherwise is 0

% % Manual double check order of neurons to properly assign pharynx/sensory/interneuron/motoneuron classification
% NeuronNames_threshold3_ordered_unweighted = expression_threshold3_ordered_unweighted.Properties.VariableNames';
% NeuronNames_threshold3_ordered_unweighted(1) = [];
% 
% NeuronNames_threshold4_ordered_unweighted = expression_threshold4_ordered_unweighted.Properties.VariableNames';
% NeuronNames_threshold4_ordered_unweighted(1) = [];
% 
% NeuronNames_anatomicalclass = anatomy.Neuron;
% 
% % Compare this visually:
% dataTable = [NeuronNames_anatomicalclass, NeuronNames_threshold3_ordered_unweighted, NeuronNames_threshold4_ordered_unweighted];

%% Equalize number of NPP-GPCR connections 'rows' across EC50 thresholds
% Do this to be able to make co-expression matrices with the same
% dimensions, so these can be plotted on top of eachother

% 1: EC50 1 uM vs 500 nM
C = intersect(connections_EC50_1uM,connections_EC50_500nM,'rows'); % calculate intersect of both connection lists
index_noID = find(~ismember(connections_EC50_1uM,C,'rows'));

connections_EC50_500nM_Adj = connections_EC50_1uM;
for i = 1:length(index_noID)
    connections_EC50_500nM_Adj.Var1{index_noID(i)} = 'NaN';
    connections_EC50_500nM_Adj.Var2{index_noID(i)} = 'NaN';
end

% 2: EC50 1 uM vs 100 nM
C = intersect(connections_EC50_1uM,connections_EC50_100nM,'rows'); % calculate intersect of both connection lists
index_noID = find(~ismember(connections_EC50_1uM,C,'rows'));

connections_EC50_100nM_Adj = connections_EC50_1uM;
for i = 1:length(index_noID)
    connections_EC50_100nM_Adj.Var1{index_noID(i)} = 'NaN';
    connections_EC50_100nM_Adj.Var2{index_noID(i)} = 'NaN';
end

% 3: EC50 1 uM vs 50 nM
C = intersect(connections_EC50_1uM,connections_EC50_50nM,'rows'); % calculate intersect of both connection lists
index_noID = find(~ismember(connections_EC50_1uM,C,'rows'));

connections_EC50_50nM_Adj = connections_EC50_1uM;
for i = 1:length(index_noID)
    connections_EC50_50nM_Adj.Var1{index_noID(i)} = 'NaN';
    connections_EC50_50nM_Adj.Var2{index_noID(i)} = 'NaN';
end

% 4: EC50 1 uM vs 10 nM
C = intersect(connections_EC50_1uM,connections_EC50_10nM,'rows'); % calculate intersect of both connection lists
index_noID = find(~ismember(connections_EC50_1uM,C,'rows'));

connections_EC50_10nM_Adj = connections_EC50_1uM;
for i = 1:length(index_noID)
    connections_EC50_10nM_Adj.Var1{index_noID(i)} = 'NaN';
    connections_EC50_10nM_Adj.Var2{index_noID(i)} = 'NaN';
end

clear i
clear index_noID
clear C

%% Compute matrices
% compute binary matrices for either NPP or GPCR that show whether a single
% gene (y-axis) is expressed in a particular neuron (x-axis). Do this for
% different CeNGEN and EC50 thresholds

numRows = size(connections_EC50_1uM, 1);
numCols = size(expression_threshold3_ordered_unweighted, 2)-1;

% initiate empty matrices to fill
expressionNPP_threshold4_1uM = nan(numRows, numCols);
expressionNPP_threshold4_500nM = expressionNPP_threshold4_1uM;
expressionNPP_threshold4_100nM = expressionNPP_threshold4_1uM;
expressionNPP_threshold4_50nM = expressionNPP_threshold4_1uM;
expressionNPP_threshold4_10nM = expressionNPP_threshold4_1uM;

expressionNPP_threshold3_1uM = expressionNPP_threshold4_1uM;
expressionNPP_threshold3_500nM = expressionNPP_threshold4_1uM;
expressionNPP_threshold3_100nM =expressionNPP_threshold4_1uM;
expressionNPP_threshold3_50nM =expressionNPP_threshold4_1uM;
expressionNPP_threshold3_10nM = expressionNPP_threshold4_1uM;

expressionGPCR_threshold4_1uM = expressionNPP_threshold4_1uM;
expressionGPCR_threshold4_500nM = expressionNPP_threshold4_1uM;
expressionGPCR_threshold4_100nM = expressionNPP_threshold4_1uM;
expressionGPCR_threshold4_50nM = expressionNPP_threshold4_1uM;
expressionGPCR_threshold4_10nM = expressionNPP_threshold4_1uM;

expressionGPCR_threshold3_1uM = expressionNPP_threshold4_1uM;
expressionGPCR_threshold3_500nM = expressionNPP_threshold4_1uM;
expressionGPCR_threshold3_100nM =expressionNPP_threshold4_1uM;
expressionGPCR_threshold3_50nM =expressionNPP_threshold4_1uM;
expressionGPCR_threshold3_10nM = expressionNPP_threshold4_1uM;

%-------------------------------
% CeNGEN threshold 3 - EC50 1 uM
%-------------------------------
for row = 1:numRows
    nameNPP = lower(connections_EC50_1uM.Var2{row}); % get name of NPP, and put it to lower case
    [~,Index] = intersect(string(expression_threshold3_ordered_unweighted.Gene), nameNPP);
    expressionNPP_threshold3_1uM(row, :) = table2array(expression_threshold3_ordered_unweighted(Index,2:end));
 
    nameGPCR = connections_EC50_1uM.Var1{row}; % get name of GPCR
    [~,Index] = intersect(string(expression_threshold3_ordered_unweighted.Gene), nameGPCR);
    expressionGPCR_threshold3_1uM(row, :) = table2array(expression_threshold3_ordered_unweighted(Index,2:end));
end
 
%---------------------------------
% CeNGEN threshold 3 - EC50 500 nM
%---------------------------------
for row = 1:numRows
    try
    nameNPP = lower(connections_EC50_500nM_Adj.Var2{row}); % get name of NPP, and put it to lower case
    [~,Index] = intersect(string(expression_threshold3_ordered_unweighted.Gene), nameNPP);
    expressionNPP_threshold3_500nM(row, :) = table2array(expression_threshold3_ordered_unweighted(Index,2:end));
 
    nameGPCR = connections_EC50_500nM_Adj.Var1{row}; % get name of GPCR
    [~,Index] = intersect(string(expression_threshold3_ordered_unweighted.Gene), nameGPCR);
    expressionGPCR_threshold3_500nM(row, :) = table2array(expression_threshold3_ordered_unweighted(Index,2:end));
    catch
        expressionNPP_threshold3_500nM(row, :) = zeros(1, numCols);
        expressionGPCR_threshold3_500nM(row, :) = zeros(1, numCols);
    end
end

%---------------------------------
% CeNGEN threshold 3 - EC50 100 nM
%---------------------------------
for row = 1:numRows
    try
    nameNPP = lower(connections_EC50_100nM_Adj.Var2{row}); % get name of NPP, and put it to lower case
    [~,Index] = intersect(string(expression_threshold3_ordered_unweighted.Gene), nameNPP);
    expressionNPP_threshold3_100nM(row, :) = table2array(expression_threshold3_ordered_unweighted(Index,2:end));
 
    nameGPCR = connections_EC50_100nM_Adj.Var1{row}; % get name of GPCR
    [~,Index] = intersect(string(expression_threshold3_ordered_unweighted.Gene), nameGPCR);
    expressionGPCR_threshold3_100nM(row, :) = table2array(expression_threshold3_ordered_unweighted(Index,2:end));
    catch
        expressionNPP_threshold3_100nM(row, :) = zeros(1, numCols);
        expressionGPCR_threshold3_100nM(row, :) = zeros(1, numCols);
    end
end

%---------------------------------
% CeNGEN threshold 3 - EC50 100 nM
%---------------------------------
for row = 1:numRows
   try
    nameNPP = lower(connections_EC50_50nM_Adj.Var2{row}); % get name of NPP, and put it to lower case
    [~,Index] = intersect(string(expression_threshold3_ordered_unweighted.Gene), nameNPP);
    expressionNPP_threshold3_50nM(row, :) = table2array(expression_threshold3_ordered_unweighted(Index,2:end));
 
    nameGPCR = connections_EC50_50nM_Adj.Var1{row}; % get name of GPCR
    [~,Index] = intersect(string(expression_threshold3_ordered_unweighted.Gene), nameGPCR);
    expressionGPCR_threshold3_50nM(row, :) = table2array(expression_threshold3_ordered_unweighted(Index,2:end));
    catch
        expressionNPP_threshold3_50nM(row, :) = zeros(1, numCols);
        expressionGPCR_threshold3_50nM(row, :) = zeros(1, numCols);
    end
end

%---------------------------------
% CeNGEN threshold 3 - EC50 10 nM
%---------------------------------
for row = 1:numRows
    try
    nameNPP = lower(connections_EC50_10nM_Adj.Var2{row}); % get name of NPP, and put it to lower case
    [~,Index] = intersect(string(expression_threshold3_ordered_unweighted.Gene), nameNPP);
    expressionNPP_threshold3_10nM(row, :) = table2array(expression_threshold3_ordered_unweighted(Index,2:end));
 
    nameGPCR = connections_EC50_10nM_Adj.Var1{row}; % get name of GPCR
    [~,Index] = intersect(string(expression_threshold3_ordered_unweighted.Gene), nameGPCR);
    expressionGPCR_threshold3_10nM(row, :) = table2array(expression_threshold3_ordered_unweighted(Index,2:end));
    catch
        expressionNPP_threshold3_10nM(row, :) = zeros(1, numCols);
        expressionGPCR_threshold3_10nM(row, :) = zeros(1, numCols);
    end
end

%-------------------------------
% CeNGEN threshold 4 - EC50 1 uM
%-------------------------------
for row = 1:numRows
    nameNPP = lower(connections_EC50_1uM.Var2{row}); % get name of NPP, and put it to lower case
    [~,Index] = intersect(string(expression_threshold4_ordered_unweighted.Gene), nameNPP);
    expressionNPP_threshold4_1uM(row, :) = table2array(expression_threshold4_ordered_unweighted(Index,2:end));
 
    nameGPCR = connections_EC50_1uM.Var1{row}; % get name of GPCR
    [~,Index] = intersect(string(expression_threshold4_ordered_unweighted.Gene), nameGPCR);
    expressionGPCR_threshold4_1uM(row, :) = table2array(expression_threshold4_ordered_unweighted(Index,2:end));
end
 
%---------------------------------
% CeNGEN threshold 4 - EC50 500 nM
%---------------------------------
for row = 1:numRows
    try
    nameNPP = lower(connections_EC50_500nM_Adj.Var2{row}); % get name of NPP, and put it to lower case
    [~,Index] = intersect(string(expression_threshold4_ordered_unweighted.Gene), nameNPP);
    expressionNPP_threshold4_500nM(row, :) = table2array(expression_threshold4_ordered_unweighted(Index,2:end));
 
    nameGPCR = connections_EC50_500nM_Adj.Var1{row}; % get name of GPCR
    [~,Index] = intersect(string(expression_threshold4_ordered_unweighted.Gene), nameGPCR);
    expressionGPCR_threshold4_500nM(row, :) = table2array(expression_threshold4_ordered_unweighted(Index,2:end));
    catch
        expressionNPP_threshold4_500nM(row, :) = zeros(1, numCols);
        expressionGPCR_threshold4_500nM(row, :) = zeros(1, numCols);
    end
end
 
%---------------------------------
% CeNGEN threshold 4 - EC50 100 nM
%---------------------------------
for row = 1:numRows
    try
    nameNPP = lower(connections_EC50_100nM_Adj.Var2{row}); % get name of NPP, and put it to lower case
    [~,Index] = intersect(string(expression_threshold4_ordered_unweighted.Gene), nameNPP);
    expressionNPP_threshold4_100nM(row, :) = table2array(expression_threshold4_ordered_unweighted(Index,2:end));
 
    nameGPCR = connections_EC50_100nM_Adj.Var1{row}; % get name of GPCR
    [~,Index] = intersect(string(expression_threshold4_ordered_unweighted.Gene), nameGPCR);
    expressionGPCR_threshold4_100nM(row, :) = table2array(expression_threshold4_ordered_unweighted(Index,2:end));
    catch
        expressionNPP_threshold4_100nM(row, :) = zeros(1, numCols);
        expressionGPCR_threshold4_100nM(row, :) = zeros(1, numCols);
    end
end

%---------------------------------
% CeNGEN threshold 4 - EC50 50 nM
%---------------------------------
for row = 1:numRows
    try
    nameNPP = lower(connections_EC50_50nM_Adj.Var2{row}); % get name of NPP, and put it to lower case
    [~,Index] = intersect(string(expression_threshold4_ordered_unweighted.Gene), nameNPP);
    expressionNPP_threshold4_50nM(row, :) = table2array(expression_threshold4_ordered_unweighted(Index,2:end));
 
    nameGPCR = connections_EC50_50nM_Adj.Var1{row}; % get name of GPCR
    [~,Index] = intersect(string(expression_threshold4_ordered_unweighted.Gene), nameGPCR);
    expressionGPCR_threshold4_50nM(row, :) = table2array(expression_threshold4_ordered_unweighted(Index,2:end));
    catch
        expressionNPP_threshold4_50nM(row, :) = zeros(1, numCols);
        expressionGPCR_threshold4_50nM(row, :) = zeros(1, numCols);
    end
end

%---------------------------------
% CeNGEN threshold 4 - EC50 10 nM
%---------------------------------
for row = 1:numRows
    try
    nameNPP = lower(connections_EC50_10nM_Adj.Var2{row}); % get name of NPP, and put it to lower case
    [~,Index] = intersect(string(expression_threshold4_ordered_unweighted.Gene), nameNPP);
    expressionNPP_threshold4_10nM(row, :) = table2array(expression_threshold4_ordered_unweighted(Index,2:end));
 
    nameGPCR = connections_EC50_10nM_Adj.Var1{row}; % get name of GPCR
    [~,Index] = intersect(string(expression_threshold4_ordered_unweighted.Gene), nameGPCR);
    expressionGPCR_threshold4_10nM(row, :) = table2array(expression_threshold4_ordered_unweighted(Index,2:end));
    catch
        expressionNPP_threshold4_10nM(row, :) = zeros(1, numCols);
        expressionGPCR_threshold4_10nM(row, :) = zeros(1, numCols);
    end
end

clear numCols
clear numRows
clear nameGPCR
clear nameNPP
clear row
clear Index

%% plot scatter for NPP and GPCR expression - MAIN FIGURE EC50 500 nM and CENGEN thrsh 4
expressionNPP = expressionNPP_threshold4_500nM;
expressionGPCR = expressionGPCR_threshold4_500nM;
coexpression = expressionNPP + expressionGPCR;
coexpression(coexpression ~= 2) = 0;

figure('Renderer', 'painters', 'Position', [10 10 2500 300])
% first, plot gray dots for expression in the neuron class
[y, x] = find(expressionNPP); %get x and y coordinates to plot scatter
scatter(x, y, 15, 'MarkerFaceColor', [0.7    0.7    0.7], 'MarkerEdgeColor', 'none'); %plot scatter of long range (scatter instead of spy to choose colour)
xlim([0 size(expressionNPP, 2)+1])
ylim([0 size(expressionNPP, 1)+1])
hold on 

% next, highlight dots for autocrine connections
[y, x] = find(coexpression == 2); %get x and y coordinates to plot scatter
scatter(x, y, 20, 'MarkerFaceColor', [0.1    0.1    0.1], 'MarkerEdgeColor', 'none'); %plot scatter of long range (scatter instead of spy to choose colour)

axis off

% add lines for each class: pharynx - sensory - interneuron - motoneuron
line([20.5 20.5],[0 size(expressionNPP, 1)+1],'Color',[0 0 0])
line([99.5 99.5],[0 size(expressionNPP, 1)+1],'Color',[0 0 0])
line([181.5 181.5],[0 size(expressionNPP, 1)+1],'Color',[0 0 0])
line([0 0],[0 size(expressionNPP, 1)+1],'Color',[0 0 0])
line([303 303],[0 size(expressionNPP, 1)+1],'Color',[0 0 0])
line([0 303],[0 0],'Color',[0 0 0])
line([0 303],[size(expressionNPP, 1)+1 size(expressionNPP, 1)+1],'Color',[0 0 0])

set(gcf,'color','w');
%export_fig 'MainFig_1A_ExpressionNPP_with_autocrine_connections_threshold4_500nm.pdf' -painters

% plot scatter for GPCR expression
figure('Renderer', 'painters', 'Position', [10 10 2500 300])
[y, x] = find(expressionGPCR); %get x and y coordinates to plot scatter
s1 = scatter(x, y, 15, 'MarkerFaceColor', [0.7    0.7    0.7], 'MarkerEdgeColor', 'none'); %plot scatter of long range (scatter instead of spy to choose colour)
xlim([0 size(expressionGPCR, 2)])
ylim([0 size(expressionGPCR, 1)+1])
hold on 

% next, highlight dots for autocrine connections
[y, x] = find(coexpression == 2); %get x and y coordinates to plot scatter
s1 = scatter(x, y, 20, 'MarkerFaceColor', [0.1    0.1    0.1], 'MarkerEdgeColor', 'none'); %plot scatter of long range (scatter instead of spy to choose colour)
axis off

% add lines for each class: pharynx - sensory - interneuron - motoneuron
line([20.5 20.5],[0 size(expressionNPP, 1)+1],'Color',[0 0 0])
line([99.5 99.5],[0 size(expressionNPP, 1)+1],'Color',[0 0 0])
line([181.5 181.5],[0 size(expressionNPP, 1)+1],'Color',[0 0 0])
line([0 0],[0 size(expressionNPP, 1)+1],'Color',[0 0 0])
line([303 303],[0 size(expressionNPP, 1)+1],'Color',[0 0 0])
line([0 303],[0 0],'Color',[0 0 0])
line([0 303],[size(expressionNPP, 1)+1 size(expressionNPP, 1)+1],'Color',[0 0 0])

set(gcf,'color','w');
%export_fig 'MainFig_1A_ExpressionGPCR_with_autocrine_connections_threshold4_500nm.pdf' -painters

%% Compute and barplot the fraction of neurons showing autocrine connections - MAIN FIGURE THRESHOLD 4 - EC50 500nM
bardata = nan(1,5); %initiate empty vector

CoExpressionMatrix = expressionGPCR_threshold4_500nM + expressionNPP_threshold4_500nM;
CoExpressionMatrix(CoExpressionMatrix ~= 2) = 0;

neuronWithOrWithoutCoExpression = sum(CoExpressionMatrix, 1);
numberNeuronsWithCoExpression = length(find(neuronWithOrWithoutCoExpression ~= 0));
bardata(1) = numberNeuronsWithCoExpression / (size(CoExpressionMatrix, 2)-1);

indexPharynx = find(contains(anatomy.FinalClassification,'Pharynx'));
neuronWithOrWithoutCoExpression_Pharynx = neuronWithOrWithoutCoExpression(indexPharynx);
bardata(2) = length(find(neuronWithOrWithoutCoExpression_Pharynx ~= 0)) / (length(neuronWithOrWithoutCoExpression_Pharynx));

indexSensory = find(contains(anatomy.FinalClassification,'sensory neuron'));
neuronWithOrWithoutCoExpression_Sensory = neuronWithOrWithoutCoExpression(indexSensory);
bardata(3) = length(find(neuronWithOrWithoutCoExpression_Sensory ~= 0)) / (length(neuronWithOrWithoutCoExpression_Sensory));

indexInter = find(contains(anatomy.FinalClassification,'interneuron'));
neuronWithOrWithoutCoExpression_Inter = neuronWithOrWithoutCoExpression(indexInter);
bardata(4) = length(find(neuronWithOrWithoutCoExpression_Inter ~= 0)) / (length(neuronWithOrWithoutCoExpression_Inter));

indexMotor = find(contains(anatomy.FinalClassification,'motor neuron'));
neuronWithOrWithoutCoExpression_Motor = neuronWithOrWithoutCoExpression(indexMotor);
bardata(5) = length(find(neuronWithOrWithoutCoExpression_Motor ~= 0)) / (length(neuronWithOrWithoutCoExpression_Motor));


figure('Renderer', 'painters', 'Position', [10 10 578 209])
b = bar(bardata);
b.FaceColor = 'flat';
%b.CData = [0.5    0.5   0.5];
yticks(0:0.25:1)
yticklabels('')
xticklabels('')
ylim([0 1])
xlim([0.5 5.5])

set(gcf,'color','w');
%export_fig 'MainFig_B_Fraction_Neurons_with_Autocrine_Connections_threshold4_500nM.pdf' -painters



%% Plot co-expression matrix for each of EC50 thresholds - THRESHOLD 4 - SUPPL FIG

CoExpressionMatrix_threshold4_1uM = expressionNPP_threshold4_1uM + expressionGPCR_threshold4_1uM;
CoExpressionMatrix_threshold4_1uM(CoExpressionMatrix_threshold4_1uM ~= 2) = 0;
CoExpressionMatrix_threshold4_1uM(CoExpressionMatrix_threshold4_1uM == 2) = 1;

CoExpressionMatrix_threshold4_500nM = expressionNPP_threshold4_500nM + expressionGPCR_threshold4_500nM;
CoExpressionMatrix_threshold4_500nM(CoExpressionMatrix_threshold4_500nM ~= 2) = 0;
CoExpressionMatrix_threshold4_500nM(CoExpressionMatrix_threshold4_500nM == 2) = 1;

CoExpressionMatrix_threshold4_100nM = expressionNPP_threshold4_100nM + expressionGPCR_threshold4_100nM;
CoExpressionMatrix_threshold4_100nM(CoExpressionMatrix_threshold4_100nM ~= 2) = 0;
CoExpressionMatrix_threshold4_100nM(CoExpressionMatrix_threshold4_100nM == 2) = 1;

CoExpressionMatrix_threshold4_50nM = expressionNPP_threshold4_50nM + expressionGPCR_threshold4_50nM;
CoExpressionMatrix_threshold4_50nM(CoExpressionMatrix_threshold4_50nM ~= 2) = 0;
CoExpressionMatrix_threshold4_50nM(CoExpressionMatrix_threshold4_50nM == 2) = 1;

CoExpressionMatrix_threshold4_10nM = expressionNPP_threshold4_10nM + expressionGPCR_threshold4_10nM;
CoExpressionMatrix_threshold4_10nM(CoExpressionMatrix_threshold4_10nM ~= 2) = 0;
CoExpressionMatrix_threshold4_10nM(CoExpressionMatrix_threshold4_10nM == 2) = 1;

figure('Renderer', 'painters', 'Position', [10 10 2500 300])
[y, x] = find(CoExpressionMatrix_threshold4_1uM); %get x and y coordinates to plot scatter
s1 = scatter(x, y, 15, 'MarkerFaceColor', [1.0000    1.0000    0.6980], 'MarkerEdgeColor', 'none');
xlim([0 size(CoExpressionMatrix_threshold4_1uM, 2)+1])
ylim([0 size(CoExpressionMatrix_threshold4_1uM, 1)+1])
hold on
axis off

[y, x] = find(CoExpressionMatrix_threshold4_500nM); %get x and y coordinates to plot scatter
s1 = scatter(x, y, 15, 'MarkerFaceColor', [0.9961    0.8000    0.3608], 'MarkerEdgeColor', 'none');

[y, x] = find(CoExpressionMatrix_threshold4_100nM); %get x and y coordinates to plot scatter
s1 = scatter(x, y, 15, 'MarkerFaceColor', [0.9922    0.5529    0.2353], 'MarkerEdgeColor', 'none');

[y, x] = find(CoExpressionMatrix_threshold4_50nM); %get x and y coordinates to plot scatter
s1 = scatter(x, y, 15, 'MarkerFaceColor', [0.9412    0.2314    0.1255], 'MarkerEdgeColor', 'none');

[y, x] = find(CoExpressionMatrix_threshold4_10nM); %get x and y coordinates to plot scatter
s1 = scatter(x, y, 15, 'MarkerFaceColor', [0.7412         0    0.1490], 'MarkerEdgeColor', 'none');

% add lines for each class: pharynx - sensory - interneuron - motoneuron
line([20.5 20.5],[0 size(CoExpressionMatrix_threshold4_1uM, 1)+1],'Color',[0 0 0])
line([99.5 99.5],[0 size(CoExpressionMatrix_threshold4_1uM, 1)+1],'Color',[0 0 0])
line([181.5 181.5],[0 size(CoExpressionMatrix_threshold4_1uM, 1)+1],'Color',[0 0 0])
line([0 0],[0 size(CoExpressionMatrix_threshold4_1uM, 1)+1],'Color',[0 0 0])
line([303 303],[0 size(CoExpressionMatrix_threshold4_1uM, 1)+1],'Color',[0 0 0])
line([0 303],[0 0],'Color',[0 0 0])
line([0 303],[size(CoExpressionMatrix_threshold4_1uM, 1)+1 size(CoExpressionMatrix_threshold4_1uM, 1)+1],'Color',[0 0 0])

set(gcf,'color','w');
%export_fig 'Suppl_CoExpressionMatrix_AllEC50Thresholds_CENGEN4.pdf' -painters

%% Plot co-expression matrix for each of EC50 thresholds - THRESHOLD 3 - SUPPL FIG

CoExpressionMatrix_threshold3_1uM = expressionNPP_threshold3_1uM + expressionGPCR_threshold3_1uM;
CoExpressionMatrix_threshold3_1uM(CoExpressionMatrix_threshold3_1uM ~= 2) = 0;
CoExpressionMatrix_threshold3_1uM(CoExpressionMatrix_threshold3_1uM == 2) = 1;
 
CoExpressionMatrix_threshold3_500nM = expressionNPP_threshold3_500nM + expressionGPCR_threshold3_500nM;
CoExpressionMatrix_threshold3_500nM(CoExpressionMatrix_threshold3_500nM ~= 2) = 0;
CoExpressionMatrix_threshold3_500nM(CoExpressionMatrix_threshold3_500nM == 2) = 1;
 
CoExpressionMatrix_threshold3_100nM = expressionNPP_threshold3_100nM + expressionGPCR_threshold3_100nM;
CoExpressionMatrix_threshold3_100nM(CoExpressionMatrix_threshold3_100nM ~= 2) = 0;
CoExpressionMatrix_threshold3_100nM(CoExpressionMatrix_threshold3_100nM == 2) = 1;
 
CoExpressionMatrix_threshold3_50nM = expressionNPP_threshold3_50nM + expressionGPCR_threshold3_50nM;
CoExpressionMatrix_threshold3_50nM(CoExpressionMatrix_threshold3_50nM ~= 2) = 0;
CoExpressionMatrix_threshold3_50nM(CoExpressionMatrix_threshold3_50nM == 2) = 1;
 
CoExpressionMatrix_threshold3_10nM = expressionNPP_threshold3_10nM + expressionGPCR_threshold3_10nM;
CoExpressionMatrix_threshold3_10nM(CoExpressionMatrix_threshold3_10nM ~= 2) = 0;
CoExpressionMatrix_threshold3_10nM(CoExpressionMatrix_threshold3_10nM == 2) = 1;

figure('Renderer', 'painters', 'Position', [10 10 2500 300])
[y, x] = find(CoExpressionMatrix_threshold3_1uM); %get x and y coordinates to plot scatter
s1 = scatter(x, y, 15, 'MarkerFaceColor', [1.0000    1.0000    0.8000], 'MarkerEdgeColor', 'none');
xlim([0 size(CoExpressionMatrix_threshold3_1uM, 2)+1])
ylim([0 size(CoExpressionMatrix_threshold3_1uM, 1)+1])
hold on
axis off

[y, x] = find(CoExpressionMatrix_threshold3_500nM); %get x and y coordinates to plot scatter
s1 = scatter(x, y, 15, 'MarkerFaceColor', [0.6314    0.8549    0.7059], 'MarkerEdgeColor', 'none');

[y, x] = find(CoExpressionMatrix_threshold3_100nM); %get x and y coordinates to plot scatter
s1 = scatter(x, y, 15, 'MarkerFaceColor', [0.2549    0.7137    0.7686], 'MarkerEdgeColor', 'none');

[y, x] = find(CoExpressionMatrix_threshold3_50nM); %get x and y coordinates to plot scatter
s1 = scatter(x, y, 15, 'MarkerFaceColor', [0.1725    0.4980    0.7216], 'MarkerEdgeColor', 'none');

[y, x] = find(CoExpressionMatrix_threshold3_10nM); %get x and y coordinates to plot scatter
s1 = scatter(x, y, 15, 'MarkerFaceColor', [0.1451    0.2039    0.5804], 'MarkerEdgeColor', 'none');

% add lines for each class: pharynx - sensory - interneuron - motoneuron
line([20.5 20.5],[0 size(CoExpressionMatrix_threshold4_1uM, 1)+1],'Color',[0 0 0])
line([99.5 99.5],[0 size(CoExpressionMatrix_threshold4_1uM, 1)+1],'Color',[0 0 0])
line([181.5 181.5],[0 size(CoExpressionMatrix_threshold4_1uM, 1)+1],'Color',[0 0 0])
line([0 0],[0 size(CoExpressionMatrix_threshold4_1uM, 1)+1],'Color',[0 0 0])
line([303 303],[0 size(CoExpressionMatrix_threshold4_1uM, 1)+1],'Color',[0 0 0])
line([0 303],[0 0],'Color',[0 0 0])
line([0 303],[size(CoExpressionMatrix_threshold4_1uM, 1)+1 size(CoExpressionMatrix_threshold4_1uM, 1)+1],'Color',[0 0 0])

set(gcf,'color','w');
%export_fig 'Suppl_CoExpressionMatrix_AllEC50Thresholds_CENGEN3.pdf' -painters
%% Look at distribution of Co-expression across neuron classes
CoExpressionMatrix = CoExpressionMatrix_threshold4_500nM;
col = [0.1451    0.2039    0.5804];
%1.0000    0.9294    0.6275
CoExpressionMatrix(CoExpressionMatrix == 2) = 1;

numberAutocrineConnections = sum(CoExpressionMatrix,1)';
numberAutocrineConnections(numberAutocrineConnections == 0) = NaN;

NeuronAC = [neuron_names, num2cell(numberAutocrineConnections), anatomy.FinalClassification];

IndexPharynx = find(count(string(NeuronAC(:,3)), 'Pharynx'));
IndexSensory = find(count(string(NeuronAC(:,3)), 'sensory neuron'));
IndexInter = find(count(string(NeuronAC(:,3)), 'interneuron'));
IndexMotor = find(count(string(NeuronAC(:,3)), 'motor neuron'));

%figure('Renderer', 'painters', 'Position', [10 10 400 300])
figure('Renderer', 'painters', 'Position', [10 10 380 380])
% first, plot histogram for all neuron classes
numbers = cell2mat(NeuronAC(:,2));
maxCounts = max(numbers);
edges = [1:maxCounts+1];
counts = histcounts(numbers, edges);
counts = counts/302; % normalize to fraction of all neurons
plot(1:maxCounts, counts, 'LineStyle', '-', 'color', col, 'LineWidth', 2) 
hold on

% then, plot each different neuron class
numbers = cell2mat(NeuronAC(IndexPharynx,2));
counts = histcounts(numbers, edges);
counts = counts/302; % normalize to fraction of all neurons
plot(1:maxCounts, counts, 'LineStyle', '-.', 'color', col, 'LineWidth', 1) 

numbers = cell2mat(NeuronAC(IndexSensory,2));
counts = histcounts(numbers, edges);
counts = counts/302; % normalize to fraction of all neurons
plot(1:maxCounts, counts, 'LineStyle', '-', 'color', col, 'LineWidth', 1)  

numbers = cell2mat(NeuronAC(IndexInter,2));
counts = histcounts(numbers, edges);
counts = counts/302; % normalize to fraction of all neurons
plot(1:maxCounts, counts, 'LineStyle', '--', 'color', col, 'LineWidth', 1) 

numbers = cell2mat(NeuronAC(IndexMotor,2));
counts = histcounts(numbers, edges);
counts = counts/302; % normalize to fraction of all neurons
plot(1:maxCounts, counts, 'LineStyle', ':', 'color', col, 'LineWidth', 1) 

ylim([0 0.2])
yticks(0:0.05:0.2)
xlim([0 20])
xticks(0:5:20)
% xlim([0 14])
% xticks(0:2:14)
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gcf,'color','w');
%export_fig 'Suppl_Distribution_coexpression_threshold3_10nM.pdf' -painters

