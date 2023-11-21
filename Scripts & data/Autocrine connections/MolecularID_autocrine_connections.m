%% Script to analyze the molecular identity of the autocrine connections 
% Jan Watteyne - jan.watteyne@kuleuven.be and/or jwatteyne@gmail.com
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
    NPP = table2array(expression_threshold4_ordered_unweighted(strcmpi(expression_threshold4_ordered_unweighted.(1), char(connections_EC50_500nM.Var1(i))), 2:303)); %get the vector for each neuropeptide
    GPCR = table2array(expression_threshold4_ordered_unweighted(strcmpi(expression_threshold4_ordered_unweighted.(1), char(connections_EC50_500nM.Var2(i))), 2:303)); %get the vector for each receptor
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

%% Rank and plot the NP - GPCR couples based on the number of times they figure in a autocrine connection
%make a table with column 1: NP - GPCR connections and column 2: number of times it is used in an autocrine connection
tableTimesInLoop = connections_EC50_500nM;
tableTimesInLoop.Properties.VariableNames = [{'GeneGPCR','GeneLigand'}];

tableTimesInLoop.TimesInLoop = cell2mat(self_loops_per_pair)';
tableTimesInLoop = sortrows(tableTimesInLoop, 3); % order the table in ascending order (how many autocrine loops)


%make an extra column to indicate if the GPCR is promiscuous or not
V = table2cell(tableTimesInLoop(:, 'GeneGPCR'));
[~, ~, idy] = unique(V); % get list of receptors that form a NPP connection and get index to count how many ligands each receptor binds
U = unique(idy);
cnt = histc(idy,U);
tableTimesInLoop.promiscuousGPCR = ismember(idy,find(cnt>1)); % locations of duplicates


%make extra columns to indicate whether the connection is flp or nlp type
IndexFLP = find(count(string(tableTimesInLoop.GeneLigand), 'flp')); % find the index of all entries in which the string 'FLP' is found
tableTimesInLoop.FLPtype = repelem(0, size(tableTimesInLoop, 1))';
tableTimesInLoop.FLPtype(IndexFLP) = 1;

IndexNLP = find(count(string(tableTimesInLoop.GeneLigand), 'nlp')); % find the index of all entries in which the string 'NLP' is found
tableTimesInLoop.NLPtype = repelem(0, size(tableTimesInLoop, 1))';
tableTimesInLoop.NLPtype(IndexNLP) = 1;

% some peptides (SNET-1 and PDF-1) are categorized as NLP
IdxFLPandNLP = sort([IndexFLP; IndexNLP]); % Some peptides (like SNET-1) do not have the FLP/NLP nomenclature
IndexOther = setdiff([1:size(tableTimesInLoop, 1)], IdxFLPandNLP);
tableTimesInLoop.NLPtype(IndexOther) = 1;

%% Plot the number of Autocrine Connections for each NPP-GPCR pair
figure('Renderer', 'painters', 'Position', [10 10 800 400])

SizeMarker = 20;
ColorFLP = [0.5529 0.8275 0.7804]; % with colorbrewer https://colorbrewer2.org/
ColorNLP = [0.7451 0.7294 0.8549];

% scatter FLP type
index = find(tableTimesInLoop.FLPtype);
scatter(index,tableTimesInLoop.TimesInLoop(index),SizeMarker,'MarkerEdgeColor',ColorFLP,...
              'MarkerFaceColor',ColorFLP,...
              'LineWidth',1.5), hold on

% scatter NLP type
index = find(tableTimesInLoop.NLPtype);
scatter(index,tableTimesInLoop.TimesInLoop(index),SizeMarker,'MarkerEdgeColor',ColorNLP,...
              'MarkerFaceColor',ColorNLP,...
              'LineWidth',1.5), 
          
          
% scatter promiscous GPCR
index = find(tableTimesInLoop.promiscuousGPCR);
scatter(index, tableTimesInLoop.TimesInLoop(index), SizeMarker,'MarkerEdgeColor',[0 0 0],...
              'LineWidth',1), 


xlim([1 height(tableTimesInLoop)])
%ylabel('Number of times in autocrine connection')
yticklabels('')
xticks('')
ylim([0 60])
yticks(0:20:60)

set(gcf,'color','w');
%export_fig 'MainFig_C_NumberLoopsForEachPair.pdf' -painters

%% Rank each individual GPCR based on its frequency in autocrine connections

V = table2cell(tableTimesInLoop(:, 'GeneGPCR'));
idx = unique(V); % get list of all unique receptors

tableGPCR_TimesInLoop = table; % initiate empty table
tableGPCR_TimesInLoop.GeneGPCR = idx;

for ii = 1:length(idx) %iterate over each GPCR and sum the number of times it figures in a loop
    NameGPCR = idx(ii);
    index = find(string(tableTimesInLoop.GeneGPCR) == NameGPCR);
    tableGPCR_TimesInLoop.TimesInLoop(ii) = sum(tableTimesInLoop.TimesInLoop(index));
    tableGPCR_TimesInLoop.promiscuousGPCR(ii) = tableTimesInLoop.promiscuousGPCR(index(1));
    tableGPCR_TimesInLoop.FLPtype(ii) = tableTimesInLoop.FLPtype(index(1));
    tableGPCR_TimesInLoop.NLPtype(ii) = tableTimesInLoop.NLPtype(index(1));
end

tableGPCR_TimesInLoop = sortrows(tableGPCR_TimesInLoop, 2); % order the table in ascending order (how many autocrine loops)
figure('Renderer', 'painters', 'Position', [10 10 900 450])

% scatter FLP type
index = find(tableGPCR_TimesInLoop.FLPtype);
scatter(index,tableGPCR_TimesInLoop.TimesInLoop(index),SizeMarker,'MarkerEdgeColor',ColorFLP,...
              'MarkerFaceColor',ColorFLP,...
              'LineWidth',1.5), hold on

% scatter NLP type
index = find(tableGPCR_TimesInLoop.NLPtype);
scatter(index,tableGPCR_TimesInLoop.TimesInLoop(index),SizeMarker,'MarkerEdgeColor',ColorNLP,...
              'MarkerFaceColor',ColorNLP,...
              'LineWidth',1.5), 
          
          
% scatter promiscous GPCR
index = find(tableGPCR_TimesInLoop.promiscuousGPCR);
scatter(index, tableGPCR_TimesInLoop.TimesInLoop(index), SizeMarker,'MarkerEdgeColor',[0 0 0],...
              'LineWidth',1), 


xlim([1 height(tableGPCR_TimesInLoop)])
ylabel('Number of times in autocrine connection')
xticks('')
ylim([0 180])
yticks(0:60:180)

set(gcf,'color','w');
%export_fig 'NumberLoopsForEachGPCR.pdf' -painters

%% Rank each individual NPP based on its frequency in autocrine connections

V = table2cell(tableTimesInLoop(:, 'GeneLigand'));
idx = unique(V); % get list of all unique NPP

tableNPP_TimesInLoop = table; % initiate empty table
tableNPP_TimesInLoop.GeneLigand = idx;

for ii = 1:length(idx) %iterate over each NPP and sum the number of times it figures in a loop
    NameNPP = idx(ii);
    index = find(string(tableTimesInLoop.GeneLigand) == NameNPP);
    tableNPP_TimesInLoop.TimesInLoop(ii) = sum(tableTimesInLoop.TimesInLoop(index));
    tableNPP_TimesInLoop.promiscuousGPCR(ii) = tableTimesInLoop.promiscuousGPCR(index(1));
    tableNPP_TimesInLoop.FLPtype(ii) = tableTimesInLoop.FLPtype(index(1));
    tableNPP_TimesInLoop.NLPtype(ii) = tableTimesInLoop.NLPtype(index(1));
end

tableNPP_TimesInLoop = sortrows(tableNPP_TimesInLoop, 2); % order the table in ascending order (how many autocrine loops)
figure('Renderer', 'painters', 'Position', [10 10 900 450])

% scatter FLP type
index = find(tableNPP_TimesInLoop.FLPtype);
scatter(index,tableNPP_TimesInLoop.TimesInLoop(index),SizeMarker,'MarkerEdgeColor',ColorFLP,...
              'MarkerFaceColor',ColorFLP,...
              'LineWidth',1.5), hold on

% scatter NLP type
index = find(tableNPP_TimesInLoop.NLPtype);
scatter(index,tableNPP_TimesInLoop.TimesInLoop(index),SizeMarker,'MarkerEdgeColor',ColorNLP,...
              'MarkerFaceColor',ColorNLP,...
              'LineWidth',1.5), 
          
          
% scatter promiscous GPCR
index = find(tableNPP_TimesInLoop.promiscuousGPCR);
scatter(index, tableNPP_TimesInLoop.TimesInLoop(index), SizeMarker,'MarkerEdgeColor',[0 0 0],...
              'LineWidth',1), 


xlim([1 height(tableNPP_TimesInLoop)])
ylabel('Number of times in autocrine connection')
xticks('')
ylim([0 120])
yticks(0:40:120)

set(gcf,'color','w');
%export_fig 'NumberLoopsForEachNPP.pdf' -painters

%% Plot number of autocrine connections per type (FLP or NLP)
figure('Renderer', 'painters', 'Position', [10 10 300 450])

% subset based type of connection
index = find(tableTimesInLoop.FLPtype);
NumLoopsFLP = tableTimesInLoop.TimesInLoop(index)';
FLPPromiscuous = tableTimesInLoop.promiscuousGPCR(index)';

index = find(tableTimesInLoop.NLPtype);
NumLoopsNLP = tableTimesInLoop.TimesInLoop(index)';
NLPPromiscuous = tableTimesInLoop.promiscuousGPCR(index)';

% format data for plotting
numrows = max([length(NumLoopsFLP), length(NumLoopsNLP)]);
DataMatrix = nan(numrows, 2);
PromiscuousMatrix = DataMatrix;
% insert data in DataMatrix
DataMatrix(1:length(NumLoopsFLP), 1) = NumLoopsFLP';
DataMatrix(1:length(NumLoopsNLP), 2) = NumLoopsNLP';

% insert data in PromiscuousMatrix
PromiscuousMatrix(1:length(NumLoopsFLP), 1) = FLPPromiscuous';
PromiscuousMatrix(1:length(NumLoopsNLP), 2) = NLPPromiscuous';

boxplot(DataMatrix, 'Notch', 'on'), hold on
    
% Delete outliers stars on boxplot (does not delete actual outliers)
h = findobj(gca,'tag','Outliers'); 
delete(h);

% Color of boxplots
h = findobj(gca,'Tag','Box');
Position = 1:length(h);
Position = fliplr(Position); %flip order to correct for weird problem
patch(get(h(Position(1)),'XData'),get(h(Position(1)),'YData'), ColorFLP,'FaceAlpha',.5);
patch(get(h(Position(2)),'XData'),get(h(Position(2)),'YData'), ColorNLP,'FaceAlpha',.5);

% Generate ScatterPlots (code adapted from NotBoxPlot())
jitter = 0.3;
X = 1:size(DataMatrix, 2);

% first, scatter FLP-type connections
thisY = DataMatrix(:,1);
thisY = thisY(~isnan(thisY));
thisX = repmat(X(1), 1, length(thisY));
J = (rand(size(thisX))-0.5) * jitter;
scatter(thisX + J, thisY, SizeMarker, 'Marker', 'o', 'MarkerEdgeColor', ColorFLP,'MarkerFaceColor',ColorFLP,'LineWidth',1.5);
indexPro = find(PromiscuousMatrix(:,1) == 1);
scatter(thisX(indexPro) + J(indexPro), thisY(indexPro), SizeMarker, 'Marker', 'o', 'MarkerEdgeColor', 'k','LineWidth',1.5);

% scatter NLP-type connections
thisY = DataMatrix(:,2);
thisY = thisY(~isnan(thisY));
thisX = repmat(X(2), 1, length(thisY));
J = (rand(size(thisX))-0.5) * jitter;
scatter(thisX + J, thisY, SizeMarker, 'Marker', 'o', 'MarkerEdgeColor', ColorNLP,'MarkerFaceColor',ColorNLP,'LineWidth',1.5);
indexPro = find(PromiscuousMatrix(:,2) == 1);
scatter(thisX(indexPro) + J(indexPro), thisY(indexPro), SizeMarker, 'Marker', 'o', 'MarkerEdgeColor', 'k','LineWidth',1.5);

set(gca,'XTickLabel',{'FLP', 'NLP'});
ylabel('Number of Self Loops')
ylim([0 60])
yticks(0:20:60)

set(gcf,'color','w');
%export_fig 'SupplFig_NumberLoopsForType.pdf' -painters
% todo: statistics in GraphPad

%% Plot number of autocrine connections depending on GPCR being promiscuous or not
% subset based type of connection
index = find(tableTimesInLoop.promiscuousGPCR == 1);
NumLoopsPromiscuous = tableTimesInLoop.TimesInLoop(index)';
FLP_Promiscuous = tableTimesInLoop.FLPtype(index)';
NLP_Promiscuous = tableTimesInLoop.NLPtype(index)';

index = find(tableTimesInLoop.promiscuousGPCR == 0);
NumLoopsNotPromiscuous = tableTimesInLoop.TimesInLoop(index)';
FLP_NotPromiscuous = tableTimesInLoop.FLPtype(index)';
NLP_NotPromiscuous = tableTimesInLoop.NLPtype(index)';

numrows = max([length(NumLoopsPromiscuous), length(NumLoopsNotPromiscuous)]);
DataMatrix = nan(numrows, 2);
% insert data in DataMatrix
DataMatrix(1:length(NumLoopsNotPromiscuous), 1) = NumLoopsNotPromiscuous';
DataMatrix(1:length(NumLoopsPromiscuous), 2) = NumLoopsPromiscuous';

figure('Renderer', 'painters', 'Position', [10 10 300 450])
boxplot(DataMatrix, 'Notch', 'on'), hold on
% Delete outliers stars on boxplot (does not delete actual outliers)
h = findobj(gca,'tag','Outliers'); 
delete(h);

% Color of boxplots
h = findobj(gca,'Tag','Box');
Position = 1:length(h);
Position = fliplr(Position); %flip order to correct for weird problem
patch(get(h(Position(1)),'XData'),get(h(Position(1)),'YData'), [0.7 0.7 0.7],'FaceAlpha',.5);
patch(get(h(Position(2)),'XData'),get(h(Position(2)),'YData'), [0.7 0.7 0.7],'FaceAlpha',.5);

% Generate ScatterPlots (code adapted from NotBoxPlot())
jitter = 0.3;
X = 1:size(DataMatrix, 2);

% scatter NumLoopsNotPromiscuous
thisY = DataMatrix(:,1);
thisY = thisY(~isnan(thisY));
thisX = repmat(X(1), 1, length(thisY));
J = (rand(size(thisX))-0.5) * jitter;

indexFLP = find(FLP_NotPromiscuous == 1);
scatter(thisX(indexFLP) + J(indexFLP), thisY(indexFLP), SizeMarker, 'Marker', 'o', 'MarkerFaceColor', ColorFLP,'MarkerEdgeColor', ColorFLP,'LineWidth',1);
indexNLP = find(NLP_NotPromiscuous == 1);
scatter(thisX(indexNLP) + J(indexNLP), thisY(indexNLP), SizeMarker, 'Marker', 'o', 'MarkerFaceColor', ColorNLP,'MarkerEdgeColor', ColorNLP,'LineWidth',1);

% scatter NumLoopsPromiscuous
thisY = DataMatrix(:,2);
thisY = thisY(~isnan(thisY));
thisX = repmat(X(2), 1, length(thisY));
J = (rand(size(thisX))-0.5) * jitter;

indexFLP = find(FLP_Promiscuous == 1);
scatter(thisX(indexFLP) + J(indexFLP), thisY(indexFLP), SizeMarker, 'Marker', 'o', 'MarkerFaceColor', ColorFLP,'MarkerEdgeColor', 'k','LineWidth',1);
indexNLP = find(NLP_Promiscuous == 1);
scatter(thisX(indexNLP) + J(indexNLP), thisY(indexNLP), SizeMarker, 'Marker', 'o', 'MarkerFaceColor', ColorNLP,'MarkerEdgeColor', 'k','LineWidth',1);

ylim([0 60])
yticks(0:20:60)
xticklabels('')
set(gcf,'color','w');
%export_fig 'SupplFig_Number_AutocrineConnections_promiscuous or not.pdf' -painters
%% Plot fraction of promiscuous receptors having autocrine connections
FractionAutocrineLoops_Notpromicuous = length(find(NumLoopsNotPromiscuous ~= 0))/length(NumLoopsNotPromiscuous);
FractionAutocrineLoops_promicuous = length(find(NumLoopsPromiscuous ~= 0))/length(NumLoopsPromiscuous);

figure('Renderer', 'painters', 'Position', [10 10 300 450])
b = bar([FractionAutocrineLoops_Notpromicuous, FractionAutocrineLoops_promicuous]);
b.FaceColor = 'flat';
% b.CData(1,:) = [0.8    0.8    0.8];
% b.CData(2,:) = [0.8    0.8    0.8];
yticks(0:0.25:1)
yticklabels('')
xticklabels('')
ylim([0 1])
xlim([0.5 2.5])

set(gcf,'color','w');
%export_fig 'SupplFig_Fraction_AutocrineConnections_promiscuous or not.pdf' -painters