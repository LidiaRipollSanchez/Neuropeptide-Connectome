%% Plot co-expression in O2 sensory circuit and VNC
% Jan Watteyne
% 20230608

clear all
close all
workspace

%% Import data
%load in distance file - from Petra
load('Worm279.mat');
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
% CeNGEN threshold 4 - EC50 500 nM
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

%% PLOT NUMBER OF AC ON XY COORDINATES
% % % 2D position
% figure,
plot(Worm279_positions(:,1), Worm279_positions(:,2),'k.')
axis equal
axis off
set(gcf,'color','w');
% %export_fig 'CorrelationNumberLoopsVsNumberNPPs.pdf' -painters

%% 2D positions of neurons with number of self loops as size of dots
%------------------------------------------------
% Retrieve neuron coordinates for all the neurons
%------------------------------------------------
NumberAC_coordinate = [cell2mat(NeuronAC(:,2)), nan(size(NeuronAC, 1),2)];
% col 1: number of AC / col 2: x coordinate / col 3: y coordinate

for ii = 1:size(NumberAC_coordinate, 1)
    try
        Index = find(contains(Worm279_labels, NeuronAC{ii,1}));
        NumberAC_coordinate(ii, 2) = Worm279_positions(Index, 1);
        NumberAC_coordinate(ii, 3) = Worm279_positions(Index, 2);
    catch 
    end
end
%%
%-------------------------------------------------------------
% Scatter number of self loops on top of total nervous system
%-------------------------------------------------------------
coordinatesNoSelfLoops = NumberAC_coordinate(isnan(NumberAC_coordinate(:,1)), 2:3);
coordinatesSelfLoops = NumberAC_coordinate(~isnan(NumberAC_coordinate(:,1)), 1:3);

col = [0.2 0.2 0.2];
figure('Renderer', 'painters', 'Position', [10 10 1300 1300])
scatter(coordinatesNoSelfLoops(:, 1), coordinatesNoSelfLoops(:, 2), '.', 'MarkerEdgeColor', col) % first, scatter coordinates neurons without SL
hold on
scatter(coordinatesSelfLoops(:, 2), coordinatesSelfLoops(:, 3),coordinatesSelfLoops(:,1)*10, 'o', 'MarkerEdgeColor', col)
axis equal
axis off
set(gcf,'color','w');

% Additional: color certain neurons with a color
% find the index of all entries in which a certain string is found

% colors chosen using colorBrewer
colA = [0.5529    0.8275    0.7804];
colB = [0.9843    0.5020    0.4471];
colAS = [1.0000    1.0000    0.7020];
colD = [0.5020    0.6941    0.8275];
colOsens = [0.7843    0.5451    0.7490];

Index = [find(count(string(NeuronAC(:,1)), 'URX')); find(count(string(NeuronAC(:,1)), 'PQR')); find(count(string(NeuronAC(:,1)), 'PVR'))];
scatter(NumberAC_coordinate(Index, 2), NumberAC_coordinate(Index, 3), NumberAC_coordinate(Index,1)*10, 'o','MarkerEdgeColor', colOsens)

Index = find(count(string(NeuronAC(:,1)), 'DA0'));
scatter(NumberAC_coordinate(Index, 2), NumberAC_coordinate(Index, 3), NumberAC_coordinate(Index,1)*10, 'o','MarkerEdgeColor', colA)

Index = find(count(string(NeuronAC(:,1)), 'DB0'));
scatter(NumberAC_coordinate(Index, 2), NumberAC_coordinate(Index, 3), NumberAC_coordinate(Index,1)*10, 'o','MarkerEdgeColor', colB)

Index = [find(count(string(NeuronAC(:,1)), 'VA0')); find(count(string(NeuronAC(:,1)), 'VA1'))];
scatter(NumberAC_coordinate(Index, 2), NumberAC_coordinate(Index, 3), NumberAC_coordinate(Index,1)*10, 'o','MarkerEdgeColor', colA)

Index = [find(count(string(NeuronAC(:,1)), 'VB0')); find(count(string(NeuronAC(:,1)), 'VB1'))];
scatter(NumberAC_coordinate(Index, 2), NumberAC_coordinate(Index, 3), NumberAC_coordinate(Index,1)*10, 'o','MarkerEdgeColor', colB)

Index = [find(count(string(NeuronAC(:,1)), 'DD0')); find(count(string(NeuronAC(:,1)), 'VD0')); find(count(string(NeuronAC(:,1)), 'VD1'))];
scatter(NumberAC_coordinate(Index, 2), NumberAC_coordinate(Index, 3), '.', 'MarkerEdgeColor',colD, 'MarkerFaceColor',colD);

Index = [find(count(string(NeuronAC(:,1)), 'AS0')); find(count(string(NeuronAC(:,1)), 'AS10')); find(count(string(NeuronAC(:,1)), 'AS11'))];
scatter(NumberAC_coordinate(Index, 2), NumberAC_coordinate(Index, 3), '.', 'MarkerEdgeColor',colAS, 'MarkerFaceColor',colAS);

%export_fig 'MainFig_2DcoordinatesNumberAC.pdf' -painters

%% plot number of AC for each motoneuron class
IndexDA = find(count(string(NeuronAC(:,1)), 'DA0'));
IndexDB = find(count(string(NeuronAC(:,1)), 'DB0'));
IndexVA = [find(count(string(NeuronAC(:,1)), 'VA0')); find(count(string(NeuronAC(:,1)), 'VA1'))];
IndexVB = [find(count(string(NeuronAC(:,1)), 'VB0')); find(count(string(NeuronAC(:,1)), 'VB1'))];
IndexDD = find(count(string(NeuronAC(:,1)), 'DD0'));
IndexVD = [find(count(string(NeuronAC(:,1)), 'VD0')); find(count(string(NeuronAC(:,1)), 'VD1'))];
IndexAS = [find(count(string(NeuronAC(:,1)), 'AS0')); find(count(string(NeuronAC(:,1)), 'AS10')); find(count(string(NeuronAC(:,1)), 'AS11'))];

DataMatrix = nan(max([length(IndexDA), length(IndexDB),length(IndexVA),length(IndexVB),length(IndexDD),length(IndexVD),length(IndexDA),length(IndexAS)]), 7);

DataMatrix(1:length(IndexDA), 1) = cell2mat(NeuronAC(IndexDA,2));
DataMatrix(1:length(IndexVA), 2) = cell2mat(NeuronAC(IndexVA,2));

a = cell2mat(NeuronAC(IndexDD,2));
a(isnan(a)) = 0;
DataMatrix(1:length(IndexDD), 3) = a;

a = cell2mat(NeuronAC(IndexVD,2));
a(isnan(a)) = 0;
DataMatrix(1:length(IndexVD), 4) = a;

DataMatrix(1:length(IndexDB), 5) = cell2mat(NeuronAC(IndexDB,2));
DataMatrix(1:length(IndexVB), 6) = cell2mat(NeuronAC(IndexVB,2));

a = cell2mat(NeuronAC(IndexAS,2));
a(isnan(a)) = 0;
DataMatrix(1:length(IndexAS), 7) = a;


figure('Renderer', 'painters', 'Position', [10 10 800 450]) %initiate figure with dimensions specified
boxplot(DataMatrix, 'Notch', 'on'), hold on %plot boxplot, with 'Notch' setting, which I like most
SizeMarker = 20;

% Delete outliers 'red' stars on boxplot (if there are any). This does not
% delete actual outliers, only their scatterplot
h = findobj(gca,'tag','Outliers'); 
delete(h);

% Specify the color of each individual boxplot
% h = findobj(gca,'Tag','Box');
% Position = 1:length(h);
% Position = fliplr(Position); %flip order to correct for order switching
% patch(get(h(Position(1)),'XData'),get(h(Position(1)),'YData'), Color,'FaceAlpha',.5);
% patch(get(h(Position(2)),'XData'),get(h(Position(2)),'YData'), Color,'FaceAlpha',.5);
% patch(get(h(Position(3)),'XData'),get(h(Position(3)),'YData'), Color,'FaceAlpha',.5); 
% patch(get(h(Position(4)),'XData'),get(h(Position(4)),'YData'), Color,'FaceAlpha',.5); 
% patch(get(h(Position(5)),'XData'),get(h(Position(5)),'YData'), Color,'FaceAlpha',.5); 
% patch(get(h(Position(6)),'XData'),get(h(Position(6)),'YData'), Color,'FaceAlpha',.5); 
% patch(get(h(Position(7)),'XData'),get(h(Position(7)),'YData'), [1 0.5 0],'FaceAlpha',.5); 

% Generate ScatterPlots (code adapted from NotBoxPlot() function of FEX)
jitter = 0.5; %specify 'jitter'
X = 1:size(DataMatrix, 2);

% FIRST, A-CLASS MOTONEURONS 
thisY = DataMatrix(:,1); % colum 1 %remark: can be automated but this gives you easy control of each boxplot
thisY = thisY(~isnan(thisY));
thisX = repmat(X(1), 1, length(thisY));
J = (rand(size(thisX))-0.5) * jitter;
scatter(thisX + J, thisY, SizeMarker, 'Marker', 'o', 'MarkerEdgeColor', colA,'MarkerFaceColor',colA,'LineWidth',1.5);

thisY = DataMatrix(:,2); % colum 2
thisY = thisY(~isnan(thisY));
thisX = repmat(X(2), 1, length(thisY));
J = (rand(size(thisX))-0.5) * jitter;
scatter(thisX + J, thisY, SizeMarker, 'Marker', 'o', 'MarkerEdgeColor', colA,'MarkerFaceColor',colA,'LineWidth',1.5);

% D -CLASS MOTONEURONS - RED
thisY = DataMatrix(:,3); % colum 3
thisY = thisY(~isnan(thisY));
thisX = repmat(X(3), 1, length(thisY));
J = (rand(size(thisX))-0.5) * jitter;
scatter(thisX + J, thisY, SizeMarker, 'Marker', 'o', 'MarkerEdgeColor', colD,'MarkerFaceColor',colD,'LineWidth',1.5);

thisY = DataMatrix(:,4); % colum 4
thisY = thisY(~isnan(thisY));
thisX = repmat(X(4), 1, length(thisY));
J = (rand(size(thisX))-0.5) * jitter;
scatter(thisX + J, thisY, SizeMarker, 'Marker', 'o', 'MarkerEdgeColor', colD,'MarkerFaceColor',colD,'LineWidth',1.5);

% B -CLASS MOTONEURONS - BLUE
thisY = DataMatrix(:,5); 
thisY = thisY(~isnan(thisY));
thisX = repmat(X(5), 1, length(thisY));
J = (rand(size(thisX))-0.5) * jitter;
scatter(thisX + J, thisY, SizeMarker, 'Marker', 'o', 'MarkerEdgeColor', colB,'MarkerFaceColor',colB,'LineWidth',1.5);

thisY = DataMatrix(:,6);
thisY = thisY(~isnan(thisY));
thisX = repmat(X(6), 1, length(thisY));
J = (rand(size(thisX))-0.5) * jitter;
scatter(thisX + J, thisY, SizeMarker, 'Marker', 'o', 'MarkerEdgeColor', colB,'MarkerFaceColor',colB,'LineWidth',1.5);

% AS -CLASS MOTONEURONS - ORANGE
thisY = DataMatrix(:,7); 
thisY = thisY(~isnan(thisY));
thisX = repmat(X(7), 1, length(thisY));
J = (rand(size(thisX))-0.5) * jitter;
scatter(thisX + J, thisY, SizeMarker, 'Marker', 'o', 'MarkerEdgeColor', colAS,'MarkerFaceColor',colAS,'LineWidth',1.5);


%set(gca,'XTickLabel',{'colum 1', 'colum 2', 'colum 3'});
ylabel('label y axis')
ylim([0 8])
yticks(0:2:8)

set(gcf,'color','w'); %white background of figure for export
%export_fig 'NumberAC_VNC.pdf' -painters %use export_fig function (on FEX) to export figure as PDF, which I then open in Inkscape to manually tweak
