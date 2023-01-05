% COMPARE 6-OHDA IMAGES TO MANUAL SCORING DONE BY KATE MORTON
% Andrew Clark
% 06-22-2022
% Updated: 08-30-2022

clc;
clear;
close all
pxlSize = 0.38/4; % Divided by 4 because of the resize

% Load all the feature data from the rotenone experiments
load('ohda-FeatureTable.mat');
load('ohda-CroppedImages.mat');
load('ohda-BreakTable.mat');
load('ManualScoring.mat');


%% Dendrite only Data
close all

% Organize data
tempMetric = dendriteBreakProp;
groupData = [];
xGroup = [];
keyPlot = [];
xkey = [];
xExpGroup = [];
groupDataCat = [];
xGroupCat = [];


for i = 1:size(tempMetric,2)

    % Organize groups into experimental groups
    if ismember(i,1:3)
        k = 1;
    elseif ismember(i,4:6)
        k = 2;
    elseif ismember(i,7:9)
        k = 3;
    elseif ismember(i,10:12)
        k = 4;
    end

    % Place into metrics
    dMetric = 100.*(1-cell2mat(tempMetric(:,i))); % place into matrix/vector
    dMetric = reshape(dMetric,[],1);
    xMetric = ones(length(dMetric),1)*i;
    
    % Organize groups into experiments
    if ismember(i,1:3:10)
        exp = 1;
    elseif ismember(i,2:3:11)
        exp = 2;
    elseif ismember(i,3:3:12)
        exp = 3;
    end

    % Needed for experimental grouping
    xExp = ones(length(dMetric),1)*exp;

    % Place into ANOVA/boxchart
    groupData= [groupData;dMetric];
    xGroup = [xGroup;xMetric];

    groupDataCat= [groupDataCat;dMetric];
    xGroupCat = [xGroupCat;xMetric];

    xExpGroup = [xExpGroup;xExp];

end


%% Adjust feature table

% Create a new table to put into k-means/k-medoids
clusterTable = [];
clusterTable.Area = MasterFeatureTable.Area;
clusterTable.MajorAxisLength = MasterFeatureTable.MajorAxisLength;
clusterTable.MinorAxisLength = MasterFeatureTable.MinorAxisLength;
clusterTable.Circularity = MasterFeatureTable.Circularity;
clusterTable.Eccentricity = MasterFeatureTable.Eccentricity;
clusterTable.Extent = MasterFeatureTable.Extent;
clusterTable.Orientation = MasterFeatureTable.Orientation;
clusterTable.ConvexArea = MasterFeatureTable.ConvexArea;
clusterTable.Solidity = MasterFeatureTable.Solidity;
clusterTable.MeanIntensity = MasterFeatureTable.MeanIntensity;
clusterTable.MaxFeretDiameter = MasterFeatureTable.MaxFeretDiameter;
clusterTable.MinFeretDiameter = MasterFeatureTable.MinFeretDiameter;

% Location
clusterTable.BlebLocation = MasterFeatureTable.Centroid(:,1)./MasterFeatureTable.DendriteLength;
clusterTable.BlebLocation(clusterTable.BlebLocation>1)=1; 
MasterFeatureTable.BlebLocation = MasterFeatureTable.Centroid(:,1)./MasterFeatureTable.DendriteLength;
MasterFeatureTable.BlebLocation(MasterFeatureTable.BlebLocation>1)=1;

% RAW
clusterTable = struct2table(clusterTable);
cluster = table2cell(clusterTable);
clusterMatrix = cellfun(@double, cluster, 'uni', false);
clusterMatrix = cell2mat(clusterMatrix);

% Try filtering, if wanted
filterCluster = clusterTable;
filterCluster(filterCluster.MinorAxisLength < 1,:) = [];
filterCluster(filterCluster.Area < 1,:) = [];
filtMasterTable = MasterFeatureTable;
filtMasterTable(filtMasterTable.MinorAxisLength < 1,:) = [];
filtMasterTable(filtMasterTable.Area < 1,:) = [];
filterClusterCell = table2cell(filterCluster);
filterClusterMatrix = cellfun(@double, filterClusterCell, 'uni', false);
filterClusterMatrix = cell2mat(filterClusterMatrix);


%% Get counts per dendrite to "categorize" dendrites

% Blebs per dendrite length
close all

% Initialize
ftPerDend = [];
ftPerDendLen = [];
tempMetric = [];
tempXmetric = [];
countData = [];
xCount = [];
clustTable = filtMasterTable;
zz = 0; % counter for storing feature counts

% For loop of experiments 
for i = 1:3 % WHICH EXPERIMENTS TO RUN
    k=0;
    for j = 1:4 % WHICH GROUPS TO RUN
        k=k+1;

        % Initialize
        ftPerDend = [];
        ftPerDendLen = [];
        ftSpread = [];
        tempMetric = [];
        tempXmetric = [];
        tempIdx = [i,i+3,i+6,i+9];

        maxImgIdx= find(clustTable.Folder==tempIdx(j));
        maxImgTable = clustTable(maxImgIdx,:);
        maxImg = max(maxImgTable.ImageNumber);
        % Third for loop for image
        for img = 1:maxImg
            for dend = 1:4
                ftSpread = [];
                tempClustIdx = find(clustTable.Dendrite == dend & filtMasterTable.Folder == tempIdx(j) & filtMasterTable.ImageNumber == img);
                expData = clustTable(tempClustIdx,:);
                expCount = size(expData,1);                
                % Feature or bleb dependent
                if expCount > 0
                    dendLength = expData.DendriteLength(1)*pxlSize;
                    dendBreak = (1-expData.BreakLength(1));
                    totDendLength = dendLength.*dendBreak; % blebs per dendrite length actually there
                    % Get spread of features
                    if expCount > 3
                        ftSpread = std(expData.BlebLocation)*pxlSize;
                    end
                else
                    dendLength = 1*pxlSize;
                    totDendLength = 1*pxlSize;
                end

                ftPerDend(img,dend) = expCount;
                ftPerDendLen(img,dend) = expCount./totDendLength;
                weightFeature = expCount./dendBreak;
                
                % feature counts
                countData = [countData; expCount];
                xCount = [xCount;ones(length(expCount),1)*tempIdx(j)]; 
            end
        end
    end
end

%% This section puts feature counts and break percentages together
% Initialize and preallocated
scoreMatrix = zeros(12,5);

maxBreak = 80; % Dendrites over this percent are considered 'full'
cutBreak = 60; % Dendrites under this are consered really broken (score 4)
minFt = 0; % Less than this ft count is 'good'
cutFt = 5; % Cut off ft count to differentiate score 1 and score 2

% Score all dendrites
score0 = find(groupDataCat >= maxBreak & countData <= minFt);
score1 = find(groupDataCat >= maxBreak & countData > minFt & countData <= cutFt);
score2 = find(groupDataCat >= maxBreak & countData > cutFt);
score3 = find(groupDataCat < maxBreak & groupDataCat >= cutBreak);
score4 = find(groupDataCat < cutBreak);

% Place into experimental groups
for i = 1:size(score0,1)
    gIdx = xGroupCat(score0(i));
    scoreMatrix(gIdx,1) = scoreMatrix(gIdx,1)+1;
end
for i = 1:size(score1,1)
    gIdx = xGroupCat(score1(i));
    scoreMatrix(gIdx,2) = scoreMatrix(gIdx,2)+1;
end
for i = 1:size(score2,1)
    gIdx = xGroupCat(score2(i));
    scoreMatrix(gIdx,3) = scoreMatrix(gIdx,3)+1;
end
for i = 1:size(score3,1)
    gIdx = xGroupCat(score3(i));
    scoreMatrix(gIdx,4) = scoreMatrix(gIdx,4)+1;
end
for i = 1:size(score4,1)
    gIdx = xGroupCat(score4(i));
    scoreMatrix(gIdx,5) = scoreMatrix(gIdx,5)+1;
end
% Test flip
scoreMatFlip = scoreMatrix([1:3:10,2:3:11,3:3:12],:);

% Make scores into percentages
expN = sum(scoreMatrix,2);
scoreRate = (scoreMatrix./expN).*100;

% Rearrange matrix to fit with experimental groups
scoreFlip =  scoreRate([1:3:10,2:3:11,3:3:12],:);
exp1Score = scoreFlip(1:4,:);
exp2Score = scoreFlip(5:8,:);
exp3Score = scoreFlip(9:12,:);

% Put manual scoring all together together
manAll = [manual1;manual2;manual3];

%% PLOT DATA AGAINST MANUAL SCORING %%
close all

% Organize data for pltBarStackGroups
stackData = zeros(4,2,5);
groupLabels = {'Control','10 mM','25 mM','50 mM'};

stackData(:,1,:) = manual1;
stackData(:,2,:) = exp1Score;

% Plot with plotBarStackGroups
h = plotBarStackGroups(stackData,groupLabels);

% Adjust size
set(gcf,'Position',[2 71 1105 906]);
% Adjust visualizations
title('Manual to Code Comparison BY200','FontSize',40)
ax=gca;
ax.FontSize = 30;
xlabel('Experimental Group','FontSize',30)
ylabel('Scoring Rate','FontSize',30)
ylim([0 100])

% Adjust colors
set(h(1,:), 'FaceColor', 'Flat')
colors =  mat2cell(turbo(numel(h(1,:))),ones(numel(h(1,:)),1), 3);
set(h(1,:), {'CData'}, colors)
set(h(2,:), 'FaceColor', 'Flat')
set(h(1,:), 'EdgeColor', 'k');
set(h(2,:), 'EdgeColor', 'k');
colors =  mat2cell(turbo(numel(h(2,:))),ones(numel(h(2,:)),1), 3);
set(h(2,:), {'CData'}, colors)


hatchfill2(h(2,:),'single','HatchAngle',45,'HatchColor','k','HatchLineWidth',3);

%% Exp 2
% Organize data for pltBarStackGroups
stackData = zeros(4,2,5);
groupLabels = {'Control','10 mM','25 mM','50 mM'};

stackData(:,1,:) = manual2;
stackData(:,2,:) = exp2Score;

% Plot with plotBarStackGroups
h = plotBarStackGroups(stackData,groupLabels);

% Adjust size
set(gcf,'Position',[2 71 1105 906]);
% Adjust visualizations
title('Tub-1 Mutants','FontSize',40)
ax=gca;
ax.FontSize = 30;
xlabel('Experimental Group','FontSize',30)
ylabel('Scoring Rate','FontSize',30)
ylim([0 100])

% Adjust colors
set(h(1,:), 'FaceColor', 'Flat')
colors =  mat2cell(turbo(numel(h(1,:))),ones(numel(h(1,:)),1), 3);
set(h(1,:), {'CData'}, colors)
set(h(2,:), 'FaceColor', 'Flat')
colors =  mat2cell(turbo(numel(h(2,:))),ones(numel(h(2,:)),1), 3);
set(h(2,:), {'CData'}, colors)
hatchfill2(h(2,:),'single','HatchAngle',45,'HatchColor','k','HatchLineWidth',3);

%% Exp 3
% Organize data for pltBarStackGroups
stackData = zeros(4,2,5);
groupLabels = {'Control','10 mM','25 mM','50 mM'};

stackData(:,1,:) = manual3;
stackData(:,2,:) = exp3Score;

% Plot with plotBarStackGroups
h = plotBarStackGroups(stackData,groupLabels);

% Adjust size
set(gcf,'Position',[2 71 1105 906]);
% Adjust visualizations
title('Tub-2 Mutants','FontSize',40)
ax=gca;
ax.FontSize = 30;
xlabel('Experimental Group','FontSize',30)
ylabel('Scoring Rate','FontSize',30)
ylim([0 100])

% Adjust colors
set(h(1,:), 'FaceColor', 'Flat')
colors =  mat2cell(turbo(numel(h(1,:))),ones(numel(h(1,:)),1), 3);
set(h(1,:), {'CData'}, colors)
set(h(2,:), 'FaceColor', 'Flat')
colors =  mat2cell(turbo(numel(h(2,:))),ones(numel(h(2,:)),1), 3);
set(h(2,:), {'CData'}, colors)
hatchfill2(h(2,:),'single','HatchAngle',45,'HatchColor','k','HatchLineWidth',3);







%% LOAD DATA FROM DUKE RESULTS - Not needed since loaded%%
% % No Tub
% manual1 = [28.57142857	7.142857143	5.882352941	0;
% 62.5	55.71428571	51.47058824	26.08695652;
% 5.357142857	2.857142857	7.352941176	8.695652174;
% 3.571428571	21.42857143	22.05882353	21.73913043;
% 0	12.85714286	13.23529412	43.47826087];
% manual1=manual1';
% 
% % Tub 1
% manual2=[0	13.15789474	14.92537313	14.1025641;
% 84.375	72.36842105	62.68656716	60.25641026;
% 7.8125	2.631578947	2.985074627	5.128205128;
% 7.8125	11.84210526	8.955223881	7.692307692;
% 0	0	10.44776119	12.82051282];
% manual2=manual2';
% 
% % Tub 2
% manual3=[13.95348837	20	25.64102564	12.08791209;
% 72.09302326	77.14285714	65.38461538	42.85714286;
% 5.813953488	1.428571429	0	10.98901099;
% 8.139534884	1.428571429	3.846153846	16.48351648;
% 0	0	5.128205128	17.58241758];
% manual3=manual3';




 
