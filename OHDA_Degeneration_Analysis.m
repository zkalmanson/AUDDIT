% ANALYZE 6-OHDA IMAGES
% Andrew Clark
% 05-12-2022
% Updated: 08-30-2022

clc;
clear;
close all
pxlSize = 0.38/4; % Divided by 4 because of the resize

% Load all the feature data from the 6-OHDA experiments
load('ohda-FeatureTable.mat');
load('ohda-CroppedImages.mat');
load('ohda-BreakTable.mat'); 
% Get width info
for i = 1:size(dendriteWidths,2)
    folderWidth = dendriteWidths(:,i);
    for ii = 1:size(dendriteWidths(:,i),1)
        tempWidth = cell2mat(folderWidth(ii));
        medWidth{ii,i} = median(tempWidth,1)';
        meanWidth{ii,i} = mean(tempWidth,1)';
    end
end

%% Dendrite only Data

% Initialization for figure production
cMap = {'k',[.5 0 0],[0 0 .75],'k',[.5 0 0],[0 0 .75],'k',[.5 0 0],[0 0 .75],'k',[.5 0 0],[0 0 .75]};
xLoc = [1:3,5:7,9:11,13:15];
figure();

% Organize data
tempMetric = dendriteLengths; %%% METRIC SELECTION HERE %%%%
groupData = [];
xGroup = [];
keyPlot = [];
xkey = [];
xExpGroup = [];

for i = 1:size(tempMetric,2)

    % Place into metrics
    dMetric = cell2mat(tempMetric(:,i)).*pxlSize; % place into matrix/vector
    if size(dMetric,2)>1
        dMetric = dMetric';
        dMetric = reshape(dMetric,[],1);% only necessary for breakage
    end
    xMetric = ones(length(dMetric),1)*xLoc(i);
    
    % Organize groups into experiments
    if ismember(i,1:3:10)
        exp = 1;
    elseif ismember(i,2:3:11)
        exp = 2;
    elseif ismember(i,3:3:12)
        exp = 3;
    end

    % Swarm chat with all of the data
    swarmchart(xMetric,dMetric,50,cMap{i},'filled','MarkerFaceAlpha',0.95)
    hold on 

    % Needed for experimental grouping
    xExp = ones(length(dMetric),1)*exp;

    % Place into ANOVA/boxchart
    groupData = [groupData;dMetric];
    xGroup = [xGroup;xMetric];
    xExpGroup = [xExpGroup;xExp];

end

% Overlay box chart
cMap = {'k',[.5 0 0],[0 0 .5],'c','k',[.5 0 0],[0 0 .5],'c','k',[.5 0 0],[0 0 .5],'c','k',[.5 0 0],[0 0 .5],'c'};
for i =1:max(xGroup)
        idx = find(xGroup==i);
        x = xGroup(idx);
        y = groupData(idx);
        boxchart(x,y,'MarkerStyle','none','BoxFaceColor',cMap{i},'LineWidth',4,'WhiskerLineColor',cMap{i})
end
legend off

% Adjust size
set(gcf,'Position',[2 71 1105 906]);

% Adjust visualizations
title('Dendrite Remaining','FontSize',40)
ax=gca;
ax.FontSize = 30;
xticks([2,6,10,14]);
xticklabels({'Control','10mM','25 mM','50 mM'})
xlabel('Experimental Group','FontSize',30)
ylabel('Dendrite Remaining (au)','FontSize',30)
% ylim([0 255])
        
% ANOVA
[~, ~, stats] = anova1(groupData,xGroup);
results = multcompare(stats,'CType','bonferroni');



%% CLUSTERING AND FILTERING

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
clusterTable.BlebLocation = (1-(MasterFeatureTable.Centroid(:,1)./MasterFeatureTable.DendriteLength));
clusterTable.BlebLocation(clusterTable.BlebLocation>1)=1; % Deal with short dendrites for now
MasterFeatureTable.BlebLocation = (1-(MasterFeatureTable.Centroid(:,1)./MasterFeatureTable.DendriteLength));
MasterFeatureTable.BlebLocation(MasterFeatureTable.BlebLocation>1)=1;

% Create cluster matrix
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

%% Feature results (area, min width, feature intensity, etc)

% Define experiments
groupLabel = {'Control';'10 mM';'25 mM';'50 mM'};
cMap = {'k',[.5 0 0],[0 0 .75],'k',[.5 0 0],[0 0 .75],'k',[.5 0 0],[0 0 .75],'k',[.5 0 0],[0 0 .75]};

% Initialize
clustTable = filtMasterTable;
k = 0;
groupData = [];
xGroup = [];
xLoc = [1:3,5:7,9:11,13:15];
figure();

% For loop of experiments 
for i = 1:3 % THESE ARE
    k=0;  
    for j = 1:4
        xMet = [];
        dMet = [];
        k=k+1;
        tempIdx = [i,i+3,i+6,i+9];
        expIdx= find(clustTable.Folder==tempIdx(j));
        expData = clustTable(expIdx,:);
        expData = expData.MeanIntensity;  %%% METRIC SELECTION %%%

        % Swarm plot gather
        dMet = [dMet;expData];
        xMet = [xMet;ones(length(expData),1)*xLoc(tempIdx(j))];

        % Swarm chat with all of the data
        swarmchart(xMet,dMet,50,cMap{i},'filled','MarkerFaceAlpha',0.95)
        hold on 

        % Gather for ANOVA
        groupData = [groupData;dMet];
        xGroup = [xGroup;xMet];

    end
end

% Overlay box chart
cMap = {'k',[.5 0 0],[0 0 .5],'c','k',[.5 0 0],[0 0 .5],'c','k',[.5 0 0],[0 0 .5],'c','k',[.5 0 0],[0 0 .5],'c'};
for i =1:max(xGroup)
        idx = find(xGroup==i);
        x = xGroup(idx);
        y = groupData(idx);
        boxchart(x,y,'MarkerStyle','none','BoxFaceColor',cMap{i},'LineWidth',4,'WhiskerLineColor',cMap{i})
end

% Make pretty
title('Feature Intensity','fontsize',40)
ax = gca;
ax.FontSize = 30;
xticks(2:4:14)
xticklabels(groupLabel)
ylabel('Intensity (au)','FontSize',30)
xlabel('Experimental Group','FontSize',30)
ylim([0 300])

% Adjust size
set(gcf,'Position',[2 71 1105 906]);

% Now do one way ANOVA
[~, ~, stats] = anova1(groupData,xGroup);
results = multcompare(stats,'CType','bonferroni');

%% Feature per dendrite results (bleb count, etc)

% Blebs per dendrite length
close all

% Define experiments
groupLabel = {'Control';'10 mM';'25 mM';'50 mM'};
cMap = {'k',[.5 0 0],[0 0 .75]};

% Initialize
ftPerDend = [];
ftPerDendLen = [];
ftSpread = [];
tempMetric = [];
tempXmetric = [];
ftCount = [];
ftCountLen = [];
groupDataPlot = [];
xGroupPlot = [];

clustTable = filtMasterTable;
zz = 0; % counter for storing feature counts
dendIdx = 0; % counter for comparing dendrite statistics
xCount = [1:3,5:7,9:11,13:15];

groupDataMetric = [];
xGroupMetric = [];

% For loop of experiments 
for i = 1:3 % WHICH EXPERIMENTS TO RUN

    groupData = [];
    xGroup = [];
    groupDataLen = [];
    xGroupLen = [];
    groupDataSpread = [];
    xGroupSpread = [];
    groupDataW = [];
    xGroupW = [];
    k=0;
    for j = 1:4 % WHICH GROUPS TO RUN
        k=k+1;
        zz=zz+1;

        % Initialize
        ftPerDend = [];
        ftPerDendLen = [];
        ftSpread = [];
        tempMetric = [];
        tempXmetric = [];
        dMet = [];
        xMet = [];
        tempIdx = [i,i+3,i+6,i+9];
        
        % Find how many images to loop over
        maxImgIdx= find(clustTable.Folder==tempIdx(j));
        maxImgTable = clustTable(maxImgIdx,:);
        maxImg = max(maxImgTable.ImageNumber);

        % Third for loop for each image in each folder
        for img = 1:maxImg
            for dend = 1:4
                dendIdx = dendIdx+1;
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
                    dendBreak = 1;
                end

                ftPerDend(img,dend) = expCount;
                ftPerDendLen(img,dend) = expCount./totDendLength;

                % Weight feature count by amount of dendrite remaining
                weightFeature = expCount./dendBreak;

                % Gather for ANOVA - features per dendrite
                groupData = [groupData;expCount];
                xGroup = [xGroup;ones(length(expCount),1)*k]; 

                % Gather for ANOVA - features per dendrite length
                groupDataLen = [groupDataLen;ftPerDendLen(img,dend)];
                xGroupLen = [xGroupLen;ones(length(expCount),1)*k];

                % Gather for ANOVA - features per dendrite length
                groupDataW = [groupDataW;weightFeature];
                xGroupW = [xGroupW;ones(length(weightFeature),1)*k]; 

                % Gather for ANOVA - features spread
                groupDataSpread = [groupDataSpread;ftSpread];
                xGroupSpread = [xGroupSpread;ones(length(ftSpread),1)*k];

                % Gather for ANOVA - dendrite metric
                groupDataMetric = [groupDataMetric;weightFeature];
                xGroupMetric = [xGroupMetric;ones(length(weightFeature),1)*xCount(tempIdx(j))]; 
                
                %%%% METRIC SELECTION FOR PLOTTING %%%
                % Other options shown above %
                dMet = [dMet;ftPerDendLen(img,dend)];
                xMet = [xMet;ones(length(ftPerDendLen(img,dend)),1).*xCount(tempIdx(j))];
                 
            end
        end
        % Now plot
        swarmchart(xMet,dMet,50,cMap{i},'filled','MarkerFaceAlpha',0.95)
        hold on
    
        % Place into ANOVA/boxchart
        groupDataPlot = [groupDataPlot;dMet];
        xGroupPlot = [xGroupPlot;xMet];
    
        % Store data for binning and comparisons
        ftCount{zz} = ftPerDend;
        ftCountLen{zz} = ftPerDendLen;
    end
end


% Overlay box chart
cMap = {'k',[.5 0 0],[0 0 .5],'c','k',[.5 0 0],[0 0 .5],'c','k',[.5 0 0],[0 0 .5],'c','k',[.5 0 0],[0 0 .5],'c'};
for i =1:max(xGroupPlot)
        idx = find(xGroupPlot==i);
        x = xGroupPlot(idx);
        y = groupDataPlot(idx);
        boxchart(x,y,'MarkerStyle','none','BoxFaceColor',cMap{i},'LineWidth',4,'WhiskerLineColor',cMap{i})
end


% Make pretty    
title('Weighted Feature Count','fontsize',40)
ax = gca;
ax.FontSize = 30;
xticks([2:4:14])
xticklabels(groupLabel)
ylabel('Weighted Feature Count','FontSize',30)
xlabel('Experiment Group','FontSize',30)
ylim([0 20])

% Adjust size
set(gcf,'Position',[2 71 1105 906]);

% Now do one way ANOVA over for loop
[p, tbl, stats] = anova1(groupDataPlot,xGroupPlot);
results = multcompare(stats,'CType','bonferroni');
