% Code used to analyze and obtain metrics from feature tables produced by
% the neurodegeneration master code.
% Date: 06-13-2022
% Updated: 08-30-2022
%% LOAD AND PREPARE DATA
clc;
clear;
close all

% Find the pixel size to get actual measurements
pxlSize = 0.1032; % Get from camera pixel size
pxlSize = pxlSize./4;% Divided by 4 because of the resize earlier
minWidth = 0;
minArea = 0; 

% Load data
load('cold-FeatureTable.mat');
load('cold-CroppedImages.mat');
load('cold-BreakTable.mat');

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
close all

% Initialization for figure production
cMap = {'k',[.5 0.5 0.5],[0 0 .75]};
xLoc = [1:3,5:7,9:11,13:15];
figure();

% Organize data
tempMetric = dendriteBreakProp;
groupData = [];
xGroup = [];
keyPlot = [];
xkey = [];
xExpGroup = [];

for i = 1:size(tempMetric,2)

    % Organize groups into experimental groups
    dMetric = [];

    % Place into metrics
    dMetric = 100-100.*cell2mat(tempMetric(:,i)); % place into matrix/vector

    % Check how data was inputed into cell array
    if size(dMetric,2)>1
        dMetric = dMetric';
        dMetric = reshape(dMetric,[],1);% only necessary for breakage
    end
    xMetric = ones(length(dMetric),1)*xLoc(i);
    
    % Organize groups into experiments
    if ismember(i,1)
        exp = 1;

    elseif ismember(i,2)
        exp = 2;

    elseif ismember(i,3)
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
cMap = {'k',[.5 0.5 0.5],[0 0 .75]};
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
xticks([1,2,3]);
xticklabels({'Pre-Cold','Control','Cold Shocked'})
xlabel('Experimental Group','FontSize',30)
ylabel('Dendrite Length (µm)','FontSize',30)
ylim([0 100])
        
% ANOVA
[~, ~, stats] = anova1(groupData,xGroup);
results = multcompare(stats,'CType','bonferroni');

%% Metric adjustment
% Adjust metrics (if necessary)
MasterFeatureTable.BreakLength = (1-MasterFeatureTable.BreakLength).*100;
MasterFeatureTable.DendriteLength = MasterFeatureTable.DendriteLength.*pxlSize;
MasterFeatureTable.Area = MasterFeatureTable.Area.*pxlSize.^2;

% Put into new table to permit filtering and/or clustering
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
clusterTable.Perimeter = MasterFeatureTable.Perimeter;
clusterTable.MeanIntensity = MasterFeatureTable.MeanIntensity;
clusterTable.DendriteIntensity = MasterFeatureTable.DendriteIntensity;
clusterTable.BreakLength = MasterFeatureTable.BreakLength;
clusterTable.MaxFeretDiameter = MasterFeatureTable.MaxFeretDiameter;
clusterTable.MinFeretDiameter = MasterFeatureTable.MinFeretDiameter;

% Add in location of feature into table
clusterTable.BlebLocation = MasterFeatureTable.Centroid(:,1)./(MasterFeatureTable.DendriteLength./pxlSize);
clusterTable.BlebLocation(clusterTable.BlebLocation>1)=1; % Deal with short dendrites for now
MasterFeatureTable.BlebLocation = MasterFeatureTable.Centroid(:,1)./(MasterFeatureTable.DendriteLength./pxlSize);
MasterFeatureTable.BlebLocation(MasterFeatureTable.BlebLocation>1)=1;

% Make everything correct data type for clustering
clusterTable = struct2table(clusterTable);
cluster = table2cell(clusterTable);
clusterMatrix = cellfun(@double, cluster, 'uni', false);
clusterMatrix = cell2mat(clusterMatrix);

% Filter table to remove 'small' features
filterCluster = clusterTable;
filterCluster(filterCluster.MinFeretDiameter < minWidth,:) = [];
filterCluster(filterCluster.Area < minArea,:) = [];
filtMasterTable = MasterFeatureTable;
filtMasterTable(filtMasterTable.MinFeretDiameter < minWidth,:) = [];
filtMasterTable(filtMasterTable.Area < minArea,:) = [];
filterClusterCell = table2cell(filterCluster);
filterClusterMatrix = cellfun(@double, filterClusterCell, 'uni', false);
filterClusterMatrix = cell2mat(filterClusterMatrix);

%% VISUALIZE FEATURE ONLY DATA (Area, Location, etc)
close all

% Define experiments
groupLabel = {'Control';'Control Aged';'Cold Shocked'};
cMap = {'k',[.5 0.5 0.5],[0 0 .75]};

% Initialize
clustTable = filtMasterTable;
k = 0;
groupData = [];
xGroup = [];
xLoc = [1:3];

figure();
% For loop of experiments 
for i = 1:3 % THESE ARE
    k=0;  
    xMet = [];
    dMet = [];
    k=k+1;
    tempIdx = [i];
    expIdx= find(clustTable.Folder==tempIdx(1));
    expData = clustTable(expIdx,:);
    expData = 1-expData.BlebLocation;  %%% METRIC SELECTION %%%

    % Swarm plot gather
    dMet = [dMet;expData];
    xMet = [xMet;ones(length(expData),1)*xLoc(tempIdx(1))];

    % Swarm chat with all of the data
    swarmchart(xMet,dMet,50,cMap{i},'filled','MarkerFaceAlpha',0.95)
    hold on 

    % Gather for ANOVA
    groupData = [groupData;dMet];
    xGroup = [xGroup;xMet];
end

% Overlay box chart
cMap = {'k',[.5 .5 0.5],[0 0 .5],'c','k',[.5 0 0],[0 0 .5],'c','k',[.5 0 0],[0 0 .5],'c','k',[.5 0 0],[0 0 .5],'c'};
for i =1:max(xGroup)
        idx = find(xGroup==i);
        x = xGroup(idx);
        y = groupData(idx);
        boxchart(x,y,'MarkerStyle','none','BoxFaceColor',cMap{i},'LineWidth',4,'WhiskerLineColor',cMap{i})
end

% Make pretty
title('Normalized Location','fontsize',36)
ax = gca;
ax.FontSize = 18;
xticks(1:3)
xticklabels(groupLabel)
ylabel('Normalized Location','FontSize',20)
ylim([0 1])
% yticks([0 1])
xlabel('Exp Group','FontSize',20)

% Adjust size
set(gcf,'Position',[2 71 1105 906]); % total 

% Now do one way ANOVA
[~, ~, stats] = anova1(groupData,xGroup);
results = multcompare(stats,'CType','bonferroni');

%% Feature per dendrite metrics

% Blebs per dendrite length
close all

% Define experiments
groupLabel = {'Control';'Control Aged';'Cold Shocked'};
cMap = {'k',[.5 0.5 0.5],[0 0 .75]};

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
xCount = [1:3];

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
    tempIdx = i;

    maxImgIdx= find(clustTable.Folder==tempIdx(1));
    maxImgTable = clustTable(maxImgIdx,:);
    maxImg = max(maxImgTable.ImageNumber);

    % Third for loop for image
    for img = 1:maxImg
        for dend = 1:4
            dendIdx = dendIdx+1;
            ftSpread = [];
            tempClustIdx = find(clustTable.Dendrite == dend & filtMasterTable.Folder == tempIdx(1) & filtMasterTable.ImageNumber == img);
            expData = clustTable(tempClustIdx,:);
            expCount = size(expData,1);

            dendMetric = expData.DendriteLength*pxlSize; % This is where you choose what you want to look at
            
            % Feature or bleb dependent
            if expCount > 0
                dendLength = expData.DendriteLength(1)*pxlSize;
                dendBreak = (expData.BreakLength(1));
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
            weightFeature = expCount./dendBreak.*100;

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
            xGroupMetric = [xGroupMetric;ones(length(weightFeature),1)*xCount(tempIdx(1))]; 
            
            % Organize data for plotting
            dMet = [dMet;ftPerDendLen(img,dend)];
            xMet = [xMet;ones(length(ftPerDendLen(img,dend)),1).*xCount(tempIdx(1))];              
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

% Overlay box chart
cMap = {'k',[.5 0.5 0.5],[0 0 .75]};
for i =1:max(xGroupPlot)
        idx = find(xGroupPlot==i);
        x = xGroupPlot(idx);
        y = groupDataPlot(idx);
        boxchart(x,y,'MarkerStyle','none','BoxFaceColor',cMap{i},'LineWidth',4,'WhiskerLineColor',cMap{i})
end

% Make pretty    
title('Feature Count per Dendrite Length','fontsize',36)
ax = gca;
ax.FontSize = 18;
xticks([1:3])
xticklabels(groupLabel)
ylabel('Feature Count (#/µm)','FontSize',20)
xlabel('Experiment Group','FontSize',20)
ylim([0 .15])

% Adjust size
set(gcf,'Position',[2 71 1105 906]);

% Now do one way ANOVA over for loop
[~, ~, stats] = anova1(groupDataPlot,xGroupPlot);
results = multcompare(stats,'CType','bonferroni');


