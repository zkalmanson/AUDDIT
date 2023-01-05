% Code used to analyze and obtain metrics from feature tables produced by
% the neurodegeneration master code.
% Date: 06-13-2022
% Updated: 08-30-2022


%% LOAD AND PREPARE DATA
clc;
clear;
close all

% Find the pixel size to get actual measurements
pxlSize = 0.1032; 
pxlSize = pxlSize./4;% Divided by 4 because of the resize
minWidth = 0;
minArea = 0; 

% Load data
load('rot-FeatureTable.mat');
load('rot-CroppedImages.mat');
load('rot-BreakTable.mat');


% Get width info
for i = 1:size(dendriteWidths,2)
    folderWidth = dendriteWidths(:,i);
    for ii = 1:size(dendriteWidths(:,i),1)
        tempWidth = cell2mat(folderWidth(ii));
        tempWidth(tempWidth==0) = NaN; % try to remove breaks from thickness calculations
        medWidth{ii,i} = nanmedian(tempWidth,1)';
        meanWidth{ii,i} = mean(tempWidth,1,'omitnan')';
    end
end

% Width info of only the controls
for c = 1:3
    controlWidth = cell2mat(medWidth(:,c));
    controlBreak = cell2mat(dendriteBreakProp(:,c)')';
    rmvIdx = find(controlBreak > 0.1); 
    controlWidth(rmvIdx) = []; % Removes dendrites with too many breaks
    refWidth{c} = controlWidth;
    controlMedWidth(c) = median(controlWidth);
    controlMeanWidth(c) = mean(controlWidth);
end


%% Dendrite only Data
close all

% Initialization for figure production
figure();

% Organize data
tempMetric = dendriteLengths; %%% METRIC SELECTION %%%
groupData = [];
xGroup = [];
keyPlot = [];
xkey = [];

for i = 1:size(tempMetric,2)
    % Organize groups
    if i <= 3
        k = 1;
    elseif i <= 6
        k = 2;
    elseif i <= 9
        k = 3;
    end
    dMetric = [];

    % Place into metrics
    dMetric = pxlSize.*cell2mat(tempMetric(:,i)); % place into matrix/vector
    
    % Check how data was inputed into cell array
    if size(dMetric,2)>1
        dMetric = dMetric';
        dMetric = reshape(dMetric,[],1);% only necessary for breakage
    end
    xMetric = ones(length(dMetric),1)*k;
    
    % Overall swarm plot
    swarmchart(xMetric,dMetric,50,'k','filled','MarkerFaceAlpha',0.95)
    hold on
    
    % Highlight key points
    for j = 1:size(dMetric,1)
        img = ceil(j/4);
        % Use this to make critical points larger
        if i == 3 && img == 6 
            keyPlot = [keyPlot;dMetric(j)];
            xkey = [xkey;xMetric(j)];

        elseif i == 5 && img == 1
            keyPlot = [keyPlot;dMetric(j)];
            xkey = [xkey;xMetric(j)];

        elseif i == 7 && img == 12
            keyPlot = [keyPlot;dMetric(j)];
            xkey = [xkey;xMetric(j)];

        end
    end

    % Place into ANOVA/boxchart
    groupData = [groupData;dMetric];
    xGroup = [xGroup;xMetric];

end

% Overlay box chart
cMap = {[0 .25 0],[.5 .5 0],[0.25 0 0]};
for i =1:max(xGroup)
    idx = find(xGroup==i);
    x = xGroup(idx);
    y = groupData(idx);
    boxchart(x,y,'MarkerStyle','none','BoxFaceColor',cMap{i},'LineWidth',4,'WhiskerLineColor',cMap{i})
end

% Overlay key points and change symbol by dendrite
d1 = 1:4:9;
d2 = 2:4:10;
d3 = 3:4:11;
d4 = 4:4:12;

for symbol = 1:12
    if ismember(symbol,d1) == 1
        if symbol < 5
            swarmchart(xkey(symbol)-.25,keyPlot(symbol),1500,'pg','filled','MarkerFaceAlpha',1,'MarkerEdgeColor','k','LineWidth',4)
        elseif symbol < 9 && symbol > 4
            swarmchart(xkey(symbol)-.25,keyPlot(symbol),1500,'py','filled','MarkerFaceAlpha',1,'MarkerEdgeColor','k','LineWidth',4)
        elseif symbol > 8 
            swarmchart(xkey(symbol)-.25,keyPlot(symbol),1500,'pr','filled','MarkerFaceAlpha',1,'MarkerEdgeColor','k','LineWidth',4)
        end
    elseif ismember(symbol,d2) == 1
        if symbol < 5
            swarmchart(xkey(symbol)-.1,keyPlot(symbol),1000,'sg','filled','MarkerFaceAlpha',1,'MarkerEdgeColor','k','LineWidth',4)
        elseif symbol < 9 && symbol > 4
            swarmchart(xkey(symbol)-.1,keyPlot(symbol),1000,'sy','filled','MarkerFaceAlpha',1,'MarkerEdgeColor','k','LineWidth',4)
        elseif symbol > 8 
            swarmchart(xkey(symbol)-.1,keyPlot(symbol),1000,'sr','filled','MarkerFaceAlpha',1,'MarkerEdgeColor','k','LineWidth',4)
        end
    elseif ismember(symbol,d3) == 1
        if symbol < 5
            swarmchart(xkey(symbol)+.1,keyPlot(symbol),1000,'^g','filled','MarkerFaceAlpha',1,'MarkerEdgeColor','k','LineWidth',4)
        elseif symbol < 9 && symbol > 4
            swarmchart(xkey(symbol)+.1,keyPlot(symbol),1000,'^y','filled','MarkerFaceAlpha',1,'MarkerEdgeColor','k','LineWidth',4)
        elseif symbol > 8 
            swarmchart(xkey(symbol)+.1,keyPlot(symbol),1000,'^r','filled','MarkerFaceAlpha',1,'MarkerEdgeColor','k','LineWidth',4)
        end
    elseif ismember(symbol,d4) == 1
        if symbol < 5
            swarmchart(xkey(symbol)+.25,keyPlot(symbol),1000,'og','filled','MarkerFaceAlpha',1,'MarkerEdgeColor','k','LineWidth',4)
        elseif symbol < 9 && symbol > 4
            swarmchart(xkey(symbol)+.25,keyPlot(symbol),1000,'oy','filled','MarkerFaceAlpha',1,'MarkerEdgeColor','k','LineWidth',4)
        elseif symbol > 8 
            swarmchart(xkey(symbol)+.25,keyPlot(symbol),1000,'or','filled','MarkerFaceAlpha',1,'MarkerEdgeColor','k','LineWidth',4)
        end
    end
end

% Adjust visualizations
ax=gca;
ax.FontSize = 30;
title('Dendrite Lengths','FontSize',40);
xticks([1:max(xGroup)]);
xticklabels({'Control','0.03 µM Rotenone','0.5 µM Rotenone'})
xlabel('Experimental Group','FontSize',30)
ylabel('Dendrite Breakage (%)','FontSize',30)
ylim([0 100])

% Adjust size
set(gcf,'Position',[1369 71 1105 906]);

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
clusterTable.BlebLocation(clusterTable.BlebLocation>1)=1;
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
% Initialize
metric = 'Area'; %%% METRIC SELECTION %%%
plotTitle = 'Min Feret Diameter';
groupData = [];
xGroup = [];
clustTable = filtMasterTable;
k = 0;
figure();

% For loop of experiments 
for i = 1:3 % Which experiment to run
    cMap = {'k',[.5 .5 .5],[.25 .25 .25],'k',[.5 .5 .5],[.25 .25 .25],'k',[.5 .5 .5],[.25 .25 .25]};
    k=0;
    for j = 1:3 % Which group
        k=k+1;
        tempIdx = [i,i+3,i+6];
        expIdx= find(clustTable.Folder==tempIdx(j));
        expData = clustTable(expIdx,:);
        expData = expData.(metric);  %%% METRIC SELECTION %%%
  
        % Gather for ANOVA/bar chart
        groupData = [groupData;expData];
        xGroup = [xGroup;ones(length(expData),1)*k]; 

        % Plotting data - swarm chart
        groupPlot = expData;
        xPlot = ones(length(expData),1)*k;
        swarmchart(xPlot,groupPlot,25,cMap{i},'filled','MarkerFaceAlpha',0.5)
        hold on
    end
end

% Overlay box plot
cMap = {[0 .25 0],[.5 .5 0],[0.25 0 0]};
for i =1:max(xGroup)
    idx = find(xGroup==i);
    x = xGroup(idx);
    y = groupData(idx);
    boxchart(x,y,'MarkerStyle','none','BoxFaceColor',cMap{i},'LineWidth',4,'WhiskerLineColor',cMap{i})
end

% Adjust visualizations
ax=gca;
ax.FontSize = 30;
title(plotTitle,'FontSize',40);
xticks([1:max(xGroup)]);
xticklabels({'Control','0.03 µM Rotenone','0.5 µM Rotenone'})
xlabel('Experimental Group','FontSize',30)
ylabel(plotTitle,'FontSize',30)
% ylim([0 100])

% Adjust size
set(gcf,'Position',[1369 71 1105 906]);

% ANOVA
[~, ~, stats] = anova1(groupData,xGroup);
results = multcompare(stats,'CType','bonferroni');

%% Get feature per dendrite results

% Blebs per dendrite length
close all
% Define experiments
exp1 = [1,4,7];
exp2 = [2,5,8];
exp3 = [3,6,9];
totExp = [exp1;exp2;exp3];
groupLabel = {'Control';'0.03 uM';'0.5 uM'};

% Initialize
groupData = [];
xGroup = [];
groupDataLen = [];
xGroupLen = [];
groupDataSpread = [];
xGroupSpread = [];
groupDataMetric = [];
xGroupMetric = [];
groupDataW = [];
xGroupW = [];
keyPlot = [];
xkey = [];
normMet = [];
normX = [];
fldr = 0;
metric = 'MinFeretDiameter'; %%% METRIC SELECTION %%%
clustTable = filtMasterTable;
figure();

% For loop of experiments 
for i = 1:3 % WHICH EXPERIMENTS TO RUN
    k=0;
    cMap = {'k',[.5 .5 .5],[.25 .25 .25]};
    for j = 1:3 % WHICH GROUPS TO RUN
        k=k+1;
        % Initialize
        ftPerDend = [];
        ftPerDendLen = [];
        ftSpread = [];
        tempMetric = [];
        tempXmetric = [];
        tempIdx = [i,i+3,i+6];
        maxImgIdx= find(clustTable.Folder==tempIdx(j));
        maxImgTable = clustTable(maxImgIdx,:);
        maxImg = max(maxImgTable.ImageNumber);
        groupPlot = [];
        xPlot = [];

        fldr = fldr+1;

        % Third for loop for image
        for img = 1:maxImg
%         for img = 3
            for dend = 1:4
                ftSpread = [];
                tempClustIdx = find(clustTable.Dendrite == dend & filtMasterTable.Folder == tempIdx(j) & filtMasterTable.ImageNumber == img);
                expData = clustTable(tempClustIdx,:);
                expCount = size(expData,1);
                dendMetric = expData.(metric); % This is where you choose what you want to look at
                
                % Get dendrite width of current dendrite
                tempWidth = cell2mat(medWidth(img,tempIdx(j)));
                tempWidth = tempWidth(dend);
                if tempWidth <= 1
                    tempWidth = 1;
                end
                
                % Normalize metric to dendrite width, if wanted
                normMetric = dendMetric./controlMeanWidth(i);
                normMet = [normMet;normMetric];
                normX = [normX;ones(length(normMetric),1).*k];
                

                % Feature or bleb dependent
                if expCount > 0
                    dendLength = expData.DendriteLength(1);
                    dendBreak = (expData.BreakLength(1))./100;
                    totDendLength = dendLength.*dendBreak; % blebs per dendrite length actually there

                    % Get spread of features only if more than 4
                    if expCount > 4
                        ftSpread = std(expData.BlebLocation)*pxlSize;
                    end

                else
                    % If no dendrite detected, remove div by 0 issue
                    dendLength = 1*pxlSize;
                    totDendLength = 1;
                    dendBreak = 1;
                end
                
                % Feature count per dend
                ftPerDend(img,dend) = expCount;

                % Feature count per dendrite length
                ftPerDendLen(img,dend) = expCount./totDendLength;

                % Weight feature count by dendrite remaining
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
                groupDataMetric = [groupDataMetric;dendMetric];
                xGroupMetric = [xGroupMetric;ones(length(dendMetric),1)*k]; 

                % Gather data for swarm chart
                groupPlot = [groupPlot;normMet];
                xPlot = [xPlot;normX];

                    % Highlight key points
                    % Use this to make critical points larger
%                     if tempIdx(j) == 3 && img == 6
%                         keyPlot = [keyPlot;ftPerDendLen(img,dend)];
%                         xkey = [xkey;ones(length(ftPerDendLen(img,dend)),1)*1];
% 
%                     elseif tempIdx(j) == 5 && img == 1
%                         keyPlot = [keyPlot;ftPerDendLen(img,dend)];
%                         xkey = [xkey;ones(length(ftPerDendLen(img,dend)),1)*2];
% 
%                     elseif tempIdx(j) == 7 && img == 12
%                         keyPlot = [keyPlot;ftPerDendLen(img,dend)];
%                         xkey = [xkey;ones(length(ftPerDendLen(img,dend)),1)*3];
%                     end
                
            end
        end

        % Plotting data - swarm chart
        swarmchart(normX,normMet,50,cMap{i},'filled','MarkerFaceAlpha',0.15)
        hold on

    end
end

% Overlay box plot
cMap = {[0 .25 0],[.5 .5 0],[0.25 0 0]};
for i =1:max(normX)
    idx = find(normX==i);
    x = normX(idx);
    y = normMet(idx);
    boxchart(x,y,'MarkerStyle','none','BoxFaceColor',cMap{i},'LineWidth',4,'WhiskerLineColor',cMap{i})
end

ax = gca;
% Make pretty
title('Normalized Feature Width','fontsize',40)
ax = gca;
ax.FontSize = 30;
xticks(1:max(xGroupMetric))
xticklabels(groupLabel)
ylabel('Feature Size to Dendrite Width','FontSize',30)
xlabel('Experiment Group','FontSize',30)
% ylim([0 10])

% Adjust size
set(gcf,'Position',[1369 71 1105 906]);

% ANOVA
[~, ~, stats] = anova1(normMet,normX);
results = multcompare(stats,'CType','bonferroni');

