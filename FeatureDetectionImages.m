% Show feature detection with rectangles
% Andrew Clark
% 04-29-2022
% Updated: 08-30-2022

clc;
clear;
close all

% Load all the feature data from the rotenone experiments
load('cold-FeatureTable.mat');
load('cold-CroppedImages.mat');

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

% Make a matrix for the data
clusterTable = struct2table(clusterTable);
cluster = table2cell(clusterTable);
clusterMatrix = cellfun(@double, cluster, 'uni', false);
clusterMatrix = cell2mat(clusterMatrix);

% Try filtering, if wanted 
filterCluster = clusterTable;
filterCluster(filterCluster.Area < 1,:) = [];
filtMasterTable = MasterFeatureTable;
filtMasterTable(filtMasterTable.Area < 1,:) = [];
filterClusterCell = table2cell(filterCluster);
filterClusterMatrix = cellfun(@double, filterClusterCell, 'uni', false);
filterClusterMatrix = cell2mat(filterClusterMatrix);

% Assume only one cluster so make every feature in cluster 1
filterCluster.Cluster = 1;

%% Initialization and user input
close all

% INSERT PIXEL SIZE FOR SCALING
pxlSize = 0.1032; % FOR SAN MIGUEL CAMERA (0.1032) MEYER CAMERA: 0.38
refPxlSize = 0.38;
sFactor = refPxlSize/pxlSize;

myDir = dir(pwd);
myDirIdx = [myDir.isdir;];
myDir = myDir(myDirIdx);
folderNames = {myDir.name};
folderNames = folderNames(~ismember(folderNames,{'.','..','.DS_Store','Before Scale','beforeFix','Feature Tables','Not Used'}));


%% NOW VISUALIZE
% FOR LOOP OVER NUMBER OF IMAGES

% Which folder is the image in   
for fldr = 2
    cd(folderNames{fldr}) % change directory to folder of interest
    close all
    parent = pwd;
    image_names = dir(fullfile(parent,'*.tif')); 
    image_names = {image_names.name}; % Sorts files by number (ex: 1,2,3,...10,11)
    image_names = natsortfiles(image_names);

    % Which image to look at
    for img_num = 2
        close all
        %% SORT BLEBS BY CATEGORY
        figure();
        i = croppedImages{img_num,fldr};
        imshow(croppedImages{img_num,fldr});
        hold on
        cMap = {'r',[.78 0 0 ],'c','m',[.25 0.25 .75],'g','c','y','m',[1,0,1]};
        
        % Loop over every feature
        for ii = 1:featureLength
            tempBlebLoc = [];
            featureTypeIdx{ii} = find(filtMasterTable.Folder == fldr & filtMasterTable.ImageNumber == img_num & filterCluster.Cluster == ii);
            
            % Find feature bounding box
            for jj = 1:length(featureTypeIdx{ii})
                tempBlebPixel = filtMasterTable.PixelList;
                tempBlebLoc = filtMasterTable.Centroid;
                tempBlebLoc = tempBlebLoc(featureTypeIdx{ii}(jj),:);
                tempBlebLocBB = filtMasterTable.BoundingBox;
                tempBlebLocBB = tempBlebLocBB(featureTypeIdx{ii}(jj),:);

                % Draw shape over detected feature
                i = insertShape(i,"Rectangle",tempBlebLocBB,"LineWidth",15,"Color","yellow");
            end
            
        end

        % NO TRACKED BLEBS
        noBlebIdx = find(blebsNoneTable.Folder == fldr & blebsNoneTable.ImageNumber == img_num);
        if isempty(noBlebIdx)
            
        else
            for nn = 1:length(noBlebIdx)
                tempNoBleb = blebsNoneTable.BoundingBox;
                tempNoBleb = tempNoBleb(noBlebIdx,:);
                rectangle('Position',tempNoBleb(nn,:),'LineWidth',2,'EdgeColor','r')
            end
        end
        figure(); imshow(i)

        % Export bleb categories to image for later analysis
        cd ..
%         f = gcf;
%         exportgraphics(f,'coldDetection.pdf','Resolution',150,'BackgroundColor','none','Append',true)
        cd(folderNames{fldr})
    end         
    cd .. % Change back to main folder    
end

