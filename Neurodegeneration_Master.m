
%% MAIN MASTER: NEURODEGENERATION CODE
% Written: ASC + ZK
% Date: 2022-04-27
% Updated: 2022-08-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization and user input
close all
clear;
clc;

% Add path of code location (CHANGE TO WHERE IMAGES ARE)
addpath(pwd)

% INSERT PIXEL SIZE FOR SCALING
pxlSize = 0.1032; % San Miguel Lab: 0.1032 µm; Meyer Lab: 0.38 µm

% Reference pixel size for scaling purposes
refPxlSize = 0.38;
sFactor = refPxlSize/pxlSize;

% Finds the folder names in the directory
myDir = dir(pwd);
myDirIdx = [myDir.isdir;];
myDir = myDir(myDirIdx);
folderNames = {myDir.name};

% Removes unwanted folders from the directory
folderNames = folderNames(~ismember(folderNames,{'.','..','.DS_Store',...
    'Before Scale','beforeFix','Feature Tables','Not Used','Overlap','Wrong','New Folder'}));

% Inititalize feature tables
MasterFeatureTable = [];
blebsNoneTable = [];

% For loop to analyze each folder of images
for fldr = 1:length(folderNames)
% for fldr = 5

    cd(folderNames{fldr}) % change directory to folder of interest

    % Loads images
    parent = pwd;
    image_names = dir(fullfile(parent,'*.tif')); 
    image_names = {image_names.name}; % Sorts files by number (ex: 1,2,3,...10,11)
    image_names = natsortfiles(image_names);

    % Initialize table variables
    blebTable = [];
    reDoImg = zeros(size(image_names,2),1);

    % Loop over total images in folder
    for img_num = 1:length(image_names)
%     for img_num = 2

        close all
        catchImg = 0; % Preallocate error check
        
        % Reads image
        file=fullfile(parent,image_names{img_num});
        im = (imread(file));

        % Finds dendrites and crops image to just dendrites
        [cropped,imOrig] = findDendrites(im,1,sFactor); %Crop down the image 

        % Change image to 8 bit (from 16 bit)
        cropped2 = im2uint8(cropped); 
        imOrig2 = im2uint8(imOrig);

        % Finds features and returns contrasted adjusted image
        [imLocCont, imFeatures, imDend, imOrigBig] = findBlebs(imOrig2,sFactor);

        % Tracks dendrites without machine learning
        [nuer1mask,nuer2mask,nuer3mask,nuer4mask,dendInterp,n1,n2,n3,n4,dendWidth]...
            = trackDendrites(imLocCont,imFeatures,sFactor); % Returns masks with each tracked dendrite

        % Error check:
        % If too much of image is binarized, there are no dendrites.
        iCheck = imbinarize(imLocCont);
        iFor = length(find(iCheck));
        iTot = size(iCheck,1)*size(iCheck,2);
        iProp = iFor/iTot;
    
        if iProp > .21 % If more than 21% of image is binarized, wrong side
            catchImg = 1;
        end
        if size(cropped,2)*pxlSize < 8 % If the image is less than 8 µm wide, wrong side
            catchImg = 1;
        end
        if size(cropped,1)*pxlSize < 30 % If image is less than 30 µm long, wrong side
            catchImg = 1;
        end
        
        % Wrong side of dendrite was detected, rerun with other side
        if catchImg == 1
            reDoImg(img_num) = 1; % keep track of which image was wrong
            [cropped,imOrig] = findDendrites(im,0,sFactor); % Choose other side
            cropped2 = im2uint8(cropped); % Make 8-bit
            imOrig2 = im2uint8(imOrig);
            [imLocCont, imFeatures, imDend, imOrigBig] = findBlebs(imOrig2,sFactor); % Find features
            [nuer1mask,nuer2mask,nuer3mask,nuer4mask,dendInterp,n1,n2,n3,n4,dendWidth]...
            = trackDendrites(imLocCont,imFeatures,sFactor); % Track dendrites
        end
        
        % Combine tracks
        n = n1|n2|n3|n4;

        % Find locations of breaks
        [dendBreak,breakImg1,breakImg2,breakImg3,breakImg4,b1start,b2start,b3start,b4start,dendIntensity] = ...
            breakAnalysis(n1,n2,n3,n4,imDend,imOrigBig);    
    
        % Find lengths of breaks
        [break1start,break1len,~] = breakLengths(dendBreak(1,:)); 
        [break2start,break2len,~] = breakLengths(dendBreak(2,:));
        [break3start,break3len,~] = breakLengths(dendBreak(3,:));
        [break4start,break4len,~] = breakLengths(dendBreak(4,:));
        breakLen = {break1len,break2len,break3len,break4len};

        % Get lengths of individual dendrites
        d1length = break1start(end)-b1start;
        d2length = break2start(end)-b2start;
        d3length = break3start(end)-b3start;
        d4length = break4start(end)-b4start;
        dLength = [d1length;d2length;d3length;d4length]; 
        dLength(dLength==0) = 1;
        
        % Overall break proportions
        totBreak = sum(dendBreak,2);
        totBreakPrp = totBreak./dLength;


    
        % Function to gather all degeneration data into a table
        [blebs1data,blebs2data,blebs3data,blebs4data,blebsNoneData] =...
            blebAnalysis(n1,n2,n3,n4,imFeatures,sFactor,totBreakPrp,dLength,imOrigBig,dendIntensity);

        %% Place bleb data into a master table for clustering
        %--- Variable names if empty array ---%
        FeatureNames = {'Area'	'Centroid'	'BoundingBox'	'SubarrayIdx'	'MajorAxisLength'...
            'MinorAxisLength'	'Eccentricity'	'Orientation'	'ConvexHull'	'ConvexImage'...
            'ConvexArea'	'Circularity'	'Image'	'FilledImage'	'FilledArea'	'EulerNumber'...
            'Extrema'	'EquivDiameter'	'Solidity'	'Extent'	'PixelIdxList'	'PixelList'	'Perimeter'...
            'PerimeterOld'	'PixelValues'	'WeightedCentroid'	'MeanIntensity'	'MinIntensity'	'MaxIntensity'...
            'MaxFeretDiameter'	'MaxFeretAngle'	'MaxFeretCoordinates'	'MinFeretDiameter'	'MinFeretAngle'...
            'MinFeretCoordinates'	'BlebNumber' 'BreakLength' 'DendriteLength' 'Dendrite' 'DendriteIntensity'...
        	'ImageNumber'	'Folder' 'ImageSize'};

        FeatureNames2 = {'Area'	'Centroid'	'BoundingBox'	'SubarrayIdx'	'MajorAxisLength'...
            'MinorAxisLength'	'Eccentricity'	'Orientation'	'ConvexHull'	'ConvexImage'...
            'ConvexArea'	'Circularity'	'Image'	'FilledImage'	'FilledArea'	'EulerNumber'...
            'Extrema'	'EquivDiameter'	'Solidity'	'Extent'	'PixelIdxList'	'PixelList'	'Perimeter'...
            'PerimeterOld'	'PixelValues'	'WeightedCentroid'	'MeanIntensity'	'MinIntensity'	'MaxIntensity'...
            'MaxFeretDiameter'	'MaxFeretAngle'	'MaxFeretCoordinates'	'MinFeretDiameter'	'MinFeretAngle'...
            'MinFeretCoordinates'	'BlebNumber'	'ImageNumber'	'Folder' 'ImageSize'};

        % Add image number to table
        blebs1data.ImageNumber = img_num.*ones(size(blebs1data,1),1);
        blebs2data.ImageNumber = img_num.*ones(size(blebs2data,1),1);
        blebs3data.ImageNumber = img_num.*ones(size(blebs3data,1),1);
        blebs4data.ImageNumber = img_num.*ones(size(blebs4data,1),1);
        blebsNoneData.ImageNumber = img_num.*ones(size(blebsNoneData,1),1);

        % Add folder to table
        blebs1data.Folder = fldr.*ones(size(blebs1data,1),1);
        blebs2data.Folder = fldr.*ones(size(blebs2data,1),1);
        blebs3data.Folder = fldr.*ones(size(blebs3data,1),1);
        blebs4data.Folder = fldr.*ones(size(blebs4data,1),1);
        blebsNoneData.Folder = fldr.*ones(size(blebsNoneData,1),1);
        
        % Add image size to table
        blebs1data.ImageSize = size(cropped,1).*ones(size(blebs1data,1),1);
        blebs2data.ImageSize = size(cropped,1).*ones(size(blebs2data,1),1);
        blebs3data.ImageSize = size(cropped,1).*ones(size(blebs3data,1),1);
        blebs4data.ImageSize = size(cropped,1).*ones(size(blebs4data,1),1);
        blebsNoneData.ImageSize = size(cropped,1).*ones(size(blebsNoneData,1),1);

        % Correct if part of the table is empty
        if isempty(blebs1data)
            blebs1data = array2table(zeros(0,43));
            blebs1data.Properties.VariableNames = FeatureNames;
        end
        if isempty(blebs2data)
            blebs2data = array2table(zeros(0,43));
            blebs2data.Properties.VariableNames = FeatureNames;
        end
        if isempty(blebs3data)
            blebs3data = array2table(zeros(0,43));
            blebs3data.Properties.VariableNames = FeatureNames;
        end
        if isempty(blebs4data)
            blebs4data = array2table(zeros(0,43));
            blebs4data.Properties.VariableNames = FeatureNames;
        end
        if isempty(blebsNoneData)
            blebsNoneData = array2table(zeros(0,39));
            blebsNoneData.Properties.VariableNames = FeatureNames2;
        end

        % Combine feature data for later analysis
        MasterFeatureTable = [MasterFeatureTable;blebs1data;blebs2data;blebs3data;blebs4data];
        blebsNoneTable = [blebsNoneTable;blebsNoneData];

        % Dendrite only results (not in table)     
        dendriteBreakProp{img_num,fldr} = totBreakPrp';
        dendriteBreakSeqs{img_num,fldr} = breakLen;
        dendriteLengths{img_num,fldr} = dLength;
        dendriteIntensities{img_num,fldr} = dendIntensity';
        dendriteWidths{img_num,fldr} = dendWidth;
        
        % Image data
        imgSize{img_num,fldr} = size(cropped,1); 
        croppedImages{img_num,fldr} = imLocCont;
        croppedImagesOrig{img_num,fldr} = imOrigBig;
    
        %% SAVE PDF IMAGES FOR VIZUAL INSPECTION - if wanted
        close all
        cd ..
    
        %---- BLEB IMAGES ----%
%         figure();
%         imshow(imOrigBig)
%         imshowpair(imLocCont,imFeatures)
%         title(sprintf('Folder: %d Image: %d',fldr,img_num),'FontSize',26)
%         f = gcf;
%         exportgraphics(f,sprintf('rotenoneFeatures.pdf'),'Resolution',100,'Append',true)
    
     
        %--- BREAK IMAGES ----%
%         breakImg = breakImg1|breakImg2|breakImg3|breakImg4;
%         figure();
%         imshowpair(imLocCont,breakImg)
%         title(sprintf('Folder: %d Image: %d',fldr,img_num),'FontSize',26)
%         f = gcf;
%         exportgraphics(f,sprintf('rotenoneBreaks.pdf'),'Resolution',100,'Append',true)
     
    
        %--------- TRACK IMAGES ---------%
%         figure();
%         imshowpair(imLocCont,n)
%         title(sprintf('Folder: %d Image: %d',fldr,img_num),'FontSize',26)
%         f = gcf;
%         exportgraphics(f,sprintf('rotenoneTracks.pdf'),'Resolution',100,'Append',true) 
%     
    
        % Change Back
        cd(folderNames{fldr})

    end
    
    cd .. % Change back to main folder
end 

% Save data as .mat file for later analysis (add what is wanted)

% Feature results
% save('coldFeatureTable.mat','MasterFeatureTable','blebsNoneTable');

% Images
% save('rotenone_CroppedImages.mat','croppedImages'); 

% Dendrite info
% save('coldBreakTable.mat','dendriteBreakProp','dendriteBreakSeqs','dendriteLengths','imgSize','dendriteIntensities','dendriteWidths');
