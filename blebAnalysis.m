% Function used to place features with their respective dendrites and
% create an over neurodegeneration feature table.

% Updated: 08-30-2022

function [blebs1data,blebs2data,blebs3data,blebs4data,blebsNoneData] = ...
    blebAnalysis(n1,n2,n3,n4,imFeatures,sFactor,totBreakProp,dLength,imOrigBig,dendIntensity)

% Combine dendrite tracks
n = n1 | n2 | n3 | n4;

% Remove small 'blebs'
imFeatures = bwareaopen(imFeatures,4*ceil(sFactor));

% Try to place blebs - get bleb pixel list
stats = regionprops('table',imFeatures,'all');
b1 = stats.PixelIdxList;

% Initialize 
bleb1 = zeros(size(imFeatures));
bleb2 = zeros(size(imFeatures));
bleb3 = zeros(size(imFeatures));
bleb4 = zeros(size(imFeatures));
blebNone = zeros(size(imFeatures));
blebLoc = zeros(size(b1,1),1);
blebSplit = zeros(size(b1,1),1);

        % For loop to place each feature with its respective dendrite
        for ii = 1:size(b1,1)
            bleb = b1{ii};
            bw1 = zeros(size(imFeatures));
            tempbw = bw1;
            tempbw(bleb) = 1; 
        
            % Check if bleb fits in any dendrite track
            tempLoc(1) = nnz(tempbw&n1);
            tempLoc(2) = nnz(tempbw&n2);
            tempLoc(3) = nnz(tempbw&n3);
            tempLoc(4) = nnz(tempbw&n4);
            tempLocNone = nnz(tempbw&~n); % Location of features not on tracks
            
            % Check to make sure not more than one track
            [m,idx] = max(tempLoc);
            if m == 0
                blebLoc(ii) = 0;
            elseif length(find(tempLoc)) > 1
        
                % Find which track has the majority (>75%) of the feature
                tempAprop = tempLoc./sum(tempLoc);
                tempAIdx = find(tempAprop>0.75, 1);
        
                if ~isempty(tempAIdx)
                    blebLoc(ii) = tempAIdx;
                else
                    tempTrack = b1{ii};
                    tempImg = zeros(size(imFeatures));
                    tempOlap = zeros([size(imFeatures),4]);
                    tempImg(tempTrack) = 1;
        
                    % Locations of overlaps
                    tempOlap(:,:,1) = tempImg&n1;
                    tempOlap(:,:,2) = tempImg&n2;
                    tempOlap(:,:,3) = tempImg&n3;
                    tempOlap(:,:,4) = tempImg&n4;
                    tempImgIdx = tempOlap(:,:,idx);
        
                    blebLoc(ii) = idx; % Index for needing to be split  
                    blebSplit(ii) = 1;
                end
            % It found the right track
            elseif length(find(tempLoc)) == 1
                blebLoc(ii) = idx;
            end
            
            % Create labeled feature images
            if blebLoc(ii) == 1
                bleb1(bleb) = 1;
            elseif blebLoc(ii) == 2
                bleb2(bleb) = 1;
            elseif blebLoc(ii) == 3
                bleb3(bleb) = 1;
            elseif blebLoc(ii) == 4
                bleb4(bleb) = 1;
            elseif blebLoc(ii) == 0
                blebNone(bleb) = 1;
            end
        end
    % Make into binary images
    bleb1 = logical(bleb1);
    bleb2 = logical(bleb2);
    bleb3 = logical(bleb3);
    bleb4 = logical(bleb4);
    blebNone = bwareaopen(logical(blebNone),2*round(sFactor));
    
    % Get regionprops
    blebs1data = regionprops('table',bleb1,imOrigBig,'all');
    blebs2data = regionprops('table',bleb2,imOrigBig,'all');
    blebs3data = regionprops('table',bleb3,imOrigBig,'all');
    blebs4data = regionprops('table',bleb4,imOrigBig,'all');
    blebsNoneData = regionprops('table',blebNone,imOrigBig,'all');
    
    % In case they are empty add the other columns
    
    % Add bleb number so number doesn't get lost
    blebs1data.BlebNumber = [1:size(blebs1data,1)]';
    blebs2data.BlebNumber = [1:size(blebs2data,1)]';
    blebs3data.BlebNumber = [1:size(blebs3data,1)]';
    blebs4data.BlebNumber = [1:size(blebs4data,1)]';
    blebsNoneData.BlebNumber = [1:size(blebsNoneData,1)]';
    
    % Add in which dendrite and its break/length
    blebs1data.BreakLength = totBreakProp(1).*ones(size(blebs1data,1),1);
    blebs2data.BreakLength = totBreakProp(2).*ones(size(blebs2data,1),1);
    blebs3data.BreakLength = totBreakProp(3).*ones(size(blebs3data,1),1);
    blebs4data.BreakLength = totBreakProp(4).*ones(size(blebs4data,1),1);
    
    % Add dendrite length to the feature table
    blebs1data.DendriteLength = dLength(1).*ones(size(blebs1data,1),1);
    blebs2data.DendriteLength = dLength(2).*ones(size(blebs2data,1),1);
    blebs3data.DendriteLength = dLength(3).*ones(size(blebs3data,1),1);
    blebs4data.DendriteLength = dLength(4).*ones(size(blebs4data,1),1);
    
    % Label which dendrite the feature was on
    blebs1data.Dendrite = 1.*ones(size(blebs1data,1),1);
    blebs2data.Dendrite = 2.*ones(size(blebs2data,1),1);
    blebs3data.Dendrite = 3.*ones(size(blebs3data,1),1);
    blebs4data.Dendrite = 4.*ones(size(blebs4data,1),1);
    
    % Label the dendrite intensity the feature was one
    blebs1data.DendriteIntensity = dendIntensity(1).*ones(size(blebs1data,1),1);
    blebs2data.DendriteIntensity = dendIntensity(2).*ones(size(blebs2data,1),1);
    blebs3data.DendriteIntensity = dendIntensity(3).*ones(size(blebs3data,1),1);
    blebs4data.DendriteIntensity = dendIntensity(4).*ones(size(blebs4data,1),1);

end