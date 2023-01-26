function [nuer1mask,nuer2mask,nuer3mask,nuer4mask,dendInterp,n1,n2,n3,n4,dendWidth] = trackDendrites(dendCrop,imFeatures,sFactor)
% Function to track dendrites from image. Takes in cropped dendrite image
% and tracks the four individual dendrites for future analysis

% Updated: 08/30/2022


%% Initialization for speed
close all
max_pt = zeros(size(dendCrop));
max_loc = zeros(size(dendCrop,1),4);
max_pic = zeros(size(dendCrop));

split = 5*ceil(sFactor*4); % Number of sections to split tracking into 
meanrange = 5*ceil(sFactor*4); % Number of sections to split tracking into (foor moving mean calculations)
minprom = 20; % Min prominence needed for a pixel to be seen as a local maxima point
obsRng = 3*ceil(sFactor*2); % How wide should the tracks be


%% Step 1: Find the 4 max intensity values for each row/y-value
    
    for ii = 1:size(dendCrop,1)
        % Find four local maxima points at each row of pixels
        max_pt(ii,:) = islocalmax(dendCrop(ii,:),'maxnumextrema',4,'minseparation',6*round(sFactor*2), ...
            'minprominence',minprom); 

        % Count how many pixels were detected as maxima
        length_max = length(find(max_pt(ii,:)== 1));
        length_max_up = 5-length_max;

        % If less than 4 max points, make the other index values on first
        % column
        if length_max == 0
            max_loc(ii,:) = [1 1 1 1];
        elseif length_max < 4
            max_loc(ii,1:4-length_max) = 1;
            max_loc(ii,length_max_up:4) = find(max_pt(ii,:) == 1);
        else
        % All 4 max values found, index their column location
            max_loc(ii,:) = find(max_pt(ii,:) == 1);
        end
        
        % Create a binary image of the four local maxima points
        max_pic(ii,max_loc(ii,:)) = 1;
    end
    
    % Create image using just local maxima at each row (removes single
    % points not near other dendrites)
    max_pic = bwareaopen(logical(max_pic(:,2:size(max_pic,2))),2);
    
    %% Step 2: Attempt assigning each index value with correct dendrite
    
    % Create temporary variables to not overwrite original max points
    max_loc2 = max_loc;
    max_loc3 = max_loc;
    
    % Replace unknown local max locations with a guess location. The guess
    % is 1/8, 3/8, 5/8, and 7/8 for the 1-4 dendrites, respectively.
    for i = 1:size(max_loc,1)
        tIdx = find(max_loc2(i,:)==1);
        for j = 1:length(tIdx)
            if tIdx(j) == 1
                max_loc2(i,1) = 1.*round(size(dendCrop,1)./8);
            elseif tIdx(j) == 2
                max_loc2(i,2) = 3.*round(size(dendCrop,1)./8);
            elseif tIdx(j) == 3
                max_loc2(i,3) = 5.*round(size(dendCrop,1)./8);
            elseif tIdx(j) == 4
                max_loc2(i,4) = 7.*round(size(dendCrop,1)./8);
            end

        end
    end
        
    % Find the median location in the whole image of each four dendrites  
    max_loc_avg = round(median(max_loc2,1));
    
    % Creates ranges where you would expect to see each dendrite (1 to 4)
    dend_loc_avg(1) = round(mean(max_loc_avg(1:2)));
    dend_loc_avg(2) = round(mean(max_loc_avg(2:3)));
    dend_loc_avg(3) = round(mean(max_loc_avg(3:4)));
    max_avg_rng(1,:) = [2 dend_loc_avg(1)];
    max_avg_rng(2,:) = [dend_loc_avg(1) dend_loc_avg(2)];
    max_avg_rng(3,:) = [dend_loc_avg(2) dend_loc_avg(3)];
    max_avg_rng(4,:) = [dend_loc_avg(3) size(max_pic,2)];
    
    % Another initialization
    max_pic2 = zeros(size(dendCrop));
    neur_range = zeros(size(max_pic,1),2,4);
    nuer_idx = zeros(size(max_loc));
        
    % Places each of 4 maxima into correct dendrite index
    for kk = 1:size(max_pic,1)

        % No dendrites were detected, place in the average bin
        if range(max_loc3(kk,:)) == 0
            for mm = 1:4
                neur_range(kk,:,mm) = [max_loc_avg(mm)-obsRng max_loc_avg(mm)+obsRng];
            end

        % Less than 4 dendrites were detected. Try to match each index
        % location with its repective dendrite. nuer_idx is a m_rows by
        % 4 column matrix if the pixel location of each dendrite.   
        elseif min(max_loc3(kk,:)) == 1 
            for nn = 1:4
               loc_max_idx = max_loc(kk,nn);
               for oo = 1:4
                   % Determine which dendrite the local max index should be
                   % assigned to.
                   if loc_max_idx < max_avg_rng(oo,2) && loc_max_idx >= max_avg_rng(oo,1)
                       if nuer_idx(kk,oo) == 0
                            nuer_idx(kk,oo) = loc_max_idx;
                       elseif nuer_idx(kk,oo) > loc_max_idx
                           nuer_idx(kk,oo-1) = loc_max_idx;
                       elseif nuer_idx(kk,oo) < loc_max_idx && oo == 1 && nuer_idx(kk,oo+1) == 0
                           nuer_idx(kk,oo+1) = loc_max_idx;
                       elseif nuer_idx(kk,oo) < loc_max_idx && oo == 1 && nuer_idx(kk,oo+1) ~= 0
                           if nuer_idx(kk,oo+2) == 0
                               nuer_idx(kk,oo+2) = loc_max_idx;
                           elseif nuer_idx(kk,oo+3) == 0
                               nuer_idx(kk,oo+3) = loc_max_idx;   
                           end
                       elseif nuer_idx(kk,oo) < loc_max_idx && nuer_idx(kk,oo-1) == 0
                           nuer_idx(kk,oo-1) = nuer_idx(kk,oo);
                           nuer_idx(kk,oo) = loc_max_idx;
                       elseif nuer_idx(kk,oo) < loc_max_idx && nuer_idx(kk,oo-1) ~= 0 && oo > 2
                           nuer_idx(kk,oo-2) = nuer_idx(kk,oo-1);
                           nuer_idx(kk,oo-1) = nuer_idx(kk,oo);
                           nuer_idx(kk,oo) = loc_max_idx;
                       elseif nuer_idx(kk,oo) < loc_max_idx && nuer_idx(kk,oo-1) ~= 0 && oo <= 2 && nuer_idx(kk,oo+1) == 0
                           nuer_idx(kk,oo+1) = loc_max_idx;
                       elseif nuer_idx(kk,oo) < loc_max_idx && nuer_idx(kk,oo-1) ~= 0 && oo <= 2 && nuer_idx(kk,oo+1) ~= 0
                           if nuer_idx(kk,oo+2) == 0
                           nuer_idx(kk,oo+2) = loc_max_idx;
                           elseif nuer_idx(kk,oo+2) ~= 0 && oo == 1
                               nuer_idx(kk,oo+3) = loc_max_idx;
                           end
                       end
                   end
               end
            end
        % All four local maxima were detected
        else
            nuer_idx(kk,:) = max_loc(kk,:);
        end
    end

    % Obtain an updated index average with newly sorted dendrite indicies.
    neur_idx_avg = round([sum(nuer_idx(:,1))./nnz(nuer_idx(:,1)),...
        sum(nuer_idx(:,2))./nnz(nuer_idx(:,2)), ...
        sum(nuer_idx(:,3))./nnz(nuer_idx(:,3)), ...
        sum(nuer_idx(:,4))./nnz(nuer_idx(:,4))]);

    % Replace all instances where tracing location is 0 with the overall
    % average of the dendrite to get moving mean
    nuer_idx3 = nuer_idx;
    for rr = 1:size(nuer_idx,1)
        for tt = 1:4
            if nuer_idx(rr,tt) == 0
                nuer_idx3(rr,tt) = neur_idx_avg(tt);
            end
        end
    end
    neur_mvmean = movmean(nuer_idx3,meanrange);

    % Create binary image to show index locations after first round of
    % filtering
    nuer_idx2 = nuer_idx;
    nuer_idx2(nuer_idx2<=0) = 1;
    nuer_idx_pic = zeros(size(max_pic2));
    for yy = 1:size(max_pic,1)
         nuer_idx_pic(yy,nuer_idx2(yy,:)) = 1;
    end
    nuer_idx_pic = logical(nuer_idx_pic);


    %% Step 3: Veritcally close nuer-idx to obtain rough dendrite image
    neur_close = imclose(nuer_idx_pic,strel('rectangle',[11*round(sFactor)*2 1])); 
    neur_close = neur_close(:,2:end);
    neur_close = bwareaopen(neur_close,4*round(sFactor)*2); 

    % Use 'closed' objects to bin them in respective dendrite locations
    % Get centroids of each object as a way to match objects with
    % dendrites.
    neur_close_props = regionprops(neur_close,'all');
    neur_close_cent = regionprops(neur_close,'centroid');
    close_cent = cell2mat(struct2cell(neur_close_cent));
    x_cent = close_cent(1:2:end);
    y_cent = close_cent(2:2:end);

    % Create range for binning from binned averages in previous section
    neur_loc_avg(1) = round(mean(neur_idx_avg(1:2)));
    neur_loc_avg(2) = round(mean(neur_idx_avg(2:3)));
    neur_loc_avg(3) = round(mean(neur_idx_avg(3:4)));
    neur_avg_rng(1,:) = [2 neur_loc_avg(1)];
    neur_avg_rng(2,:) = [neur_loc_avg(1) neur_loc_avg(2)];
    neur_avg_rng(3,:) = [neur_loc_avg(2) neur_loc_avg(3)];
    neur_avg_rng(4,:) = [neur_loc_avg(3) size(nuer_idx_pic,2)];

    % Initialization for closed object binning
    bin_avg = zeros(split,4);
    bin_loc_avg = zeros(split,3);
    bin_avg_rng = zeros(4,2,split);
    y_bin = zeros(split,2);
    binNan = zeros(split,4);

    % Splits image into smaller y portions to account for overall
    % intensity differences along y.
    img_split = round(size(nuer_idx,1)/split); 

    % For loop to create a controlled moving average for object binning
    for run = 1:split
    
        % Find pixel rows to analyze
        top = (run-1)*img_split+1:run*img_split;
        topover = top>size(nuer_idx,1);
        top(topover) = [];

        % Deals with too large image split and rounding keeps making
        % smaller until it will work
        topCheck = 2;
        while isempty(top) == 1
            top = (run-topCheck)*img_split+1:(run-1)*img_split;
            topover = top>size(nuer_idx,1);
            top(topover) = [];
            topCheck = topCheck+1;
        end
        
        % Find the top and bottom rows of the bins
        y_bin(run,:) = [top(1)-.99,top(length(top))+.01];
        
        % Moving bin averages
        bin_avg(run,:) = round([sum(nuer_idx(top,1))./nnz(nuer_idx(top,1)), sum(nuer_idx(top,2))./nnz(nuer_idx(top,2)), ...
        sum(nuer_idx(top,3))./nnz(nuer_idx(top,3)), sum(nuer_idx(top,4))./nnz(nuer_idx(top,4))]);

        % If NaNs included, find them and replace with previous average
        binNan(run,:) = isnan(bin_avg(run,:));
        for kk = 1:4
           if binNan(run,kk) == 1
               bin_avg(run,kk) = neur_idx_avg(kk);
           end
        end

        % Create range for binning
        bin_loc_avg(run,1) = round(mean(bin_avg(run,1:2)));
        bin_loc_avg(run,2) = round(mean(bin_avg(run,2:3)));
        bin_loc_avg(run,3) = round(mean(bin_avg(run,3:4)));
        bin_avg_rng(1,:,run) = [2 bin_loc_avg(run,1)];
        bin_avg_rng(2,:,run) = [bin_loc_avg(run,1) bin_loc_avg(run,2)];
        bin_avg_rng(3,:,run) = [bin_loc_avg(run,2) bin_loc_avg(run,3)];
        bin_avg_rng(4,:,run) = [bin_loc_avg(run,3) size(nuer_idx_pic,2)];

    end

    % Edit moving mean to get out of outliers (to be used later)
    for rr = 1:size(neur_mvmean,1)
        for zone = 1:split
            if rr < y_bin(zone,2) && rr >= y_bin(zone,1)
                xx = floor(y_bin(zone,1)):ceil(y_bin(zone,2));
                xxunder = xx<1;
                xx(xxunder) = [];
                xxover = xx>size(neur_mvmean,1);
                xx(xxover) = [];
                neur_zone_avg = mean(round(neur_mvmean(xx,:)),1);
                zone_idx = zone;
            end
        end   
        for tt = 1:4
           absdif = abs(neur_mvmean(rr,tt)-neur_zone_avg(tt));
           if absdif > 8*round(sFactor)*2
               neur_mvmean(rr,tt) = neur_zone_avg(tt);
           end
        end
    end

    % Place object in bins based on x value of centroid
    binIdx = zeros(length(x_cent),4);
    for ii = 1:length(x_cent)
        for zone = 1:split
            if y_cent(ii) < y_bin(zone,2) && y_cent(ii) >= y_bin(zone,1)
                zone_idx = zone;
            end
        end
        for jj = 1:4
            if x_cent(ii) < bin_avg_rng(jj,2,zone_idx) && x_cent(ii) >= bin_avg_rng(jj,1,zone_idx)
                binIdx(ii,jj) = ii;
            end
        end 
    end
    
    % Start keeping dendrites indicies separate
    dend1idx = binIdx(:,1);
    dend1idx(dend1idx==0) = [];
    dend2idx = binIdx(:,2);
    dend2idx(dend2idx==0) = [];
    dend3idx = binIdx(:,3);
    dend3idx(dend3idx==0) = [];
    dend4idx = binIdx(:,4);
    dend4idx(dend4idx==0) = [];

    % Create images with binned objects
    dend1close = neur_close;
    dend2close = neur_close;
    dend3close = neur_close;
    dend4close = neur_close;
    
    % Remove pixels from other dendrites
    for uu = 1:length(x_cent)
        if uu ~= dend1idx
           dend1close(neur_close_props(uu).PixelIdxList) = 0;
        end
        if uu ~= dend2idx
            dend2close(neur_close_props(uu).PixelIdxList) = 0;
        end
        if uu ~= dend3idx
            dend3close(neur_close_props(uu).PixelIdxList) = 0;
        end
        if uu ~= dend4idx
            dend4close(neur_close_props(uu).PixelIdxList) = 0;
        end
    end
      
    %% Step 4: Skeletonize dendrites to find center points of each dendrite

    % Begin skeletonizing previous closed dendrite images
    dend1skel = bwskel(dend1close,'MinBranchLength',3*round(sFactor)*2);
    dend2skel = bwskel(dend2close,'MinBranchLength',3*round(sFactor)*2);
    dend3skel = bwskel(dend3close,'MinBranchLength',3*round(sFactor)*2);
    dend4skel = bwskel(dend4close,'MinBranchLength',3*round(sFactor)*2);

    % Label the detected objects in each binary dendrite image
    dend1_lbl = bwlabel(dend1skel);
    dend2_lbl = bwlabel(dend2skel);
    dend3_lbl = bwlabel(dend3skel);
    dend4_lbl = bwlabel(dend4skel);

    % Find the centroids of each binary object for each dendrite
    dend1cent = regionprops(dend1skel,'centroid');
    dend1cent = cell2mat(struct2cell(dend1cent));
    dend1cent = dend1cent(1:2:end);
    dend2cent = regionprops(dend2skel,'centroid');
    dend2cent = cell2mat(struct2cell(dend2cent));
    dend2cent = dend2cent(1:2:end);
    dend3cent = regionprops(dend3skel,'centroid');
    dend3cent = cell2mat(struct2cell(dend3cent));
    dend3cent = dend3cent(1:2:end);
    dend4cent = regionprops(dend4skel,'centroid');
    dend4cent = cell2mat(struct2cell(dend4cent));
    dend4cent = dend4cent(1:2:end);

    %% Step 5: Rebin skeletonized points to ensure they are labeled as the correct dendrite 
    

    % Start with dendrite 1
    %%%% ---DENDRITE 1--- %%%%

    % Initialize
    obj2mov1up = zeros(size(dendCrop));

    % For loop over each pixel row
    for yy = 2:size(dendCrop,1)-1

        % Find how many detected points were in this pixel row
        dend1_temp = dend1skel(yy,:);
        tempidx = find(dend1_temp);

        % If more than one point was found on this row, determine which one
        % should stay labeled as dendrite 1,
        if length(tempidx) > 1
            
            % Determine which zone you are in (which y-section)
            for zone = 1:split
                if yy < y_bin(zone,2) && yy >= y_bin(zone,1)
                    zone_idx = zone;
                end
            end

            % Label object
            tempobj = dend1_lbl(yy,tempidx);
            tempdist = zeros(length(tempidx),1);
            
            % Find the difference between object and previous moving
            % average
            for zz = 1:length(tempidx)
                tempdist(zz) = abs(tempidx(zz)-bin_avg(zone_idx,1));
            end

            % Find and keep object closer to the moving average
            [~, mindist_idx] = min(tempdist);
            real_idx = tempidx(mindist_idx);
            real_obj = tempobj(mindist_idx);
            
            % Determine whether the other object should be placed in
            % another dendrite image.
            for i = 1:length(tempidx)
                tempup = zeros(1,3);
                tempdown = zeros(1,3);
                tempup(1,:) = dend1_lbl(yy-1,tempidx(i)-1:tempidx(i)+1);
                tempdown(1,:) = dend1_lbl(yy+1,tempidx(i)-1:tempidx(i)+1);
                tempup(tempup<=0) = [];
                tempdown(tempdown<=0) = [];
                if isempty(tempup) ~= 1 || isempty(tempdown) ~= 1
                    if isempty(tempup) == 1
                        tempup = 0;
                    end
                    if isempty(tempdown) == 1
                        tempdown = 0;
                    end
                    if tempobj(i) == tempup(1) || tempobj(i) == tempdown(1)
                        temp_cent = dend1cent(tempobj(i));
                        % If the onject is not the real object, and it is
                        % out of bounds of the moving average, label it as
                        % an object to move up (to dendrite 2)
                        if temp_cent < neur_avg_rng(1,1) || temp_cent > neur_avg_rng(1,2)
                            if tempidx(i) ~= real_idx && tempidx(i) > real_idx && tempobj(i) ~= real_obj
                                obj2mov1up(yy,tempidx(i)) = dend1_lbl(yy,tempidx(i));
                                dend1skel(yy,tempidx(i)) = 0;
                            end
                        end

                    end
                elseif tempidx(i) ~= real_idx && tempidx(i) > real_idx && tempobj(i) ~= real_obj
                        obj2mov1up(yy,tempidx(i)) = dend1_lbl(yy,tempidx(i));
                        dend1skel(yy,tempidx(i)) = 0;
                end            

            end
        end
    end
    % Find the rows and columns of the objects to move and add them to the
    % dendrite 2 image.
    [obj1uprow, obj1upcol] = find(obj2mov1up);
    for i = 1:length(obj1uprow)
        dend2skel(obj1uprow(i), obj1upcol(i)) = 1;
    end

    %%%% ---DENDRITE 2--- %%%%  
    % It is the same step, but now objects can be moved back down to
    % dendrite 1.

    % Initialize
    obj2mov2up = zeros(size(dendCrop));
    obj2mov2dn = zeros(size(dendCrop));

    for yy = 2:size(dendCrop,1)-1
        dend2_temp = dend2skel(yy,:);
        tempidx = find(dend2_temp);
        if length(tempidx) > 1
            for zone = 1:split
                if yy < y_bin(zone,2) && yy >= y_bin(zone,1)
                    zone_idx = zone;
                end
            end

            % Label object
            tempobj = dend2_lbl(yy,tempidx);
            tempdist = zeros(length(tempidx),1);

            for zz = 1:length(tempidx)
                tempdist(zz) = abs(tempidx(zz)-bin_avg(zone_idx,2));
            end

            % Determine which object is the real one based on moving
            % average
            [~, mindist_idx] = min(tempdist);
            real_idx = tempidx(mindist_idx);
            real_obj = tempobj(mindist_idx);

            for i = 1:length(tempidx)
                tempup = zeros(1,3);
                tempdown = zeros(1,3);
                tempup(1,:) = dend2_lbl(yy-1,tempidx(i)-1:tempidx(i)+1);
                tempdown(1,:) = dend2_lbl(yy+1,tempidx(i)-1:tempidx(i)+1);
                tempup(tempup<=0) = [];
                tempdown(tempdown<=0) = [];
                if isempty(tempup) ~= 1 || isempty(tempdown) ~= 1
                    if isempty(tempup) == 1
                        tempup = 0;
                    end
                    if isempty(tempdown) == 1
                        tempdown = 0;
                    end
                    if tempobj(i) == tempup(1) || tempobj(i) == tempdown(1)
                        temp_cent = dend2cent(tempobj(i));
                        if temp_cent < neur_avg_rng(2,1) || temp_cent > neur_avg_rng(2,2)
                            if tempidx(i) ~= real_idx && tempidx(i) > real_idx && tempobj(i) ~= real_obj
                                obj2mov2up(yy,tempidx(i)) = dend2_lbl(yy,tempidx(i));
                                dend2skel(yy,tempidx(i)) = 0;
                            elseif tempidx(i) ~= real_idx && tempidx(i) < real_idx && tempobj(i) ~= real_obj
                                obj2mov2dn(yy,tempidx(i)) = dend2_lbl(yy,tempidx(i));
                                dend2skel(yy,tempidx(i)) = 0;
                            end
                        end

                    end

                % Label as object that needs to be moved up to dendrite 3
                elseif tempidx(i) ~= real_idx && tempidx(i) > real_idx && tempobj(i) ~= real_obj
                    obj2mov2up(yy,tempidx(i)) = dend2_lbl(yy,tempidx(i));
                    dend2skel(yy,tempidx(i)) = 0;

                % Label as object that needs to move to dendrite 1
                elseif tempidx(i) ~= real_idx && tempidx(i) < real_idx && tempobj(i) ~= real_obj
                    obj2mov2dn(yy,tempidx(i)) = dend2_lbl(yy,tempidx(i));
                    dend2skel(yy,tempidx(i)) = 0;
                end            
            end
        end
    end

    % Find indicies for objects that need to be moved up/down
    [obj2uprow, obj2upcol] = find(obj2mov2up);
    [obj2dnrow, obj2dncol] = find(obj2mov2dn);

    % Adjust dendrite 1 and/or dendrite 2 images
    for i = 1:length(obj2uprow)
        dend3skel(obj2uprow(i), obj2upcol(i)) = 1;
    end
    for i = 1:length(obj2dnrow)
        dend1skel(obj2dnrow(i), obj2dncol(i)) = 1;
    end

    %%%% ---DENDRITE 3--- %%%%
    % Same as dendrite 2.

    obj2mov3up = zeros(size(dendCrop));
    obj2mov3dn = zeros(size(dendCrop));
    for yy = 2:size(dendCrop,1)-1
        dend3_temp = dend3skel(yy,:);
        tempidx = find(dend3_temp);
        if length(tempidx) > 1
            for zone = 1:split
                if yy < y_bin(zone,2) && yy >= y_bin(zone,1)
                    zone_idx = zone;
                end
            end
            tempobj = dend3_lbl(yy,tempidx);
            tempdist = zeros(length(tempidx),1);
            for zz = 1:length(tempidx)
                tempdist(zz) = abs(tempidx(zz)-bin_avg(zone_idx,3));
            end
                [~, mindist_idx] = min(tempdist);
                real_idx = tempidx(mindist_idx);
                real_obj = tempobj(mindist_idx);
            for i = 1:length(tempidx)
                tempup = zeros(1,3);
                tempdown = zeros(1,3);
                tempup(1,:) = dend3_lbl(yy-1,tempidx(i)-1:tempidx(i)+1);
                tempdown(1,:) = dend3_lbl(yy+1,tempidx(i)-1:tempidx(i)+1);
                tempup(tempup<=0) = [];
                tempdown(tempdown<=0) = [];
                if isempty(tempup) ~= 1 || isempty(tempdown) ~= 1
                    if isempty(tempup) == 1
                        tempup = 0;
                    end
                    if isempty(tempdown) == 1
                        tempdown = 0;
                    end
                    if tempobj(i) == tempup(1) || tempobj(i) == tempdown(1)
                        temp_cent = dend3cent(tempobj(i));
                        if temp_cent < neur_avg_rng(3,1) || temp_cent > neur_avg_rng(3,2)
                            if tempidx(i) ~= real_idx && tempidx(i) > real_idx && tempobj(i) ~= real_obj
                                obj2mov3up(yy,tempidx(i)) = dend3_lbl(yy,tempidx(i));
                                dend3skel(yy,tempidx(i)) = 0;
                            elseif tempidx(i) ~= real_idx && tempidx(i) < real_idx && tempobj(i) ~= real_obj
                                obj2mov3dn(yy,tempidx(i)) = dend3_lbl(yy,tempidx(i));
                                dend3skel(yy,tempidx(i)) = 0;
                            end
                        end

                    end
                elseif tempidx(i) ~= real_idx && tempidx(i) > real_idx && tempobj(i) ~= real_obj
                    obj2mov3up(yy,tempidx(i)) = dend3_lbl(yy,tempidx(i));
                    dend3skel(yy,tempidx(i)) = 0;
                elseif tempidx(i) ~= real_idx && tempidx(i) < real_idx && tempobj(i) ~= real_obj
                    obj2mov3dn(yy,tempidx(i)) = dend3_lbl(yy,tempidx(i));
                    dend3skel(yy,tempidx(i)) = 0;
                end            

            end
        end
    end
    [obj3uprow, obj3upcol] = find(obj2mov3up);
    [obj3dnrow, obj3dncol] = find(obj2mov3dn);
    for i = 1:length(obj3uprow)
        dend4skel(obj3uprow(i), obj3upcol(i)) = 1;
    end
    for i = 1:length(obj3dnrow)
        dend2skel(obj3dnrow(i), obj3dncol(i)) = 1;
    end

    %%%% ---DENDRITE 4--- %%%%
    % Same as dendrite 1 (see above for more comments), except objects can
    % only be sent back down to dendrite 3.

    obj2mov4dn = zeros(size(dendCrop));
    for yy = 2:size(dendCrop,1)-1
        dend4_temp = dend4skel(yy,:);
        tempidx = find(dend4_temp);
        if length(tempidx) > 1
            for zone = 1:split
                if yy < y_bin(zone,2) && yy >= y_bin(zone,1)
                    zone_idx = zone;
                end
            end
            tempobj = dend4_lbl(yy,tempidx);
            tempdist = zeros(length(tempidx),1);
            for zz = 1:length(tempidx)
                tempdist(zz) = abs(tempidx(zz)-bin_avg(zone_idx,4));
            end
                [~, mindist_idx] = min(tempdist);
                real_idx = tempidx(mindist_idx);
                real_obj = tempobj(mindist_idx);
            for i = 1:length(tempidx)  
                tempup = zeros(1,3);
                tempdown = zeros(1,3);
                tempup(1,:) = dend4_lbl(yy-1,tempidx(i)-1:tempidx(i)+1);
                tempdown(1,:) = dend4_lbl(yy+1,tempidx(i)-1:tempidx(i)+1);
                tempup(tempup<=0) = [];
                tempdown(tempdown<=0) = [];
                if isempty(tempup) ~= 1 || isempty(tempdown) ~= 1
                    if isempty(tempup) == 1
                        tempup = 0;
                    end
                    if isempty(tempdown) == 1
                        tempdown = 0;
                    end
                    if tempobj(i) == tempup(1) || tempobj(i) == tempdown(1)
                        temp_cent = dend4cent(tempobj(i));
                        if temp_cent < neur_avg_rng(4,1) || temp_cent > neur_avg_rng(4,2)
                            if tempidx(i) ~= real_idx && tempidx(i) < real_idx && tempobj(i) ~= real_obj
                                obj2mov4dn(yy,tempidx(i)) = dend4_lbl(yy,tempidx(i));
                                dend4skel(yy,tempidx(i)) = 0;
                            end
                        end
                    end
                elseif tempidx(i) ~= real_idx && tempidx(i) < real_idx && tempobj(i) ~= real_obj
                    obj2mov4dn(yy,tempidx(i)) = dend4_lbl(yy,tempidx(i));
                    dend4skel(yy,tempidx(i)) = 0;
                end            

            end
        end
    end

    % Find object to move down and adjust dendrite 3 image.
    [obj4dnrow, obj4dncol] = find(obj2mov4dn);
    for i = 1:length(obj4dnrow)
        dend3skel(obj4dnrow(i), obj4dncol(i)) = 1;
    end


    %% Step 6: Remove multiple points labeled as the same dendrite
    
    % Label objects in updated dendrite images
    dend1_lbl = bwlabel(dend1skel);
    dend2_lbl = bwlabel(dend2skel);
    dend3_lbl = bwlabel(dend3skel);
    dend4_lbl = bwlabel(dend4skel);
    
    % Find size and location of labeled objects for each dendrite
    dend1area = regionprops(dend1skel,'area');
    dend1area = [dend1area.Area];
    [~,dend1max] = max(dend1area);
    dend1cent = regionprops(dend1skel,'centroid');
    dend1cent = cell2mat(struct2cell(dend1cent));
    dend1cent = dend1cent(1:2:end);

    dend2area = regionprops(dend2skel,'area');
    dend2area = [dend2area.Area];
    [~,dend2max] = max(dend2area);
    dend2cent = regionprops(dend2skel,'centroid');
    dend2cent = cell2mat(struct2cell(dend2cent));
    dend2cent = dend2cent(1:2:end);

    dend3area = regionprops(dend3skel,'area');
    dend3area = [dend3area.Area];
    [~,dend3max] = max(dend3area);
    dend3cent = regionprops(dend3skel,'centroid');
    dend3cent = cell2mat(struct2cell(dend3cent));
    dend3cent = dend3cent(1:2:end);

    dend4area = regionprops(dend4skel,'area');
    dend4area = [dend4area.Area];
    [~,dend4max] = max(dend4area);
    dend4cent = regionprops(dend4skel,'centroid');
    dend4cent = cell2mat(struct2cell(dend4cent));
    dend4cent = dend4cent(1:2:end);

    % SinglePointBin Function (Finally decided functions would be better)
    dend1skel = SinglePointBin(dend1skel,dend1_lbl,dend1cent,dend1max,1,neur_mvmean,neur_idx_avg);
    dend2skel = SinglePointBin(dend2skel,dend2_lbl,dend2cent,dend2max,2,neur_mvmean,neur_idx_avg); 
    dend3skel = SinglePointBin(dend3skel,dend3_lbl,dend3cent,dend3max,3,neur_mvmean,neur_idx_avg);
    dend4skel = SinglePointBin(dend4skel,dend4_lbl,dend4cent,dend4max,4,neur_mvmean,neur_idx_avg);


    %% Step 7: Using moving moving mean to remove weird jumps in dendrite skeletons and replace empties with moving mean
    
    % Temp variables to not mess with original dendrite skeletons
    dend1skel2 = dend1skel;
    dend2skel2 = dend2skel;
    dend3skel2 = dend3skel;
    dend4skel2 = dend4skel;

    % Initialize for indexing and speed
    dend1skel_idx = zeros(size(dend1skel2,1),1);
    dend2skel_idx = zeros(size(dend2skel2,1),1);
    dend3skel_idx = zeros(size(dend3skel2,1),1);
    dend4skel_idx = zeros(size(dend4skel2,1),1);
    dend1id = zeros(size(dend1skel2,1),1);
    dend2id = zeros(size(dend2skel2,1),1);
    dend3id = zeros(size(dend3skel2,1),1);
    dend4id = zeros(size(dend4skel2,1),1);

    % Re-check each pixel row for weird jumps and/or empty values
    for rr = 1:size(dend4skel2,1)
        dend1temp = find(dend1skel2(rr,:));

        % Replace empty indicies with moving mean
        if isempty(dend1temp) == 1
            dend1skel_idx(rr) = neur_mvmean(rr,1);
            dend1id(rr) = 0;

        % Keep detected point if only 1 is found
        elseif length(dend1temp) ==1
            dend1skel_idx(rr) = find(dend1skel2(rr,:));
            dend1id(rr) = dend1temp;

        % If more than one found, replace index with the mean value (single
        % point bin failed)
        elseif length(dend1temp) > 1
            dend1skel_idx(rr) = mean(find(dend1skel2(rr,:)));
            dend1id(rr) = mean(dend1temp);
        end

        % Repeat process for dendrite 2-4.

        % Dendrite 2.
        dend2temp = find(dend2skel2(rr,:));
        if isempty(dend2temp) == 1
            dend2skel_idx(rr) = neur_mvmean(rr,2);
            dend2id(rr) = 0;
        elseif length(dend2temp) ==1
            dend2skel_idx(rr) = find(dend2skel2(rr,:));
            dend2id(rr) = dend2temp;
        elseif length(dend2temp) > 1
            dend2skel_idx(rr) = mean(find(dend2skel2(rr,:)));
            dend2id(rr) = mean(dend2temp);
        end
        
        % Dendrite 3.
        dend3temp = find(dend3skel2(rr,:));
        if isempty(dend3temp) == 1
            dend3skel_idx(rr) = neur_mvmean(rr,3);
            dend3id(rr) = 0;
        elseif length(dend3temp) ==1
            dend3skel_idx(rr) = find(dend3skel2(rr,:));
            dend3id(rr) = dend3temp;
        elseif length(dend3temp) > 1
            dend3skel_idx(rr) = mean(find(dend3skel2(rr,:)));
            dend3id(rr) = mean(dend3temp);
        end

        % Dendrite 4.
        dend4temp = find(dend4skel2(rr,:));
        if isempty(dend4temp) == 1
            dend4skel_idx(rr) = neur_mvmean(rr,4);
            dend4id(rr) = 0;
        elseif length(dend4temp) ==1
            dend4skel_idx(rr) = find(dend4skel2(rr,:));
            dend4id(rr) = dend4temp;
        elseif length(dend4temp) > 1
            dend4skel_idx(rr) = mean(find(dend4skel2(rr,:)));
            dend4id(rr) = mean(dend4temp);
        end
        
        % Update NaN and 0 values with an index of 1 for each dendrite
        dend1skel_idx(isnan(dend1skel_idx)) = 1;
        dend1skel_idx(dend1skel_idx==0) = 1;
        dend2skel_idx(isnan(dend2skel_idx)) = 1;
        dend2skel_idx(dend2skel_idx==0) = 1;
        dend3skel_idx(isnan(dend3skel_idx)) = 1;
        dend3skel_idx(dend3skel_idx==0) = 1;
        dend4skel_idx(isnan(dend4skel_idx)) = 1;
        dend4skel_idx(dend4skel_idx==0) = 1;

        % Make new skeleton with updated index values from removing extras
        dend1skel2(rr,ceil(dend1skel_idx(rr))) = 1;
        dend2skel2(rr,ceil(dend2skel_idx(rr))) = 1;
        dend3skel2(rr,ceil(dend3skel_idx(rr))) = 1;
        dend4skel2(rr,ceil(dend4skel_idx(rr))) = 1;

    end

    % Create new images with updated skeleton values
    dend1 = zeros(size(dend1skel));
    dend2 = zeros(size(dend2skel));
    dend3 = zeros(size(dend3skel));
    dend4 = zeros(size(dend4skel));
    dend1id(dend1id<=1) = 1;
    dend2id(dend2id<=1) = 1;
    dend3id(dend3id<=1) = 1;
    dend4id(dend4id<=1) = 1;

    for rr = 1:length(dend1id)
        dend1(rr,round(dend1id(rr))) = 1;
        dend2(rr,round(dend2id(rr))) = 1;
        dend3(rr,round(dend3id(rr))) = 1;
        dend4(rr,round(dend4id(rr))) = 1;    
    end

    %% Step 8: Smooth skeleton plots to remove noisy tracks 

    [dend1smooth, dend1smoothidx] = NeurSmoothData(dend1id,dend1skel2,sFactor); 
    [dend2smooth, dend2smoothidx] = NeurSmoothData(dend2id,dend2skel2,sFactor);
    [dend3smooth, dend3smoothidx] = NeurSmoothData(dend3id,dend3skel2,sFactor);
    [dend4smooth, dend4smoothidx] = NeurSmoothData(dend4id,dend4skel2,sFactor);

    %% Step 9: Interpolate tracks to any remaining empty pixel rows
    
    % Initialization
    dend1_idx = zeros(size(dendCrop,1),1);
    dend2_idx = zeros(size(dendCrop,1),1);
    dend3_idx = zeros(size(dendCrop,1),1);
    dend4_idx = zeros(size(dendCrop,1),1);

    % Index each dendrite to find empty rows
    for xx = 1:size(dendCrop,1)
        if isempty(find(dend1smooth(xx,:),1)) == 1
            dend1_idx(xx) = 0;
        else
            dend1_idx(xx) = find(dend1smooth(xx,:));
        end
        if isempty(find(dend2smooth(xx,:),1)) == 1
            dend2_idx(xx) = 0;
        else
            dend2_idx(xx) = find(dend2smooth(xx,:));
        end    
        if isempty(find(dend3smooth(xx,:),1)) == 1
            dend3_idx(xx) = 0;
        else
            dend3_idx(xx) = find(dend3smooth(xx,:));
        end   
        if isempty(find(dend4smooth(xx,:),1)) == 1
            dend4_idx(xx) = 0;
        else
            dend4_idx(xx) = find(dend4smooth(xx,:));
        end    
    end
    idx1 = find(dend1_idx);
    idx2 = find(dend2_idx);
    idx3 = find(dend3_idx);
    idx4 = find(dend4_idx);
    
    %%%%--- DENDRITE 1 INTERPOLATION---%%%%
    % Double check to make sure index is not empty
    if isempty(idx1)
        idx1 = 1;
        neur1 = dend1_idx;
        neur1(isnan(neur1)) = 0;
        neur1(neur1==0) = ceil(size(dendCrop,2)./4);
        neur(:,1) = neur1;
    else    
        neur1 = dend1_idx;
        nz_neur1 = find(neur1~=0);
        nz_input1 = neur1(nz_neur1);
        if idx1(1) > length(neur1)/2 
            neur1 = interp1(nz_neur1,nz_input1,1:length(neur1),'linear',dend1smoothidx(idx1(1)));
        elseif idx1(end) < length(neur1)/2
            neur1 = interp1(nz_neur1,nz_input1,1:length(neur1),'linear',dend1smoothidx(idx1(end)));
        else
            neur1 = interp1(nz_neur1,nz_input1,1:length(neur1),'linear','extrap');
        end    
        neur1(isnan(neur1)) = 0;
        neur(:,1) = neur1;
    end
    
    %%%%--- DENDRITE 2 INTERPOLATION---%%%%
    % Double check to make sure index is not empty
    if isempty(idx2)
        idx2 = 1;
        neur2 = dend2_idx;
        neur2(isnan(neur2)) = 0;
        neur2(neur2==0) = ceil(size(dendCrop,2)./4)*2;
        neur(:,2) = neur2;
    else
        neur2 = dend2_idx;
        nz_neur2 = find(neur2~=0);
        nz_input2 = neur2(nz_neur2);
        if idx2(1) > length(neur2)/2 
            neur2 = interp1(nz_neur2,nz_input2,1:length(neur2),'linear',dend2smoothidx(idx2(1)));
        elseif idx2(end) < length(neur2)/2
            neur2 = interp1(nz_neur2,nz_input2,1:length(neur2),'linear',dend2smoothidx(idx2(end)));
        else
            neur2 = interp1(nz_neur2,nz_input2,1:length(neur2),'linear','extrap');
        end    
        neur2(isnan(neur2)) = 0;
        neur(:,2) = neur2;            
    end
    
    %%%%--- DENDRITE 3 INTERPOLATION---%%%%  
    if isempty(idx3)
         idx3 = 1;
         neur3 = dend3_idx;
         neur3(isnan(neur3)) = 0;
         neur3(neur3==0) = ceil(size(dendCrop,2)./4)*3;
         neur(:,3) = neur3;
    else
        neur3 = dend3_idx;
        nz_neur3 = find(neur3~=0);
        nz_input3 = neur3(nz_neur3);
        if idx3(1) > length(neur3)/2 
            neur3 = interp1(nz_neur3,nz_input3,1:length(neur3),'linear',dend3smoothidx(idx3(1)));
        elseif idx3(end) < length(neur3)/2
            neur3 = interp1(nz_neur3,nz_input3,1:length(neur3),'linear',dend3smoothidx(idx3(end)));
        else
            neur3 = interp1(nz_neur3,nz_input3,1:length(neur3),'linear','extrap');
        end    
        neur3(isnan(neur3)) = 0;
        neur(:,3) = neur3;
    end
    
    %%%%--- DENDRITE 4 INTERPOLATION---%%%%
    if isempty(idx4)
         idx4 = 1;
         neur4 = dend4_idx;
         neur4(isnan(neur4)) = 0;
         neur4(neur4==0) = floor(size(dendCrop,2)./4)*4-1;
         neur(:,4) = neur4;
    else
        neur4 = dend4_idx;
        nz_neur4 = find(neur4~=0);
        nz_input4 = neur4(nz_neur4);
        if idx4(1) > length(neur4)/2 
            neur4 = interp1(nz_neur4,nz_input4,1:length(neur4),'linear',dend4smoothidx(idx4(1)));
        elseif idx4(end) < length(neur4)/2
            neur4 = interp1(nz_neur4,nz_input4,1:length(neur4),'linear',dend4smoothidx(idx4(end)));
        else
            neur4 = interp1(nz_neur4,nz_input4,1:length(neur4),'linear','extrap');
        end
        neur4(isnan(neur4)) = 0;
        neur(:,4) = neur4;
    end
    
    % Combine and round indicies for final tracks   
    neur = round(neur);

    %Create dendrite ranges
    neur_range = zeros(size(nuer_idx,1),2,4);
    obsRng = 15*ceil(sFactor/2); 
    for zz = 1:size(nuer_idx,1)
        for mm = 1:4
            neur_range(zz,:,mm) = [neur(zz,mm)-obsRng neur(zz,mm)+obsRng];
        end
    end
    neur_range(neur_range<=0) = 1;
    neur_range(neur_range>size(dendCrop,2)) = size(dendCrop,2);

    % Create mask images of tracks with observation ranges

    % Initialization
    nuer1mask = zeros(size(dendCrop));
    nuer2mask = zeros(size(dendCrop));
    nuer3mask = zeros(size(dendCrop));
    nuer4mask = zeros(size(dendCrop));
    
    % Create images
    for ll = 1:size(dendCrop,1)
        nuer1mask(ll,neur_range(ll,1,1):neur_range(ll,2,1)) = 1;
        nuer2mask(ll,neur_range(ll,1,2):neur_range(ll,2,2)) = 1;
        nuer3mask(ll,neur_range(ll,1,3):neur_range(ll,2,3)) = 1;
        nuer4mask(ll,neur_range(ll,1,4):neur_range(ll,2,4)) = 1;
    end
    
    % Total track image
    nuermask = nuer1mask | nuer2mask | nuer3mask | nuer4mask;
    
    % Show final mask
%     imshowpair(dendCrop,nuermask)

    %% Step 10: Remove overlapping tracks

    % Find which rows were interpolated (for overlap reasons)
    dend1Interp = find(dend1_idx==0);
    dend2Interp = find(dend2_idx==0);
    dend3Interp = find(dend3_idx==0);
    dend4Interp = find(dend4_idx==0);
    dendInterp = {length(dend1Interp); length(dend2Interp); length(dend3Interp); length(dend4Interp)};

    % Create a mask that finds area of overlap
    oLap = nuer1mask+nuer2mask+nuer3mask+nuer4mask;
    oLapLoc = find(oLap>1);

    % Create a coded overlap system
    x1 = 1.*nuer1mask;
    x2 = 2.*nuer2mask;
    x3 = 5.*nuer3mask;
    x4 = 10.*nuer4mask;
    x = x1+x2+x3+x4; % coded overlap

    % For loop that looks at each possible way of mask overlapping. If the
    % overlap region contains a dendrite that was interpolated, the
    % interpolated dendrite will be the one with the smaller track.

    for qq = 1:length(oLapLoc)
        tempLap = x(oLapLoc(qq));
        [row, ~] = ind2sub(size(nuermask),oLapLoc(qq));
        if tempLap == 3 % 1+2
            if isempty(find(dend1Interp==row, 1)) ~= 1 && isempty(find(dend2Interp==row, 1)) == 1
                x1(oLapLoc(qq)) = 0;
            elseif isempty(find(dend1Interp==row, 1)) == 1 && isempty(find(dend2Interp==row, 1)) ~= 1
                x2(oLapLoc(qq)) = 0;
            else 
                x2(oLapLoc(qq)) = 0;
            end
        elseif tempLap == 6 %1+3
            if isempty(find(dend1Interp==row, 1)) ~= 1 && isempty(find(dend3Interp==row, 1)) == 1
                x1(oLapLoc(qq)) = 0;
            elseif isempty(find(dend1Interp==row, 1)) == 1 && isempty(find(dend3Interp==row, 1)) ~= 1
                x3(oLapLoc(qq)) = 0;
            else 
                x3(oLapLoc(qq)) = 0;
            end
        elseif tempLap == 11 % 1+4
            if isempty(find(dend1Interp==row, 1)) ~= 1 && isempty(find(dend4Interp==row, 1)) == 1
                x1(oLapLoc(qq)) = 0;
            elseif isempty(find(dend1Interp==row, 1)) == 1 && isempty(find(dend4Interp==row, 1)) ~= 1
                x4(oLapLoc(qq)) = 0;
            else 
                x4(oLapLoc(qq)) = 0;
            end
        elseif tempLap == 7 %2+3
            if isempty(find(dend2Interp==row, 1)) ~= 1 && isempty(find(dend3Interp==row, 1)) == 1
                x2(oLapLoc(qq)) = 0;
            elseif isempty(find(dend2Interp==row, 1)) == 1 && isempty(find(dend3Interp==row, 1)) ~= 1
                x3(oLapLoc(qq)) = 0;
            else 
                x3(oLapLoc(qq)) = 0;
            end
        elseif tempLap == 12 %2+4
            if isempty(find(dend2Interp==row, 1)) ~= 1 && isempty(find(dend4Interp==row, 1)) == 1
                x2(oLapLoc(qq)) = 0;
            elseif isempty(find(dend2Interp==row, 1)) == 1 && isempty(find(dend4Interp==row, 1)) ~= 1
                x4(oLapLoc(qq)) = 0;
            else 
                x4(oLapLoc(qq)) = 0;
            end
        elseif tempLap == 15 %3+4
            if isempty(find(dend3Interp==row, 1)) ~= 1 && isempty(find(dend4Interp==row, 1)) == 1
                x3(oLapLoc(qq)) = 0;
            elseif isempty(find(dend3Interp==row, 1)) == 1 && isempty(find(dend4Interp==row, 1)) ~= 1
                x4(oLapLoc(qq)) = 0;
            else 
                x4(oLapLoc(qq)) = 0;
            end
        elseif tempLap == 8 %1+2+3
            if isempty(find(dend1Interp==row, 1)) ~= 1 && isempty(find(dend2Interp==row, 1)) ~= 1 && isempty(find(dend3Interp==row, 1)) == 1
                x1(oLapLoc(qq)) = 0;
                x2(oLapLoc(qq)) = 0;
            elseif isempty(find(dend1Interp==row, 1)) == 1 && isempty(find(dend2Interp==row, 1)) ~= 1 && isempty(find(dend3Interp==row, 1)) ~= 1
                x2(oLapLoc(qq)) = 0;
                x3(oLapLoc(qq)) = 0;
            elseif isempty(find(dend1Interp==row, 1)) ~= 1 && isempty(find(dend2Interp==row, 1)) == 1 && isempty(find(dend3Interp==row, 1)) ~= 1
                x1(oLapLoc(qq)) = 0;
                x3(oLapLoc(qq)) = 0;
            else 
                x2(oLapLoc(qq)) = 0;
                x3(oLapLoc(qq)) = 0;
            end
        elseif tempLap == 17 %2+3+4
            if isempty(find(dend4Interp==row, 1)) ~= 1 && isempty(find(dend2Interp==row, 1)) ~= 1 && isempty(find(dend3Interp==row, 1)) == 1
                x4(oLapLoc(qq)) = 0;
                x2(oLapLoc(qq)) = 0;
            elseif isempty(find(dend4Interp==row, 1)) == 1 && isempty(find(dend2Interp==row, 1)) ~= 1 && isempty(find(dend3Interp==row, 1)) ~= 1
                x2(oLapLoc(qq)) = 0;
                x3(oLapLoc(qq)) = 0;
            elseif isempty(find(dend4Interp==row, 1)) ~= 1 && isempty(find(dend2Interp==row, 1)) == 1 && isempty(find(dend3Interp==row, 1)) ~= 1
                x4(oLapLoc(qq)) = 0;
                x3(oLapLoc(qq)) = 0;
            else 
                x2(oLapLoc(qq)) = 0;
                x3(oLapLoc(qq)) = 0;
            end
        elseif tempLap == 16 % 1+3+4
            if isempty(find(dend1Interp==row, 1)) ~= 1 && isempty(find(dend4Interp==row, 1)) ~= 1 && isempty(find(dend3Interp==row, 1)) == 1
                x1(oLapLoc(qq)) = 0;
                x4(oLapLoc(qq)) = 0;
            elseif isempty(find(dend1Interp==row, 1)) == 1 && isempty(find(dend4Interp==row, 1)) ~= 1 && isempty(find(dend3Interp==row, 1)) ~= 1
                x4(oLapLoc(qq)) = 0;
                x3(oLapLoc(qq)) = 0;
            elseif isempty(find(dend1Interp==row, 1)) ~= 1 && isempty(find(dend4Interp==row, 1)) == 1 && isempty(find(dend3Interp==row, 1)) ~= 1
                x1(oLapLoc(qq)) = 0;
                x3(oLapLoc(qq)) = 0;
            else 
                x4(oLapLoc(qq)) = 0;
                x3(oLapLoc(qq)) = 0;
            end
        elseif tempLap == 13 % 1+2+4
            if isempty(find(dend1Interp==row, 1)) ~= 1 && isempty(find(dend2Interp==row, 1)) ~= 1 && isempty(find(dend4Interp==row, 1)) == 1
                x1(oLapLoc(qq)) = 0;
                x2(oLapLoc(qq)) = 0;
            elseif isempty(find(dend1Interp==row, 1)) == 1 && isempty(find(dend2Interp==row, 1)) ~= 1 && isempty(find(dend4Interp==row, 1)) ~= 1
                x2(oLapLoc(qq)) = 0;
                x4(oLapLoc(qq)) = 0;
            elseif isempty(find(dend1Interp==row, 1)) ~= 1 && isempty(find(dend2Interp==row, 1)) == 1 && isempty(find(dend4Interp==row, 1)) ~= 1
                x1(oLapLoc(qq)) = 0;
                x4(oLapLoc(qq)) = 0;
            else 
                x2(oLapLoc(qq)) = 0;
                x4(oLapLoc(qq)) = 0;
            end
        elseif tempLap == 18 % all overlap
                x2(oLapLoc(qq)) = 0;
                x3(oLapLoc(qq)) = 0;
                x4(oLapLoc(qq)) = 0;
        end
    end
    
    % Make into logical image
    n1 = logical(x1);
    n2 = logical(x2);
    n3 = logical(x3);
    n4 = logical(x4);
    n = n1 | n2 | n3 | n4;
    
    %% Step 11: Mask out dendrite images with tracks to obtain dendrite widths 
    close all

    % Mask out cropped image
    dend1img = bsxfun(@times, dendCrop, cast(n1, 'like', dendCrop)); % mask out image
    dend2img = bsxfun(@times, dendCrop, cast(n2, 'like', dendCrop));
    dend3img = bsxfun(@times, dendCrop, cast(n3, 'like', dendCrop));
    dend4img = bsxfun(@times, dendCrop, cast(n4, 'like', dendCrop));
    dendimg = bsxfun(@times, dendCrop, cast(n, 'like', dendCrop));
    
    % Get dendrite widths
    width1 = findDendriteWidth(dendimg,dend1img,imFeatures,sFactor);
    width2 = findDendriteWidth(dendimg,dend2img,imFeatures,sFactor);
    width3 = findDendriteWidth(dendimg,dend3img,imFeatures,sFactor);
    width4  = findDendriteWidth(dendimg,dend4img,imFeatures,sFactor);
    dendWidth = [width1,width2,width3,width4];

end