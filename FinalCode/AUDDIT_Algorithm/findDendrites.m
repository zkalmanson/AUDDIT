function [im_dend,imOrig] = findDendrites(img,side,sFactor)
% Purpose: Imports images of neurons and crops to the size that has
% dendrites for future processing (tracking and degradation
% quantification). User chooses which number image gets sent into the
% function.

% Updated: 08-29-2022

%% Step 1: Load and rotate images
% Binarize image and rotate to make dendrites vertical

% Check to see if it is a color image then make grayscale
if size(img,3) == 3
    img = img(:,:,2);
end
% Create a new image that gets cropped with the adjusted one
imOrig = img;

% Check to see if the image needs to be adjusted (8v16bit) - different
% contrast adjustment used
if isa(img,'uint16')
    im = imadjust(img);
    im2 = imadjust(im);
elseif isa(img,'uint8')
    im2 = img;
end

% Wide-range binarization to rotate worms vertically
imbw = imbinarize(im2,0.1);
prps = regionprops(imbw,'Area','Orientation');
img_prps(:,1) = [prps.Area];
img_prps(:,2) = [prps.Orientation];
[~, max_idx] = max(img_prps(:,1)); % Use only largest object
max_agl = img_prps(max_idx,2); % Find rotation angle
imrot  = imrotate(im2,90-max_agl,'nearest'); % contrasted
imOrig = imrotate(imOrig,90-max_agl,'nearest'); % original image gets same treatment


%% Step 2: Split image into 2
% Cut image into 2 using cell body as reference point

% Binarization of only cell body to get 'center' points
imbw2 = imbinarize(imrot,0.9); % Only get the cell body
imbwe = imerode(imbw2,strel('line',3*round(sFactor),0)); % Remove extra, larger blebs
imbw2 = imdilate(imbwe,strel('rectangle',[7*round(sFactor*.5),13*round(sFactor*.5)])); % Connects cell bodies together
rot_prps = regionprops(bwareaopen(imbw2,10*round(sFactor)),'All'); % Gets regionprops data from the binarized image 
rot_area(:,1) = [rot_prps.Area];
[~, max_idx2] = max(rot_area(:,1));
rot_cell = struct2cell(rot_prps);
rot_bb = cell2mat(rot_cell(3,max_idx2)); % Finds the cell body bounding box
rot_center = cell2mat(rot_cell(2,max_idx2));


% Larger binarization to get bounding box of whole neuron
imbw3 = imclearborder(imbinarize(imrot,0.4));
rot_bb_prps = regionprops(bwareaopen(imbw3,10*round(sFactor)),'Area','BoundingBox'); 
rot3_cell = struct2cell(rot_bb_prps);
rot2_area(:,1) = rot3_cell(1,:);
[~, max_idx3] = max(cell2mat(rot2_area(:,1)));
x_rot_bb = cell2mat(rot3_cell(2,max_idx3)); 

% Find the cropping locations from cell body and overall bounding box
x_range = [rot_center(1)-1.25*x_rot_bb(3) rot_center(1)+1.25*x_rot_bb(3)]; % Uses 1.25x the x length of the neruon bounding box
y_range = [rot_bb(2) rot_bb(2)+1*rot_bb(4)]; % bounding box of the cell body
x_mag = x_range(2)-x_range(1); % total magnitude change of bounding box
top_y = 1; % First 'row'
bot_y = length(imrot(:,1)); % Last 'row'

% Now create 2 images from the rotated image
imtop = imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]);
imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);  % (2nd point used to be center_y)use y range to get rid of cell body/ center_y to include cell body


%% Step 3: Obtain information from each side

% Binarize each image to obtain dendrite 'objects' -- skip first and last
% rows to avoid edge artifacts from rotating and cropping
T_top = adaptthresh(imtop(5:size(imtop,1)-(4),5:size(imtop,2)-(4)),0.3); 
T_bot = adaptthresh(imbot(5:size(imbot,1)-(4),5:size(imbot,2)-(4)),0.3);

top_bw = bwareaopen(imbinarize(imtop(5:size(imtop,1)-(4),5:size(imtop,2)-(4)),T_top),400); %SCALE BROKE IT 
bot_bw = bwareaopen(imbinarize(imbot(5:size(imbot,1)-(4),5:size(imbot,2)-(4)),T_bot),400);

% figure()
% imshowpair(top_bw,bot_bw,'montage')

% Obtain regionprops data for each image
top_prps = regionprops(top_bw,'all');
bot_prps = regionprops(bot_bw,'all');

% If no objects recognized - give values of 0 for object size/orientation
if isempty(top_prps)
    top_stats(:,1) = 0;
    top_stats(:,2) = 0;
    top_stats(:,3) = 0;
    top_stats(:,4) = 0;
    
else
    top_stats(:,1) = [top_prps.Area];
    top_stats(:,2) = [top_prps.Orientation];
    top_stats(:,3) = [top_prps.Circularity];
    top_stats(:,4) = [top_prps.MajorAxisLength];

end
if isempty(bot_prps)
    bot_stats(:,1) = 0;
    bot_stats(:,2) = 0; 
    bot_stats(:,3) = 0;
    bot_stats(:,4) = 0;

else
    bot_stats(:,1) = [bot_prps.Area];
    bot_stats(:,2) = [bot_prps.Orientation];
    bot_stats(:,3) = [bot_prps.Circularity];
    bot_stats(:,4) = [bot_prps.MajorAxisLength];
end

% Keep the four largest objects from each side
[~, top_idx] = maxk(top_stats(:,1),4);
[~, bot_idx] = maxk(bot_stats(:,1),4);

% Obtain the orientation of the 4 largest objects
top_ang = top_stats(top_idx,2);
bot_ang = bot_stats(bot_idx,2);

topLen = top_stats(top_idx,4);
botLen = bot_stats(bot_idx,4);

% Obtain stats about the distribution of objects
top_std = std(abs(top_ang));
bot_std = std(abs(bot_ang));
topSize = size(top_bw,1);
botSize = size(bot_bw,1);

% Obtain axis lengths of the 4 largest
topAxis = top_stats(top_idx,4);
botAxis = bot_stats(bot_idx,4);

% Proportional sizes
topSizeProp = topSize/(topSize+botSize);
botSizeProp = botSize/(topSize+botSize);
topAxisProp = median(topAxis)/(median(topAxis)+median(botAxis));
botAxisProp = median(botAxis)/(median(topAxis)+median(botAxis));

%% Step 4: If statement to determine which side is dendrites
% Uses size and angle data to chose which side has the dendrites
% If nothing binarized, other side must be the dendrites

% First check circularity
if abs(max(top_stats(:,3))-max(bot_stats(:,3))) > 0.2 && max(top_stats(:,3)) < max(bot_stats(:,3)) && median(topLen) > median(botLen) && topSizeProp > .65
    imtop = imrotate(imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
    imOrig2 = imrotate(imcrop(imOrig,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
    im_dend = imtop; 
elseif abs(max(top_stats(:,3))-max(bot_stats(:,3))) > 0.2 && max(top_stats(:,3)) > max(bot_stats(:,3)) && median(topLen) < median(botLen) && botSizeProp > .65
    imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
    im_dend = imbot;  
    imOrig2 = imcrop(imOrig, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
elseif size(top_prps,1) == 0 && mean(abs(bot_ang)) >= 75
    imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
    im_dend = imbot;    
    imOrig2 = imcrop(imOrig, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
elseif size(bot_prps,1) == 0 && mean(abs(top_ang)) >= 75
    imtop = imrotate(imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
    im_dend = imtop;
    imOrig2 = imrotate(imcrop(imOrig,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
        % One side is empty and the other size is > 65%
elseif botSizeProp > 0.65 && ((size(top_prps,1) == 0 || size(bot_prps,1) == 0))
    imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
    im_dend = imbot;
    imOrig2 = imcrop(imOrig, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
elseif topSizeProp > 0.65 && ((size(top_prps,1) == 0 || size(bot_prps,1) == 0))
    imtop = imrotate(imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
    im_dend = imtop;
    imOrig2 = imrotate(imcrop(imOrig,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
        % One side has only one object. the other has more than 1, and the
        % orientation is close to vertical and axis is longer
elseif size(top_prps,1) == 1 && size(top_prps,1)<size(bot_prps,1) && mean(abs(bot_ang)) >= 75 ...
        && mean(topAxis) < mean(botAxis)
    imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
    im_dend = imbot;
    imOrig2 = imcrop(imOrig, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
elseif size(bot_prps,1) == 1 && size(top_prps,1)>size(bot_prps,1) && mean(abs(top_ang)) >= 75 ...
        && mean(botAxis) < mean(topAxis)
    imtop = imrotate(imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
    im_dend = imtop;
    imOrig2 = imrotate(imcrop(imOrig,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
else
    % Determine which image to use by which has the most large objects
    % closest to parallel and longer
    if top_std < bot_std && mean(abs(top_ang)) >= 80 && mean(topAxis) > mean(botAxis)
        imtop = imrotate(imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
        im_dend = imtop;
        imOrig2 = imrotate(imcrop(imOrig,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
    elseif top_std > bot_std && mean(abs(bot_ang)) >= 80 && mean(topAxis) < mean(botAxis)
        imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
        im_dend = imbot;
        imOrig2 = imcrop(imOrig, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
    elseif abs(max(top_stats(:,3))-max(bot_stats(:,3))) > 0.25 && max(top_stats(:,3)) < max(bot_stats(:,3)) && topSizeProp > 0.45 && topAxisProp > 0.4
        imtop = imrotate(imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
        im_dend = imtop; 
        imOrig2 = imrotate(imcrop(imOrig,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
    elseif abs(max(top_stats(:,3))-max(bot_stats(:,3))) > 0.25 && max(top_stats(:,3)) > max(bot_stats(:,3)) && botSizeProp > 0.45 && botAxisProp > 0.4
        imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
        im_dend = imbot; 
        imOrig2 = imcrop(imOrig, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
    elseif abs(mean(abs(top_ang))-90)  < abs(mean(abs(bot_ang)) - 90) && mean(topAxis) > mean(botAxis) && topSizeProp > 0.45
        imtop = imrotate(imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
        im_dend = imtop;
        imOrig2 = imrotate(imcrop(imOrig,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
    elseif abs(mean(abs(top_ang))-90)  > abs(mean(abs(bot_ang)) - 90) && mean(topAxis) < mean(botAxis) && botSizeProp > 0.45
        imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
        im_dend = imbot;
        imOrig2 = imcrop(imOrig, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
    elseif size(botAxis,1) == 1 && abs(max(bot_stats(:,2))) < 45
        imtop = imrotate(imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
        im_dend = imtop;
        imOrig2 = imrotate(imcrop(imOrig,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
    elseif size(topAxis,1) == 1 && abs(max(top_stats(:,2))) < 45
        imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
        im_dend = imbot; 
        imOrig2 = imcrop(imOrig, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
    elseif topSizeProp > .7
        imtop = imrotate(imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
        im_dend = imtop;
        imOrig2 = imrotate(imcrop(imOrig,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
    elseif botSizeProp > .7 
        imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
        im_dend = imbot;
        imOrig2 = imcrop(imOrig, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
    elseif abs(max(top_stats(:,3))-max(bot_stats(:,3))) > 0.25 && max(top_stats(:,3)) < max(bot_stats(:,3)) && topSizeProp > 0.45
        imtop = imrotate(imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
        im_dend = imtop; 
        imOrig2 = imrotate(imcrop(imOrig,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
    elseif abs(max(top_stats(:,3))-max(bot_stats(:,3))) > 0.25 && max(top_stats(:,3)) > max(bot_stats(:,3)) && botSizeProp > 0.45
        imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
        im_dend = imbot;
        imOrig2 = imcrop(imOrig, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
    elseif median(topAxis) > median(botAxis) 
        imtop = imrotate(imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
        im_dend = imtop;
        imOrig2 = imrotate(imcrop(imOrig,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
    elseif median(topAxis) < median(botAxis) 
        imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
        im_dend = imbot; 
        imOrig2 = imcrop(imOrig, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
        % All else failed - have user choose which side has the dendrites
    else        
        imshowpair(top_bw,bot_bw,'montage');
        imgprompt = ('Can not pick out which side has dendrites! \nPlease pick which side has dendrites: Left (1) or Right (2)');
        user_img = input(imgprompt);
        user_img = str2double(user_img);
        while user_img ~= 1 || user_img ~= 2
            if user_img ~= 1 || user_img ~= 2
                break
            else
                user_img = input(imgprompt);
                user_img = str2double(user_img); 
            end
        end
        if user_img == 1
            imtop = imrotate(imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
            im_dend = imtop; 
            imOrig2 = imrotate(imcrop(imOrig,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
        elseif user_img == 2
            imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
            im_dend = imbot;  
            imOrig2 = imcrop(imOrig, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
        end
        
    end
end 

% Double check with user that correct image was chosen
close all
if side == 1
    usercheck = 'y'; % first go
elseif side == 0
    usercheck = 'n'; % Tracking incorrect so must be wrong
end
if isempty(usercheck)
    usercheck = 'y';
end
while usercheck ~= 'y' || usercheck ~= 'n'
    if usercheck ~= 'y' || usercheck ~= 'n'
        break
    else
    userprompt = ('Is this the correct picture? (Y/N) \n');
    usercheck = input(userprompt);
    usercheck = lower(usercheck);
    end
end

% Ensure that both top and bottom images are not the same, then use the
% other image
topcheck = mean2(imtop);
botcheck = mean2(imbot);
dendcheck = mean2(im_dend);
if usercheck == 'y'
    imOrig = imOrig2;
    close all
elseif usercheck == 'n' && sum(topcheck-dendcheck) == 0
    close all
    imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
    im_dend = imbot;       
    imOrig = imcrop(imOrig, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
elseif usercheck == 'n' && sum(botcheck-dendcheck) == 0
    close all
    imtop = imrotate(imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
    im_dend = imtop;
    imOrig = imrotate(imcrop(imOrig,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
end


%% Step 5: Crop denrite image down to only dendrites - no cell body - allows for dendrite tracking

bordersToKeep = [1 0 0 0]; % Allows for imclearborder but keeping the 'top' border

% Binarize and remove any object that is not touching the top of the image
close all
imbw4 = imbinarize(im_dend,0.2);

% Clearing all but top border
imbw4edgepad = padarray(imbw4,[bordersToKeep(1), bordersToKeep(4)],'pre');
imbw4edgepad = padarray(imbw4edgepad,[bordersToKeep(3), bordersToKeep(2)],'post');
imbw4edgepad2 = imclearborder(imbw4edgepad);
topLeft = [bordersToKeep(1), bordersToKeep(4)]+1;
botRight = topLeft + [size(imbw4,1), size(imbw4,2)] - 1;
imbw4clear =imbw4edgepad2(topLeft(1):botRight(1), topLeft(2):botRight(2)); 
imbw4clear = bwareaopen(imbw4clear,5*round(sFactor)); 
% imshow(imbw4clear)

%% COMMENTED OUT FOR ACCURACY TEST
% Find column values  where image disappears due to roation
clearmax = zeros(size(imbw4clear,2),1);
for ii = 1:size(imbw4clear,2)
   clearmax(ii) = max(imbw4clear(:,ii));
end
clearidx(1) = find(clearmax,1);
clearidx(2) = find(clearmax,1,'last');

% Further processing by eroding objects - attempts to find remaining cell body
imbw5 = imerode(imbw4,strel('rectangle',[2*round(sFactor),10*round(sFactor)])); %% SCALED
imbw5clear = imclearborder(imbw5);
%imshowpair(imbw5,imbw5clear,'montage')

% Remove everything that did not touch an edge of the image
imbw5edge = imbw5;
imbw5edge(imbw5clear) = 0;
imbw5edge = imerode(imbw5edge,strel('line',20*round(sFactor),0)); % Erode it to help obtain bounding box

% Need to clear border everywhere that is not top of the image (in case any
% objects touch the other borders)
imedgepad = padarray(imbw5edge,[bordersToKeep(1), bordersToKeep(4)],'pre');
imedgepad = padarray(imedgepad,[bordersToKeep(3), bordersToKeep(2)],'post');
imedgepad2 = imclearborder(imedgepad);
topLeft = [bordersToKeep(1), bordersToKeep(4)]+1;
botRight = topLeft + [size(imbw5edge,1), size(imbw5edge,2)] - 1;
imbw5edge =imedgepad2(topLeft(1):botRight(1), topLeft(2):botRight(2));
%imshow(imbw5edge)

% Obtain regionprops data from the dendrites (imbw5edge)
close all
imbw5props = regionprops(imbw5edge,'Area','BoundingBox');
rot5_cell = struct2cell(imbw5props);
rot5_area(:,1) = rot5_cell(1,:);
[~, max_idx5] = max(cell2mat(rot5_area(:,1)));
cellbody_bb = cell2mat(rot5_cell(2,max_idx5)); % Obtain bounding box of the largest dendrite

% no need to get rid of cell body (it's not there)
if isempty(cellbody_bb) == 1
    imbw6 = imbw4;
else 
    % Crop images (both binary and grayscale) to only dendrites. (Use imbw4 for height, imbw5 for width)
    imbw6 = imcrop(imbw4,[clearidx(1)-10, cellbody_bb(2)+cellbody_bb(4)+15, clearidx(2)-clearidx(1)+20, size(imbw4,1)]); 
    im_dend = imcrop(im_dend,[clearidx(1)-10, cellbody_bb(2)+cellbody_bb(4)+15, clearidx(2)-clearidx(1)+20, size(imbw4,1)]); 
    imOrig = imcrop(imOrig,[clearidx(1)-10, cellbody_bb(2)+cellbody_bb(4)+15, clearidx(2)-clearidx(1)+20, size(imbw4,1)]);
end

% Dilate and remove any remaining objects that touch any border besides the top
im_dilpad = padarray(imbw6,[bordersToKeep(1), bordersToKeep(4)],'pre');
im_dilpad = padarray(im_dilpad,[bordersToKeep(3), bordersToKeep(2)],'post');
im_dil2 = imclearborder(im_dilpad);
%imshowpair(imbw6,im_dilpad,'montage')

topLeft = [bordersToKeep(1), bordersToKeep(4)]+1;
botRight = topLeft + [size(imbw6,1), size(imbw6,2)] - 1;
im_dil =im_dil2(topLeft(1):botRight(1), topLeft(2):botRight(2));


% Close binary image to create one big blob and get final bounding box
im_close = imclose(im_dil,strel('disk',100,8)); 

% Decrease bounding box to correctect size and finally fix orientation to make
% dendrites as vertical as possible.
close all
close_or = regionprops(bwareaopen(im_close,5*round(sFactor)),'orientation'); % Scaled
close_area = regionprops(bwareaopen(im_close,5*round(sFactor)),'area'); % Scaled
[~, close_idx] = max(cell2mat(struct2cell(close_area)));

% Fix orientation
close_o = cell2mat(struct2cell(close_or(close_idx)));
close_sign = sign(close_o);
if abs(90-abs(close_o)) > 5
    im_close = imrotate(im_close,close_sign*(90-abs(close_o)));
    im_dend = imrotate(im_dend,close_sign*(90-abs(close_o)));
    close_area = regionprops(bwareaopen(im_close,5*round(sFactor)),'area'); % SCALED
    [~, close_idx] = max(cell2mat(struct2cell(close_area)));
    imOrig = imrotate(imOrig,close_sign*(90-abs(close_o)));
end
% imshowpair(im_dend,im_close)


% Fix cropping of greyscale image
close_props = regionprops(bwareaopen(im_close,5*round(sFactor)),'boundingbox'); % SCALED
close_bb = cell2mat(struct2cell(close_props(close_idx)));
close_bb(1) = close_bb(1)-5*round(sFactor*.5); % add a little padding in x
close_bb(3) = close_bb(3)+10*round(sFactor*.5); 
close_bb(4) = close_bb(4)+10*round(sFactor*.5);
dend_crop = imcrop(im_dend,close_bb);
imOrig = imcrop(imOrig,close_bb);
im_dend = dend_crop;
im_dend(im_dend==0) = 1;

end