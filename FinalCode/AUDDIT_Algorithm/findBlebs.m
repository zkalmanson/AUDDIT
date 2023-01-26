function [imLocCont,imFeatures,imDend,imOrigBig] = findBlebs(imOrig,sFactor)   
    
    close all

    % Get background standard deviation
    t = adaptthresh(imOrig,.2); 
    imbw = imbinarize(imOrig,t); 
    imNot = ~imbw;
    bkStd = std2(imOrig(imNot));

    % Local contrast adjust cropped image. 
    % The adjustment is dependent on background standard
    % deviation of the image and/or mean of the image intensities (Darker
    % images have a stronger contrast adjustment).
    if bkStd < 0.75
        cropped2 = localcontrast(imOrig,.5,.75);
    elseif bkStd < 1
        cropped2 = localcontrast(imOrig,.5,.75);
    elseif bkStd < 2
        cropped2 = localcontrast(imOrig,1,.5);
    elseif mean2(nonzeros(imOrig)) < 2
        cropped2 = localcontrast(imOrig,0.975,1);
    elseif mean2(nonzeros(imOrig)) < 3
        cropped2 = localcontrast(imOrig,0.95,1); 
    elseif mean2(nonzeros(imOrig)) < 3.5
        cropped2 = localcontrast(imOrig,0.8,1);
    elseif mean2(nonzeros(imOrig)) < 4
        cropped2 = localcontrast(imOrig,0.75,1);
    elseif mean2(nonzeros(imOrig)) < 5 && mean2(nonzeros(imOrig))>= 4
        cropped2 = localcontrast(imOrig,0.5,1);
    elseif mean2(nonzeros(imOrig)) < 7 && mean2(nonzeros(imOrig))>= 5
        cropped2 = localcontrast(imOrig,0.4,1);
    elseif mean2(nonzeros(imOrig)) < 10 && mean2(nonzeros(imOrig))>= 7
        cropped2 = localcontrast(imOrig,0.3,1);       
    elseif mean2(nonzeros(imOrig)) < 13 && mean2(nonzeros(imOrig))>= 10
        cropped2 = localcontrast(imOrig,0.275,1);    
    elseif mean2(nonzeros(imOrig)) < 19 && mean2(nonzeros(imOrig))>= 13
        cropped2 = localcontrast(imOrig,0.25,1);
    elseif mean2(nonzeros(imOrig)) > 19 && std2(nonzeros(imOrig)) < 20 
        cropped2 = localcontrast(imOrig,0.05,.75);
    else % Bright enough with high enough contrast to not use contrast adjustment (usually does not happen)
        cropped2 = imOrig;
    end
    
    % Resize the image by 4x to obtain more pixels to play with
    img = imresize(cropped2,4);
    img(img==1) = NaN;
    imOrigBig = imresize(imOrig,4);
    
    % Filter image for uneven illumination through background subtraction
    se = strel('rectangle',[round(size(imOrig,1)/7)*round(sFactor*2) round(size(imOrig,1)/7)*round(sFactor*2)]); 
    background = imopen(img,se);
    imLocCont = img-background;

    % Global thresholding and binarization
    T = graythresh(imLocCont);

    % Extent of thresholding is dependent on original image background
    % intensity and standard deviation (similar to above)

    % For really dim or really bright images
    if mean2(imLocCont) < 15 && mean(nonzeros(imOrig)) < 6  && bkStd > 2||  mean2(imLocCont) < 15 && mean(nonzeros(imOrig)) > 20
        if mean2(nonzeros(imLocCont)) < 10
            imbwG = imbinarize(imLocCont,T*1.4); 
        elseif mean2(nonzeros(imLocCont)) < 20
            imbwG = imbinarize(imLocCont,1.6*T);
        elseif mean2(nonzeros(imLocCont)) < 30
            imbwG = imbinarize(imLocCont,1.8*T); 
        elseif mean2(nonzeros(imLocCont)) < 40
            imbwG = imbinarize(imLocCont,T*1.9);
        elseif mean2(nonzeros(imLocCont)) < 50
            imbwG = imbinarize(imLocCont,T*1.95);
        elseif mean2(nonzeros(imLocCont)) > 50
            imbwG = imbinarize(imLocCont,T*2); 
        end
    elseif bkStd < 2 
        % For low background noise images
        % Most confocal images
        if mean2(nonzeros(imLocCont)) < 10
            imbwG = imbinarize(imLocCont,T*.9); 
        elseif mean2(nonzeros(imLocCont)) < 15
            imbwG = imbinarize(imLocCont,1.1*T);
        elseif mean2(nonzeros(imLocCont)) < 20
            imbwG = imbinarize(imLocCont,1.15*T);
        elseif mean2(nonzeros(imLocCont)) < 30
            imbwG = imbinarize(imLocCont,1.2*T); 
        elseif mean2(nonzeros(imLocCont)) < 40
            imbwG = imbinarize(imLocCont,T*1.25);
        elseif mean2(nonzeros(imLocCont)) < 50
            imbwG = imbinarize(imLocCont,T*1.3);
        elseif mean2(nonzeros(imLocCont)) > 50
            imbwG = imbinarize(imLocCont,T*1.5); 
        end
    else
        % Most other images
        if mean2(nonzeros(imLocCont)) < 10
            imbwG = imbinarize(imLocCont,T*.9); 
        elseif mean2(nonzeros(imLocCont)) < 15
            imbwG = imbinarize(imLocCont,1.05*T);
        elseif mean2(nonzeros(imLocCont)) < 20
            imbwG = imbinarize(imLocCont,1.1*T);
        elseif mean2(nonzeros(imLocCont)) < 30
            imbwG = imbinarize(imLocCont,1.2*T); 
        elseif mean2(nonzeros(imLocCont)) < 40
            imbwG = imbinarize(imLocCont,T*1.3);
        elseif mean2(nonzeros(imLocCont)) < 50
            imbwG = imbinarize(imLocCont,T*1.4);
        elseif mean2(nonzeros(imLocCont)) > 50
            imbwG = imbinarize(imLocCont,T*1.5); 
        end
    end
    % figure(); imshow(imbwG)


    % Creates a binary mask, which is dependent on contrast adjusted image,
    % to help determine if dendrites are visible (attempts to roughly
    % segment dendrites).
    if mean2(imLocCont) < 10
        if bkStd > 2
            imDend = imbinarize(imLocCont,(.9*T)); 
        else
            imDend = imbinarize(imLocCont,(.5*T));
        end

        % check if too much binarized
        iProp = length(find(imDend))/(size(imDend,1)*size(imDend,2));
        z=1.25;

        % While loop that changes threshold to ensure not too much background is segmented 
        while iProp > 0.225
            z=z+.25;
            imDend = imbinarize(imLocCont,(z*T));
            iProp = length(find(imDend))/(size(imDend,1)*size(imDend,2)); 
        end

    elseif mean2(imLocCont) < 15
        imDend = imbinarize(imLocCont,(.5*T)); 

        % check if too much binarized
        iProp = length(find(imDend))/(size(imDend,1)*size(imDend,2));
        z=.75;
        while iProp > 0.225
            z=z+.25;
            imDend = imbinarize(imLocCont,(z*T));
            iProp = length(find(imDend))/(size(imDend,1)*size(imDend,2)); 
        end

   elseif mean2(imLocCont) < 20
        imDend = imbinarize(imLocCont,(.75*T)); 

        % check if too much binarized
        iProp = length(find(imDend))/(size(imDend,1)*size(imDend,2));
        z=.9;
        while iProp > 0.225
            z=z+.25;
            imDend = imbinarize(imLocCont,(z*T));
            iProp = length(find(imDend))/(size(imDend,1)*size(imDend,2)); 
        end   

    % Bright images below have shown to not need the check
    elseif mean2(imLocCont) < 25
        imDend = imbinarize(imLocCont,(1.2*T)); 
    elseif mean2(imLocCont) < 30
        imDend = imbinarize(imLocCont,(1.4*T)); 
    elseif mean2(imLocCont) < 35
        imDend = imbinarize(imLocCont,(1.6*T)); 
    elseif mean2(imLocCont) > 35
        imDend = imbinarize(imLocCont,(2*T));        
    end
    imDend = bwareaopen(imDend,3*round(sFactor));
    
    % Creates a binary mask to get image intensity near the dendrites
    % Binarization adapts for dim images with little to no dendrite by 
    % finding if mean pixels is under 5.
    if mean2(nonzeros(imLocCont)) < 15
        imbwG2 = imbinarize(imLocCont,(T)); 
    elseif mean2(nonzeros(imLocCont)) > 50
        imbwG2 = imbinarize(imLocCont,(.55*T)); 
    else
        imbwG2 = imbinarize(imLocCont,(.85*T)); 
    end
    
    % Create mask for dendrites, and surrounding pixels
    imMask = bsxfun(@times, imLocCont, cast(imbwG2, 'like', imLocCont));
    imMaskInt = mean2(imMask);

    % Erode global thresholded segemntation. Erosion is dependent on
    % background pixel noise and mean pixel intensities of the dendrites
    % and surrounding pixels (obtained with imbwG2 and imMask). Brighter
    % images required larger erosion structuring elements.
    if bkStd > 2
        if imMaskInt < 6 
            imE = bwareaopen(imerode(imbwG,strel('rectangle',[ceil(sFactor),10*ceil(sFactor*(1/2))])),20*round(sFactor/2));
        elseif imMaskInt < 8 && imMaskInt >= 6
            imE = bwareaopen(imerode(imbwG,strel('rectangle',[ceil(sFactor),14*ceil(sFactor*(1/2))])),20*round(sFactor/2));
        elseif imMaskInt < 10 && imMaskInt >= 8
            imE = bwareaopen(imerode(imbwG,strel('rectangle',[2*ceil(sFactor),16*ceil(sFactor*(1/2))])),20*round(sFactor/2));      
        elseif imMaskInt < 14 && imMaskInt >= 10
            imE = bwareaopen(imerode(imbwG,strel('rectangle',[2*ceil(sFactor),16*ceil(sFactor*(1/2))])),20*round(sFactor/2)); 
        elseif imMaskInt < 18 && imMaskInt >= 14
            imE = bwareaopen(imerode(imbwG,strel('rectangle',[2*ceil(sFactor),18*ceil(sFactor*(1/2))])),20*round(sFactor/2));      
        elseif imMaskInt < 20 && imMaskInt >= 18
            imE = bwareaopen(imerode(imbwG,strel('rectangle',[3*ceil(sFactor),18*ceil(sFactor*(1/2))])),20*round(sFactor/2)); 
        else
            imE = bwareaopen(imerode(imbwG,strel('rectangle',[4*ceil(sFactor),20*ceil(sFactor*(1/2))])),20*round(sFactor/2));
        end
    else
        if imMaskInt < 2 
            imE = bwareaopen(imerode(imbwG,strel('rectangle',[ceil(sFactor),10*ceil(sFactor*(1/2))])),20*round(sFactor/2));
        elseif imMaskInt < 4 
            imE = bwareaopen(imerode(imbwG,strel('rectangle',[ceil(sFactor),13*ceil(sFactor*(1/2))])),20*round(sFactor/2));
        elseif imMaskInt < 6
            imE = bwareaopen(imerode(imbwG,strel('rectangle',[ceil(sFactor),15*ceil(sFactor*(1/2))])),20*round(sFactor/2));
        elseif imMaskInt < 8 && imMaskInt >= 6
            imE = bwareaopen(imerode(imbwG,strel('rectangle',[ceil(sFactor),16*ceil(sFactor*(1/2))])),20*round(sFactor/2));
        elseif imMaskInt < 10 && imMaskInt >= 8
            imE = bwareaopen(imerode(imbwG,strel('rectangle',[ceil(1*sFactor),16*ceil(sFactor*(1/2))])),20*round(sFactor/2));      
        elseif imMaskInt < 14 && imMaskInt >= 10
            imE = bwareaopen(imerode(imbwG,strel('rectangle',[ceil(1*sFactor),16*ceil(sFactor*(1/2))])),20*round(sFactor/2)); 
        elseif imMaskInt < 18 && imMaskInt >= 14
            imE = bwareaopen(imerode(imbwG,strel('rectangle',[ceil(1*sFactor),16*ceil(sFactor*(1/2))])),20*round(sFactor/2));      
        elseif imMaskInt < 20 && imMaskInt >= 18
            imE = bwareaopen(imerode(imbwG,strel('rectangle',[ceil(1*sFactor),18*ceil(sFactor*(1/2))])),20*round(sFactor/2)); 
        else
            imE = bwareaopen(imerode(imbwG,strel('rectangle',[ceil(6*sFactor),36*ceil(sFactor*(1/2))])),20*round(sFactor/2));
        end        
    end


    % Create one more mask to try to make sure every feature segmented,
    % should be segmented.
    dendPixels = imLocCont(imLocCont > T.*max(imLocCont,[],'all'));
    T2 = graythresh(dendPixels);
    if mean2(nonzeros(imLocCont)) > 40
        imbw = bwareaopen(imbinarize(imLocCont,T2*1.2),3*round(sFactor));
    else
        imbw = bwareaopen(imbinarize(imLocCont,T2*.9),1*round(sFactor));
    end

    % Remove any dark pixels from feature extraction
    imE2 = imE;
    imE2(imbw==0) = 0; % remove pixels that are 0 on the max intensity

    imE2 = bwareaopen(imE2,20*round(sFactor)); %Get rid of small detected features
    imE2 = imclearborder(imE2); % Lower chances of detecting cell body as a feature

    % Active contour segementation to get improved feature shape
    % definition. Changing contraction bias and smooth factor values
    % changes definition of features.
    im2 = bwmorph(imE2,'thicken',4*ceil(sFactor/2)); 
    imFeatures = bwareaopen(activecontour(imLocCont,im2,50*round(sFactor/2),'Chan-Vese','SmoothFactor',2,'ContractionBias',.25),round(sFactor*5)); 
    
%   figure(); imshowpair(imActive,img2)

    
end