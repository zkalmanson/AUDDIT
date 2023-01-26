% Function used to determine where breaks are located
% Updated: 08-30-2022

function [dendBreak,breakImg1,breakImg2,breakImg3,breakImg4,break1start,break2start,break3start,break4start,dendIntensity]...
    = breakAnalysis(n1,n2,n3,n4,imDend,imOrigBig)

    % Initialization
    dend1 = n1 & imDend;
    dend2 = n2 & imDend;
    dend3 = n3 & imDend;
    dend4 = n4 & imDend;    
    break1 = zeros(1,size(n1,1));
    break2 = zeros(1,size(n1,1));
    break3 = zeros(1,size(n1,1));
    break4 = zeros(1,size(n1,1));
    breakImg1 = zeros(size(n1));
    breakImg2 = zeros(size(n1));
    breakImg3 = zeros(size(n1));
    breakImg4 = zeros(size(n1));
    n1 = logical(n1);
    n2 = logical(n2);
    n3 = logical(n3);
    n4 = logical(n4);

    % Finds rows where no dendrite is detected
    for ii = 1:size(n1,1)
        if isempty(find(dend1(ii,:),1)) == 1 && isempty(find(imOrigBig(ii,n1(ii,:)), 1)) ~= 1 % don't count breaks where the background is completely black 
            break1(ii) = 1;
            breakImg1(ii,n1(ii,:)) = 1;
        end
        if isempty(find(dend2(ii,:),1)) == 1 && isempty(find(imOrigBig(ii,n2(ii,:)), 1)) ~= 1
            break2(ii) = 1;
            breakImg2(ii,n2(ii,:)) = 1;
        end
        if isempty(find(dend3(ii,:),1)) == 1 && isempty(find(imOrigBig(ii,n3(ii,:)), 1)) ~= 1
            break3(ii) = 1;
            breakImg3(ii,n3(ii,:)) = 1;
        end
        if isempty(find(dend4(ii,:),1)) == 1 && isempty(find(imOrigBig(ii,n4(ii,:)), 1)) ~= 1
            break4(ii) = 1;
            breakImg4(ii,n4(ii,:)) = 1;
        end
    end

    % Combine variables and create break images
    dendBreak = [break1; break2; break3; break4];
    breakImg1 = logical(breakImg1);
    breakImg2 = logical(breakImg2);
    breakImg3 = logical(breakImg3);
    breakImg4 = logical(breakImg4);

    % Determine which pixel row the dendrite starts - avoids issues of
    % black cropped background from rotating image earlier.
    % Assume it starts at pixel row 1
    break1start = 1;
    break2start = 1;
    break3start = 1;
    break4start = 1;
    
    % Change the start but do not go further than a 1/3rd of the image
    for ii = 1:(round(size(n1,1)/3))
        if isempty(find(dend1(ii,:),1)) == 1 && isempty(find(imOrigBig(ii,n1(ii,:)), 1)) == 1 
            break1start = ii;
        end
        if isempty(find(dend2(ii,:),1)) == 1 && isempty(find(imOrigBig(ii,n2(ii,:)), 1)) == 1
            break2start = ii;
        end
        if isempty(find(dend3(ii,:),1)) == 1 && isempty(find(imOrigBig(ii,n3(ii,:)), 1)) == 1
            break3start = ii;
        end
        if isempty(find(dend4(ii,:),1)) == 1 && isempty(find(imOrigBig(ii,n4(ii,:)), 1)) == 1
            break4start = ii;
        end
    end
    
    % Get intensities of dendrites using original images
    dend1img = bsxfun(@times, imOrigBig, cast(dend1, 'like', imOrigBig));
    dend2img = bsxfun(@times, imOrigBig, cast(dend2, 'like', imOrigBig));
    dend3img = bsxfun(@times, imOrigBig, cast(dend3, 'like', imOrigBig));
    dend4img = bsxfun(@times, imOrigBig, cast(dend4, 'like', imOrigBig));
    
    % Mean intensities of the dendrites
    dend1intensity = mean2(nonzeros(dend1img));
    dend2intensity = mean2(nonzeros(dend2img));
    dend3intensity = mean2(nonzeros(dend3img));
    dend4intensity = mean2(nonzeros(dend4img));
    
    % Combine into one variable
    dendIntensity = [dend1intensity;dend2intensity;dend3intensity;dend4intensity];


end