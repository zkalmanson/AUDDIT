% Create a function to binarize and detect blebs from individual dendrite
% regions
% Now removes blebs from the dendrite thickness
% Breaks are removed by finding '0' widths

function [width1] = findDendriteWidth(dendimg,dend1img,imFeatures,sFactor)
    
    % Remove detected feature from dendrites
    imRmv = imdilate(imFeatures,strel('line',40*round(sFactor/2),0));
    dendRmv = dend1img;
    dendRmv(imRmv) = 0;

    % Make image into a double for signal processing
    dend1test = double(dendRmv);
    dend1test(dend1test<=1) = NaN;
    dend1mvmean = movmean(mean(dend1test,2,'omitnan'),100,'omitnan');
    changedist = 5; % Total change distance allowed
    minpoint = 50; 
    fallidx = [];

    % Create image thresholds and initialize for speed
    dendthresh = double(prctile(nonzeros(dendimg),25));
    dend1thresh = double(prctile(nonzeros(dend1img),25));
    dendidx = zeros(size(dend1img,1),2);
    dendtest = zeros(size(dend1img));
    dendtest2 = zeros(size(dend1img));
    width1 = zeros(size(dend1img,1),1);

    % Remove points under the threshold
    if dend1thresh < dendthresh
        dend1thresh = dendthresh;
    end
    for ii = 1:size(dend1img,1)
        dendtest(ii,:) = dend1test(ii,:)-dend1thresh;
        dendtest2(ii,:) = dend1test(ii,:)-.25*dend1mvmean(ii);

        % If no change points were found, use the observation range as
        % starting points.
        if ~isempty(find(dend1img(ii,:), 1))
            dendidx(ii,1) = find(dend1img(ii,:),1,'first');
            dendidx(ii,2) = find(dend1img(ii,:),1,'last');
        else
            dendidx(ii,1) = 1;
            dendidx(ii,2) = 1;
        end
    end
    % Remove out of bounds indicies
    dendtest(dendtest<=0) = NaN;
    
    % Find data change points along dendrite
    [chngepts, S1] = ischange(dendtest,2,'maxnumchanges',2);
    dendbw = zeros(size(chngepts));

    % Find the width of the dendrite throughout image using change points
    for ii = 1:size(dend1img,1)
        templine = chngepts(ii,:);
        tempS = S1(ii,:);
        tempS(isnan(tempS)) = 0;
        mintempS(ii) = min(tempS);
        % Determine if the change point is going up our down in intensity
        oneidxDif = diff(dendtest(ii,:));
        if mintempS(ii) >= 5
            mintempS(ii) = 5;
        end
        mintempChngeS = tempS-mintempS(ii);

        % Find how many change points there are in this row
        oneidx = find(templine);
        
        % If 2 change points are detected, make sure they are correct
        if length(oneidx) == 2
            oneidxSign = oneidxDif(oneidx-1);

            % Remove if change points are in the wrong direction
            if sign(oneidxSign(1)) == sign(oneidxSign(2))
                if abs(oneidxSign(1)) > abs(oneidxSign(2))
                    oneidx(2) = [];
                elseif abs(oneidxSign(1)) < abs(oneidxSign(2))
                    oneidx(1) = [];
                end
            end
        end
        
        % If no change points are detected, try to place at area with high
        % slope
        if isempty(oneidx)
            temponeidx = tempS;
            temponeidx(isnan(temponeidx)) = 0;
            tempidx = find(temponeidx,1);
            if isempty(tempidx) == 1
                oneidx = [1 2];
            elseif temponeidx(tempidx) > changedist
                oneidx = [tempidx tempidx];
            else
                oneidx = [1 2];
            end
        end

        % If only one was detected, based on the sign make the beginning or
        % end of the range the other point.
        if length(oneidx) == 1
           oneidxSign = oneidxDif(oneidx-1);
           oneidxS  = tempS(oneidx);
           maxoneidx = max(tempS);
           for jj = 1:length(templine)
               if tempS(jj) == oneidxS && oneidxS == maxoneidx
                  oneidx(2) = jj; 
               elseif tempS(jj) == oneidxS
                   oneidx(2) = oneidx(1);
               end
           end
          if oneidxSign < 0
              if dendtest(ii,dendidx(ii,1)) >= dend1thresh*2
               if oneidx(1) == oneidx(2)
                   tempone = oneidx(1);
                   oneidx(1) = dendidx(ii,1);
                   oneidx(2) = tempone;
               end
              end
          elseif oneidxSign > 0
              if dendtest(ii,dendidx(ii,2)) <= mean(nonzeros(dend1img),'omitnan')
               if oneidx(1) == oneidx(2)
                   tempone = oneidx(1);
                   oneidx(2) = dendidx(ii,2);
                   oneidx(1) = tempone;
               end
              end
          end
        end
        
        % Fix for logic issues
        if oneidx(2) == 1
            oneidx(2) = 2;
        end
        
        tempchnge(1) = abs(tempS(oneidx(1))-mintempS(ii));
        tempchnge(2) = tempS(oneidx(2)-1)-tempS(oneidx(2));
        
        % Try to find location where the change distance is large enough
        while tempchnge(1) < changedist
            nextPt = find(mintempChngeS>=changedist,1,'first');
            if isempty(nextPt) == 1
                if tempS(dendidx(ii,1)) > minpoint
                    oneidx(1) = dendidx(ii,1);
                    tempS(isnan(tempS)) = 0;
                    fallidx = find(tempS(oneidx(1):dendidx(ii,2))<=minpoint,1,'first');
                    oneidx(2) = fallidx;
                    break
                else
                    oneidx(1) = 1;
                    oneidx(2) = 1;
                    break
                end
            else
                oneidx(1) = nextPt;
                tempchnge(1) = abs(tempS(oneidx(1))-mintempS(ii));
            end
        end
        
        if isempty(fallidx) == 1
           fallidx = 0; 
        end

        % Adapt indicies to make sure it make logical sense
        if tempchnge(2) < 0 && oneidx(1) ~= 1
            nextPt = find(tempS(oneidx(2):dendidx(ii,2))<tempS(oneidx(2)),1,'first');
            origPt = oneidx(2);
            if origPt == 1
                origPt = 2;
            end
            if isempty(nextPt) == 1
                nextPt = dendidx(ii,2);
                oneidx(2) = dendidx(ii,2);
            end
            tempchnge(2) = tempS(origPt-1)-tempS(nextPt);
            if tempchnge(2) > changedist
                oneidx(2) = nextPt;
            else
                oneidx(2) = dendidx(ii,2);
            end
        end
        if oneidx(2) == 2
            oneidx(2) = 1;
        end
        if oneidx(2) == 1
            oneidx(2) = oneidx(1);
        end

        % Idicies are set, find the with of this row of the dendrites.
        dendbw(ii,oneidx(1):oneidx(2)) = 1;
        dendbw(ii,1) = 0;
        width1(ii) = length(find(dendbw(ii,:)));

    end
end