% Function to get 1 point per y value
% Updated: August 30, 2022

function [dendSkel] = SinglePointBin(dendSkel,dend_lbl,dendCent,dendMax,dendNum,neur_mvmean,neur_idx_avg)

    % Determine how large of deviation is permitted
    pxldev = 7;

    % For each pixel row, keep the point that is closest to the moving mean
    for yy = 1:size(dendSkel,1)
        % Find how many detected points on current pixel row
        dend2_temp = dendSkel(yy,:);
        tempidx = find(dend2_temp);
        tempobj = dend_lbl(yy,tempidx);
        tempcent = dendCent(tempobj);
        
        % If multiple found, determine how far they are from the moving
        % mean
        if length(tempidx) > 1
            tempdist = zeros(length(tempidx),1);
            objtempdist = zeros(length(tempidx),1);
            for zz = 1:length(tempidx)
                tempdist(zz) = abs(tempidx(zz)-neur_mvmean(yy,dendNum));
                objtempdist(zz) = abs(tempcent(zz)-neur_idx_avg(dendNum));
            end
            
            % Keep the point that is closer to the moving mean
            [~, mindist_idx] = min(tempdist);
            real_idx = tempidx(mindist_idx);
            [~, objmindist_idx] = min(objtempdist);
            objreal_idx = tempidx(objmindist_idx);

            w = 0;    
            for i = 1:length(tempidx)
                % Keep the object if it is larger
                if tempobj(i) == dendMax
                    dendSkel(yy,tempidx(i)) = 1;
                    zero_idx = tempidx([1:i-1,i+1:end]);
                    dendSkel(yy,zero_idx) = 0;
                    break
                elseif tempobj(i) ~= dendMax
                    % Remove if it has no relation
                    if tempidx(i) ~= real_idx && tempidx(i) ~= objreal_idx
                        dendSkel(yy,tempidx(i)) = 0;
                    elseif real_idx ~= objreal_idx
                        if tempidx(i) ~= real_idx && tempidx(i) == objreal_idx
                           objidxdist = objtempdist(i);
                           obj_id = i;
                           w = w+1;
                        elseif tempidx(i) == real_idx && tempidx(i) ~= objreal_idx                    
                           idxdist = tempdist(i);
                           id = i;
                           w=w+1;
                        end  
                    end
                    if w == 2
                        if idxdist > objidxdist 
                            dendSkel(yy,tempidx(id)) = 0;
                        elseif idxdist < objidxdist
                            dendSkel(yy,tempidx(obj_id)) = 0;          
                        end
                    end                
                end
            end

        else
            % Remove point if it deviates from moving mean by too much
            % Not commoon
            tempdist = abs(tempidx-neur_mvmean(yy,dendNum));
            if tempdist > pxldev
                dendSkel(yy,tempidx) = 0;
            end
        end
    end
end