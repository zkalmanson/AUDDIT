% Create a function to smooth data right before interpolation

% Updated: August 30, 2022

function [dendSmooth,dendSmoothIdx] = NeurSmoothData(dendId,dendSkel2,sFactor)
    % Initialize
    dendSmoothIdx = dendId;
    dendSmooth = zeros(size(dendSkel2));

    % Find NaN values or values less than 1 (out of bounds)
    dendSmoothIdx(dendSmoothIdx<=1) = NaN;
    dendnan = isnan(dendSmoothIdx);
    
    % Smooth data using built in function
    dend1smoothidx_temp = smoothdata(dendSmoothIdx,'rloess',35*round(sFactor)*2);
    dendSmoothIdx = dend1smoothidx_temp;

    % Replace NaN or out of bounds values as 0 in smoothed data
    dend1nan = isnan(dendSmoothIdx);
    dendSmoothIdx(dend1nan) = 0;
    dendSmoothIdx(dendSmoothIdx<=1) = 0;
    for zz = 1:length(dendId)
        if dendnan(zz) == 0
            dendSmooth(zz,round(dendSmoothIdx(zz))+1) = 1;
        elseif dendnan(zz) == 1
            dendSmoothIdx(zz) = 0;
        end
    end

    % Create a smoothed image
    dendSmooth = logical(dendSmooth(:,2:end));
end
