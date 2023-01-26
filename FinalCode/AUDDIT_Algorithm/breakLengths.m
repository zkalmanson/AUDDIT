% Function used to determine how long breaks are and the sequences of
% breaks

% Updated: 08-30-2022

function [oneStart,oneLen,k1] = breakLengths(b1)
    % Initialize
    n = length(b1);
    oneStart = zeros(1,n); % Where breaks starts
    oneLen = zeros(1,n); % Length of break
    k1= 0; % Number of breaks
    inOnes = 0;

    % For loop over image size
    for k = 1:n
        if b1(k) == 0 % this row is not in break
            inOnes = 0;
        elseif ~inOnes  % Start of new break
            inOnes = 1;
            k1 = k1+1;
            oneStart(k1) = k;
            oneLen(k1) = 1;
        else     % The break is continuing                   
            oneLen(k1) = oneLen(k1)+1;
        end
    end
    
    % Remove unwanted zeros
    oneLen(oneLen==0) = []; 
    oneStart(oneStart==0) = [];

    if isempty(oneStart)
        oneStart = size(b1,2);
    end

end