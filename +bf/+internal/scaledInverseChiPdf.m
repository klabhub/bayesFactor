function y = scaledInverseChiPdf(x,df,scale)
% The Scaled Inverse Chi-Squared Distribution.
% INPUT
% x = the parameter value
% df = degrees of freedom
% scale = scale (tau squared).
% OUTPUT
% y = The probaility density [nrX nrScales]
%
% BK - 2018
if nargin <3
    scale =1; % Default to scaled inverse Chi-squared.
end

[nrX,dimX] = max(size(x));
nrScalesInX = min(size(x));

if iscell(scale)
    scale = [scale{:}];    
    if dimX==1
        scale = reshape(scale,1,[]);        
    else
        scale = reshape(scale,[],1);
    end
end
nrScale = numel(scale);

if nrScale >nrScalesInX
    if dimX==1
        x = repmat(x,[1 nrScale]);
        scale = repmat(scale,[nrX 1]);
    else
        x = repmat(x,[nrScale 1]);
        scale = repmat(scale,[1 nrX]);
    end
end

z = x<0;
x(z) =NaN;

y = bf.internal.inverseGammaPdf(x,df/2,df.*scale/2);
y(z) = 0;
end
