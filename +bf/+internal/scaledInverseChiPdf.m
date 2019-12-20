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

nrX = size(x,1);
if iscell(scale)
    scale = [scale{:}];
end
nrScale = numel(scale);
if isvector(scale) && nrScale>1
    scale = scale(:)'; %Force Row
    scale =repmat(scale,[nrX 1 ]);
end
if isvector(x) && nrScale >1
    x = x(:); %Force col
    x = repmat(x,[1 nrScale]);
end
z = x<0;
x(z) =NaN;

y = bf.internal.inverseGammaPdf(x,df/2,df.*scale/2);
y(z) = 0;
end
