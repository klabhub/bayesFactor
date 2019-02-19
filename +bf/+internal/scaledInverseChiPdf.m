function y = scaledInverseChiPdf(x,df,scale)
% The Scaled Inverse Chi-Squared Distribution.
% INPUT
% x = the parameter value
% df = degrees of freedom
% scale = scale (tau squared).
% OUTPUT
% y = The probaility density
%
% BK - 2018
z = x<0;
y = zeros(size(z));
if nargin <3
    scale =1; % Default to scaled inverse Chi-squared.
end
y(~z) = bf.internal.inverseGammaPdf(x(~z),df/2,df*scale/2);
end
