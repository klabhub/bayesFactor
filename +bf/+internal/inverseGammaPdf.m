function y = inverseGammaPdf(x,alpha,beta)
% The inverse Gamma PDF.
% INPUT
% x  (>0)
% alpha - shape parameter
% beta  - scale parameter
% OUTPUT
% y = the pdf [nrX nrBeta]
% BK - 2018
if isvector(x);x=x(:);end % Force col
nrX = size(x,1);
nrBeta = numel(beta);
if isvector(beta) && nrBeta>1
    beta = beta(:)'; % Force row
    beta =repmat(beta,[nrX  1]);
end
if isvector(x)
    x = repmat(x,[1 nrBeta]);
end
z = x<0;
x(z) =NaN;
y = (beta.^alpha)./gamma(alpha).*(1./x).^(alpha+1).*exp(-beta./x);
y(z) = 0;
end
