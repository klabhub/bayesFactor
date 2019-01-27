function y = inverseGammaPdf(x,alpha,beta)
% The inverse Gamma PDF.
% INPUT
% x  (>0)
% alpha - shape parameter
% beta  - scale parameter
%
% BK - 2018
%assert(all(x>0),'The inverse gamma PDF is only defined for x>0')
z = x<0;
y = zeros(size(z));
y(~z) = (beta.^alpha)/gamma(alpha)*(1./x(~z)).^(alpha+1).*exp(-beta./x(~z));
end
