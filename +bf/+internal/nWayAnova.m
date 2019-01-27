function bf10 = nWayAnova(y,X,varargin)
% Bayes Factor analysis for an N-Way Anova.
% Don't call this directly, use bf.anova instead.
% y = data values
% X = design matrix for ANOVA (indicator vars)  no constant term
%
% Parm/Value pairs
% 'sharedPriors'  - Cell array of vectors indicating which effects (columns
% of X) share the same prior. [{1:nrEffects}]: all effects share the same prior.
% 'options' - Monte Carlo integration and parrallel computation optionss. [bf.options]
% BK 2018
nrEffects = size(X,2);

p =inputParser;
p.addParameter('sharedPriors',{},@iscell); % Which effects share a prior? A cell array with indices corresponding to columns of X
p.addParameter('options',bf.options);
p.parse(varargin{:});

if isempty(p.Results.sharedPriors)
    sharedPriors = {1:nrEffects};
else
    sharedPriors = p.Results.sharedPriors;
end

prior = @(g)(bf.internal.scaledInverseChiPdf(g,1,1));
integrand = @(varargin) (bf.internal.rouderS(cat(1,varargin{:}),y,X,sharedPriors,p.Results.options).*prod(prior(cat(1,varargin{:})),1));
nrDims = numel(sharedPriors);
if nrDims>= p.Results.options.nDimsForMC
    % Use MC Sampling to calculate the integral
    bf10 = bf.internal.mcIntegral(integrand,prior,nrDims,p.Results.options);
else
    switch (nrDims)
        case 1
            bf10 = integral(integrand,0,Inf);
        case 2
            bf10 = integral2(integrand,0,Inf,0,Inf);
        case 3
            bf10 = integral3(integrand,0,Inf,0,Inf,0,Inf);
    end
end
end
