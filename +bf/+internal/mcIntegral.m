function v = mcIntegral(fun,prior,options)
% Monte Carlo integration
%
% INPUT
% fun -  The function to integrate. This should be specified as a
%       function_handle that takes a single input (g)
% prior - the prior distribution of the g's. A function_handle.
%           used to do importance sampling
% nrDims - The number of dimensions to integrate over. [1]
% 
% options - A struct with options.  []
% OUTPUT
% v -  The value of the integral. (Typically the BF10).


%% Setup the PDF to do importance sampling
gRange =  (options.minG:options.stepG:options.maxG);
pdf = prior(gRange);
pdf = pdf./sum(pdf,2);% Normalize to use as weight
% Draw samples weighted by this prior.
nrDims = size(pdf,1);
g =nan(nrDims,options.nrSamples);
for d=1:nrDims    
    g(d,:) = randsample(gRange,options.nrSamples,true,pdf(d,:));
end
%% Evaluate the function at these g values
bf10Samples = fun(g);
v = mean(bf10Samples); % Expectation value ~ integral.
end
