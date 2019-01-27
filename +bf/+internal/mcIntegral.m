function v = mcIntegral(fun,prior,nrDims,options)
% Monte Carlo integration
%
% INPUT
% fun -  The function to integrate. This should be specified as a
%       function_handle that takes a single input (g)
% prior - the prior distribution of the g's. A function_handle.
%
% nrDims - The number of dimensions to integrate over. [1]
% options - A struct with options.  []
% OUTPUT
% v -  The value of the integral. (Typically the BF10).


%% Setup the PDF to do importance sampling
gRange =  options.minG:options.stepG:options.maxG;
pdf = prior(gRange);
pdf = pdf./sum(pdf);
% Draw samples weighted by this prior.
g =nan(nrDims,options.nrSamples);
for d=1:nrDims
    g(d,:) = randsample(gRange,options.nrSamples,true,pdf);
end
%% Evaluate the function at these g values
bf10Samples = fun(g);
pg = prod(prior(g),1);  % Probability of each g combination
v = mean(bf10Samples./pg); % Expectation value- = integral.
end
