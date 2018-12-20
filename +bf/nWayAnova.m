function bf10 = nWayAnova(y,X,varargin)
% ANOVA BF
% y = data values
% X = design matrix for ANOVA (indicator vars)  no constant term
%
% Parm/Value pairs
% 'sharedPriors'  - Cell array of vectors indicating which effects (columns
% of X) share the same prior. [{1:nrEffects}]: all effects share the same prior.
%
% BK 2018
nrEffects = size(X,2);

p =inputParser;
p.addParameter('sharedPriors',{},@iscell); % Which effects share a prior? A cell array with indices corresponding to columns of X
p.addParameter('mc',[],@isstruct);  % A struct to specify Monte Carlo integration options. See bf.mcIntegral for details
p.addParameter('nDimsForMC',4,@(x) (x<=4)); % By default 1,2,3 dimensional integrals are done with regular integration, and 4 and higher use MC integration.
p.parse(varargin{:});

if isempty(p.Results.sharedPriors)
    sharedPriors = {1:nrEffects};
else
    sharedPriors = p.Results.sharedPriors;
end

prior = @(g)(bf.scaledInverseChiPdf(g,1,1));
integrand = @(varargin) (rouderS(cat(1,varargin{:}),y,X,sharedPriors).*prod(prior(cat(1,varargin{:})),1));
nrDims = numel(sharedPriors);
if nrDims>= p.Results.nDimsForMC
    % Use MC Sampling to calculate the integral
    bf10 = bf.mcIntegral(integrand,prior,nrDims,p.Results.mc);
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


function G = gMatrix(grouping,g)
% Generate a matrix in which each row corresponds to an effect, each column a value
% that will be integrated over. The grouping cell array determines which of
% the effects share a prior (i.e. levles of the same factor) and which have
% their own.
% grouping   - Cell array with vectors that contain effect indices (i.e.
% columns of the design%matrix) that share a prior.
% g  - The values for each of the priors. Each row is an independent prior,
%               each column is a sample
%
% EXAMPLE
% gMatrix({1 2],[3 4]},[0 0.1 0.2 0.3; 0.6 0.7 0.8 0.9])
% g = [ 0 0.1 0.2 0.3;
%       0 0.1 0.2 0.3;
%       0.6 0.7 0.8 0.9;
%       0.6 0.7 0.8 0.9]


nrEffects = sum(cellfun(@numel,grouping));
assert(nrEffects>0,'The number of groups must be at least one');
nrValues = size(g,2);
G = nan(nrEffects,nrValues);
for i=1:numel(grouping)
    G(grouping{i},:) = repmat(g(i,:),[numel(grouping{i}) 1]);
end
end

function value= rouderS(g,y,X,grouping)
% The S(g) function of Eq 9 in Rouder et al.
% g = Matrix of g values, each row is an effect, each column is a value
% that we're integrating over.
% y = Data values
% X = design matrix, without a constant term, with indicator variables only

g = gMatrix(grouping,g);

nrObservations = size(X,1);
one = ones(nrObservations,1);
P0 = 1./nrObservations*(one*one');
yTilde = (eye(nrObservations)-P0)*y;
XTilde = (eye(nrObservations)-P0)*X;
nrPriorValues=size(g,2);
value= nan(1,nrPriorValues);
for i=1:nrPriorValues
    if all(g(:,i)==0)
        value(i)=0;
    else
        G = diag(g(:,i));
        invG = diag(1./g(:,i));
        Vg = XTilde'*XTilde + invG;
        yBar = one'*y/nrObservations;
        preFactor= 1./(sqrt(det(G))*sqrt(det(Vg)));
        numerator =    y'*y-nrObservations*yBar^2;
        denominator = ((yTilde'*yTilde) -yTilde'*XTilde*(Vg\XTilde'*yTilde));
        value(i)= preFactor*(numerator/denominator).^((nrObservations-1)/2);
    end
end
end