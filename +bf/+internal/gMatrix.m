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

