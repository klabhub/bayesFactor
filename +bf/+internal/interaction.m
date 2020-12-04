function X = interaction(varargin)
% Create all interaction terms from two or more dummy coded design matrices.
% (See Box II in Rouder et al. 2012)
%
% INPUT
% Xa, Xb, Xc,... = design matrices with matching number of observation (rows)
%
% OUTPUT
% X = Dummy coded design matrix [nrObservations nrA*nrB*nrC...]
%
% BK - 2019

nObs = cellfun(@(x) size(x,1),varargin);
if ~numel(unique(nObs))==1
    error('Design matricces should all have the same number of observations (%d)',nObs);
end
X = varargin{1};
for i=2:nargin
    nB = size(varargin{i},2);
    nA = size(X,2);
    Xa = repmat(X,[1 nB]);
    Xb = repmat(varargin{i},[1 nA]);
    ix = (repmat(1:nA:nA*nB,[nA 1]) + repmat((0:nA-1)',[1 nB]))';
    Xa = Xa(:,ix(:));
    X = Xa.*Xb;
end

end