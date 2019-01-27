function out = simulateLinearModel(o,lm,effectSize,N)
% Simulate data based on a linear model .
% 
% lm = A linearMixedModel that specifies the base model for a
% single sample/subject.
% effectSize = A vector with the mean of each factor and level
% (e.g. for two factors (a,b) with 2 and 3 levels respectively,
% specify the effectSize as  [a1 a2 b1 b2 b3].
% N = How many observations to generate. Note that N=1 means
% one observation for each of the combinations (i.e. 6 for the
% example).
X = designMatrix(o,lm,bf.internal.getAllTerms(lm),'zeroSumConstraint',false,'treatAsRandom',{});
X = [X{:}];
if numel(effectSize)==1
    effectSize = effectSize*ones(1,size(X,2));
end
if size(X,2) ~= numel(effectSize)
    error('Please specify effect size for the linear model as one coefficient for each column in the deisgn matrix');
end
y =  repmat(X,[N 1])*effectSize(:)+randn([N*size(X,1) 1]);

tbl = table(y,'variablenames',{lm.ResponseName});
for i=1:lm.NumPredictors
    tbl = addvars(tbl,repmat(lm.Variables.(lm.PredictorNames{i}),[N 1]),'newvariablenames',{lm.PredictorNames{i}});
end
out = {tbl,char(lm.Formula)};
end
