function [bf10,lm] = anova(x,y,varargin)
% Function to analyze an N-way ANOVA with fixed and/or random effects
%
% INPUT
%  The user can provide the data as output of fitlme (i.e. a linear
%  mixed model)
%  x = a LinearMixedModel from the stats toolbox (created by
%  fitlme).
% OR the user can provide a table and a formula
%  x = A table
%  y = A Wilkinson notation formula.
%
% Parm/Value pairs:
% 'sharedPriors' - Which columns (i.e. factors) in the table should share a
%                   their prior on their effect size. The default is that
%                   all levels within a factor share a prior.
%                   To share priors across factors, use
%                      {{'a','b'},{'c','d'}}   -> share priors for a and b
%                       and, separately, for c and d.
%                      {{'a','b','c','d'}} ->  share priors for all factors
%                      (this is the 'single g' approach in Rouder.
%                      Shortcuts:
%                       'within' - share within a fixed effect factor, not across
%                       'singleG' - share across all fixed effects.
% 'treatAsRandom' - Factors to be treated as random effects.
% 'options' - Monte Carlo Integration options. Defaults to bf.options.m
% OUTPUT
% bf10 - The Bayes Factor comparing the model to the model with intercept only.
%       To compute BF for more refined hypotheses you compute
%       a BF for the full model, and a restricted model and
%       then take the ratio. See rouderFigures for examples.
% lm  = The linear mixed model.
%
% BK -2018

if isa(x,'LinearMixedModel')
    lm = x;
    if nargin<=2
        args = {};
    else
        args = cat(2,{y},varargin); % Y must be part of the vararing
    end
elseif isa(x,'table') && ischar(y)  % Specified a table and formula
    lm = fitlme(x,y);
    args = varargin;
else
    error('bf.anova requires either a LinearMixedModel or a Table & Formula as its input');
end
p=inputParser;
p.addParameter('sharedPriors','within',@(x) ischar(x) || (iscell(x) && iscell(x{1}))); % Cell containing cells with factors(columns) that share a prior.
p.addParameter('treatAsRandom',{});
p.addParameter('options',bf.options);
p.addParameter('scale',sqrt(2)/2); 
p.parse(args{:});


f=lm.Formula;
if ~isempty(f.GroupingVariableNames)
    error('Not implemented yet');
end

allTerms = bf.internal.getAllTerms(lm);
% Construct the design matrix
[X,y] = bf.internal.designMatrix(lm,allTerms,'zeroSumConstraint',true,'treatAsRandom',p.Results.treatAsRandom);
nrAllTerms = numel(allTerms);

% Setup sharing of priors as requested
if ischar(p.Results.sharedPriors)
    switch upper(p.Results.sharedPriors)
        case 'WITHIN'
            % Share priors for each level of each factor, but not across factors
            sharedPriors = cell(1,nrAllTerms);
            [sharedPriors{:}] = deal(allTerms{:});
        case 'SINGLEG'
            sharedPriors = {allTerms};
        case {'NONE',''}
            sharedPriors ={};
    end
else
    sharedPriors = p.Results.sharedPriors;
end
sharedPriorIx = cell(1,numel(sharedPriors));
[sharedPriorIx{:}] = deal([]);

% Assign groups of effects to use the same prior on effect size.
soFar  =0;
if isempty(sharedPriors)
    sharedPriorIx = {};
else
    for i=1:numel(allTerms)
        match = cellfun(@any,cellfun(@(x)(strcmp(allTerms{i},x)),sharedPriors,'UniformOutput',false));
        if ~any(match)
            error(['Shared priors not defined for ' allTerms{i}]);
        end
        nrInThisTerm  = size(X{i},2);
        sharedPriorIx{match} = cat(2,sharedPriorIx{match},soFar+(1:nrInThisTerm));
        soFar = soFar+nrInThisTerm;
    end
end
%% Call the nWayAnova function for the actual analysis
X= [X{:}];
bf10 = bf.internal.nWayAnova(y,X,'sharedPriors',sharedPriorIx,'options',p.Results.options,'scale',p.Results.scale);
end


