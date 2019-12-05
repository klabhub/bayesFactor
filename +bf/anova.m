function [bf10,lm,lmAlternative] = anova(x,y,varargin)
% Function to analyze an N-way ANOVA with fixed and/or random effects
% Currently this function can handle all fixed effects and random intercept
% effects to model repeated measurements. Slope random effects have not 
% been implemented yet.
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
% 'alternativeModel' - The model against which the BF is to be calculated.
%                       By defaults this is the intercept-only model. For a
%                       model with random effects it is the model with only
%                       random effects. In other words, the BF quantifies
%                       the evidence that the model with all the fixed
%                       effects is better than a model with no fixed
%                       effects. 
%                       To test against a specific other model, a formula
%                       can be specified here. 
%                       For example if the main formula is rt~ori*freq, we
%                       can determine the Bayes Factor for just the main
%                       effect of orientation by specifying an alternative
%                       model without just that main effect
%                       ('rt~freq+ori:freq')
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
% 'scale' - The scale of the distribution of the prior effects. [sqrt(2)/2)
% 'randomEffectsScale' - The scale of the distribution of prior random
% effects. Random effects are typically expected to be larger than fixed
% effects so the scale defaults to sqrt(2).
%
% OUTPUT
% bf10 - The Bayes Factor comparing the model to the model with intercept
%       only or with teh model specified in 'alternativeModel'.
%       To compute BF for more refined hypotheses you compute
%       a BF for the full model, and a restricted model and
%       then take the ratio. See rouderFigures for examples.
% lm  = The linear mixed model.
%
% See gettingStarted.m for examples.
% 
% BK -2018
% Nov 2019 - Added intercept random effects for repeated measurements
% designs.

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
p.addParameter('scale',sqrt(2)/2,@isnumeric); 
p.addParameter('randomEffectsScale',sqrt(2),@isnumeric); % Wide scale by default for RE.
p.addParameter('alternativeModel','',@ischar);
p.parse(args{:});


f=lm.Formula;

%% Random effects
% Extract the formula and create a dummy lmm for each of the grouping
% variables so that we can use the same code (internal.designMatrix) to
% extract the design matrix.
% Priors are shared across all grouping variables in the same group
% For instance, (1|subject:block) + (1|day) 
% will share priors for subject:block and, separately, for all levels of day
reX = {};
reSharedPriors = {};
reTerms ={};
for grp = 1:numel(f.GroupingVariableNames)
    thisREFormula = f.RELinearFormula{grp};
    if strcmpi(thisREFormula.LinearPredictor,'1')
        %Intercept only random effect,
        thisTerms = strjoin(f.GroupingVariableNames{grp},':');
        reLm = fitlme(lm.Variables,[f.ResponseName '~ -1 + ' thisTerms]);
        [thisReX] = bf.internal.designMatrix(reLm,{thisTerms},'treatAsRandom',f.GroupingVariableNames{grp},'zeroSumConstraint',false);        
    else
        error('Slope random effects have not been implemented yet');
    end
    reX = cat(2,reX,thisReX);
    reTerms = cat(2,reTerms,{thisTerms});
    reSharedPriors = cat(2,reSharedPriors,{thisTerms});    
end
nrReTerms = numel(reTerms);
if isscalar(p.Results.randomEffectsScale)
   reScale= repmat(p.Results.randomEffectsScale,[1 nrReTerms]);
else
    reScale = p.Results.randomEffectsScale;
    assert(numel(reScale)==nrReTerms,'The number of randomEffectScales (%d) should match the number of random effects (%d) (or be a scalar)',numel(reScale),nrReTerms);
end
%% Fixed effects
feTerms = bf.internal.getAllTerms(lm);
% Construct the design matrix
[feX,y] = bf.internal.designMatrix(lm,feTerms,'zeroSumConstraint',true,'treatAsRandom',p.Results.treatAsRandom);
nrFeTerms = numel(feTerms);
% Setup sharing of priors as requested
if ischar(p.Results.sharedPriors)
    switch upper(p.Results.sharedPriors)
        case 'WITHIN'
            % Share priors for each level of each factor, but not across factors
            feSharedPriors = cell(1,nrFeTerms);
            [feSharedPriors{:}] = deal(feTerms{:});
        case 'SINGLEG'
            feSharedPriors = {feTerms};
        case {'NONE',''}
            feSharedPriors ={};
    end
else
    feSharedPriors = p.Results.sharedPriors;
end

if isscalar(p.Results.scale)
   feScale= repmat(p.Results.scale,[1 nrFeTerms]);
else
    feScale = p.Results.scale;
    assert(numel(feScale)==nrFeTerms,'The number of scales (%d) should match the number of fixed effects (%d) (or be a scalar)',numel(feScale),nrFeTerms);
end



%% Now combine FE and RE to get the BF for the full model.
X = cat(2,feX,reX);
fullSharedPriors= cat(2,feSharedPriors,reSharedPriors);
fullTerms = cat(2,feTerms,reTerms);
fullScale = [feScale reScale];
fullSharedPriorIx = bf.internal.sharedPriorIx(X,fullTerms,fullSharedPriors);
%% Call the nWayAnova function for the actual analysis
X= [X{:}];
bf10 = bf.internal.nWayAnova(y,X,'sharedPriors',fullSharedPriorIx,'options',p.Results.options,'scale',fullScale);

if ~isempty(p.Results.alternativeModel) 
    out= find(strcmpi(args(1:2:end),'AlternativeModel'));
    args([out out+1])=[];
    % Fit the alternative mdoel with same args
    [bf10Alternative,lmAlternative]= bf.anova(lm.Variables,p.Results.alternativeModel,args{:});
elseif nrReTerms >0    
    %% It there are RE, the alternative model has only the Random Effects
    reSharedPriorIx =  bf.internal.sharedPriorIx(reX,reTerms,reSharedPriors);
    bf10Alternative = bf.internal.nWayAnova(y,[reX{:}],'sharedPriors',reSharedPriorIx,'options',p.Results.options,'scale',reScale);
    lmAlternative = [];
else
    bf10Alternative =1;
    lmAlternative = [];
end
% To get the BF for the model versus the alternative we
% divide this out.
bf10 = bf10/bf10Alternative;
end


