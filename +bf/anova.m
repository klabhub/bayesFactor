function [bf10,lm,lmAlternative] = anova(x,y,varargin)
% Function to analyze an N-way ANOVA with fixed categorical and/or random effects
% as well as continuous covariates. 
% 
% Currently this function can handle all fixed effects and random intercept
% effects to model repeated measurements. Slope random effects have not 
% been implemented yet.  
% 
% For anova with only continous covariates (i.e. regression),
% the bf.bfFromR2 function is a quicker option. Use this function when
% mixing categorical and continuous covariates.
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
%                       or a cell array of formulas can be specified here. 
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
% 'options' - Monte Carlo Integration options. Defaults to bf.options
% 
% 'scale' - The scale of the distribution of the prior effects. [sqrt(2)/2)
% 'randomEffectsScale' - The scale of the distribution of prior random
% effects. Random effects are typically expected to be larger than fixed
% effects so the scale defaults to sqrt(2).
% 'continuousScale' The scale of the distribution of the prior effects for
% continuous covariates. Defaults to 1 for consistency with the
% Liang et al formula used  in bf.bfFromR2.
%
% Note on Scales: 
% Either a single scale should be specified, or one scale per group of covariates 
% that share a prior.
% 
%
% OUTPUT
% bf10 - The Bayes Factor comparing the model to the model with intercept
%       only or with teh model specified in 'alternativeModel'.
%       To compute BF for more refined hypotheses you compute
%       a BF for the full model, and a restricted model and
%       then take the ratio. See rouderFigures for examples.
% lm  = The linear mixed model.
% lmAlternative = the linear mixed model for the alternative formula
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
p.addParameter('continuousScale',1,@(x) (isscalar(x) && isnumeric(x))); % Default for continuous
p.addParameter('alternativeModel','',@(x) (ischar(x) || iscell(x)));
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
        thisReX = bf.internal.designMatrix(reLm,{thisTerms},'forceCategorical',true,'treatAsRandom',f.GroupingVariableNames{grp},'zeroSumConstraint',false);        
    else
        error('Slope random effects have not been implemented yet');
    end
    reX = cat(2,reX,thisReX);
    reTerms = cat(2,reTerms,{thisTerms});
    reSharedPriors = cat(2,reSharedPriors,{thisTerms});    
end
nrRePriors = numel(reSharedPriors);
nrReTerms = numel(reTerms);
if isscalar(p.Results.randomEffectsScale)
   reScale= cell(1,nrRePriors);
   [reScale{:}]= deal(p.Results.randomEffectsScale);% separate effects, but same scale 
else
    reScale = num2cell(p.Results.randomEffectsScale); % One scale for each effect
    assert(numel(reScale)==nrRePriors,'The number of randomEffectScales (%d) should match the number of random effects grouping variables (%d) (or be a scalar)',numel(reScale),nrRePriors);
end
%% Fixed effects
feTerms = bf.internal.getAllTerms(lm);
% Construct the design matrix
[feX,y,isContinuous] = bf.internal.designMatrix(lm,feTerms,'zeroSumConstraint',true,'treatAsRandom',p.Results.treatAsRandom);
% Handle the categorical variables first.
catFeX      = feX(~isContinuous);
catFeTerms  = feTerms(~isContinuous);
nrCatFeTerms = size(catFeX,2);
% Setup sharing of priors as requested for categorical variables, all
% continuous variables share one Zellner-Siow prior
if ischar(p.Results.sharedPriors)
    switch upper(p.Results.sharedPriors)
        case 'WITHIN'
            % Share priors for each level of each factor, but not across factors
            catFeSharedPriors = cell(1,nrCatFeTerms);
            [catFeSharedPriors{:}] = deal(catFeTerms{:});
        case 'SINGLEG'
            catFeSharedPriors = {catFeTerms};
        case {'NONE',''}
            catFeSharedPriors ={};
    end
else
    catFeSharedPriors = p.Results.sharedPriors;
end
nrCatFeSharedPriors = numel(catFeSharedPriors);
if isscalar(p.Results.scale)
   catFeScale= cell(1,nrCatFeSharedPriors);
   [catFeScale{:}] =deal(p.Results.scale);
else
    catFeScale = num2cell(p.Results.scale);
    assert(numel(catFeScale)==nrCatFeSharedPriors,'The number of scales (%d) should match the number of fixed effects (%d) (or be a scalar)',numel(catFeScale),nrCatFeSharedPriors);
end

  

%% Now combine FE and RE .
X = cat(2,catFeX,reX);
fullSharedPriors= cat(2,catFeSharedPriors,reSharedPriors);
fullTerms = cat(2,catFeTerms,reTerms);
fullScale = cat(2,catFeScale,reScale);
fullSharedPriorIx = bf.internal.sharedPriorIx(X,fullTerms,fullSharedPriors);


%% Add the continuous effects
% These always share a single g, but with a scale factor that is determined
% by their (co)-variance
if any(isContinuous)
    contIx = size(X,2)+(1:sum(isContinuous));     
    X = cat(2,X,feX(isContinuous));
    fullScale = cat(2,fullScale,{p.Results.continuousScale});
    fullSharedPriorIx = cat(2,fullSharedPriorIx, {contIx}); 
else
    contIx =[];
end
%% Call the nWayAnova function for the actual analysis
X= [X{:}];
bf10 = bf.internal.nWayAnova(y,X,'sharedPriors',fullSharedPriorIx,'options',p.Results.options,'scale',fullScale,'continuous',contIx);

if ~isempty(p.Results.alternativeModel) 
    % A specific alternative model formula was specified
    if ischar(p.Results.alternativeModel)
        alternativeModel = {p.Results.alternativeModel};
    else
       alternativeModel = p.Results.alternativeModel;
    end
    out= find(strcmpi(args(1:2:end),'AlternativeModel'));
    args([out out+1])=[];
    nrAlternatives = numel(alternativeModel);
    bf10Alternative = nan(nrAlternatives,1);    
    lmAlternative =  cell(nrAlternatives,1);
    for alt = 1:nrAlternatives 
        % Fit the alternative mdoel with same args
        [bf10Alternative(alt),lmAlternative{alt}]= bf.anova(lm.Variables,alternativeModel{alt},args{:});
    end
elseif nrReTerms >0    
    %% It there are RE, the alternative model has only the Random Effects
    reSharedPriorIx =  bf.internal.sharedPriorIx(reX,reTerms,reSharedPriors);
    bf10Alternative = bf.internal.nWayAnova(y,[reX{:}],'sharedPriors',reSharedPriorIx,'options',p.Results.options,'scale',reScale);
    lmAlternative = [];
else % No alternative model specified , compare to the Intercept only model
    bf10Alternative =1;
    lmAlternative = [];
end
% To get the BF for the model versus the alternative we
% divide this out.
bf10 = bf10./bf10Alternative;
end


