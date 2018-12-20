function [bf10,model] = lme(tbl,response,fixedEffects,randomEffects,varargin)
% Analyze table data with a linear mixed effects model.
%
% INPUT
% tbl = The table with data.
% response = name of the column that represents the response
% fixedEffects = cell array of columns that should be modeled as fixed
%                   effects
% randomEffects = cell array of columns that should be modeled as random
%                   effects.
%
% Parm/Value pairs:
% 'encoding'  - How to encode the different levels of the model.
%               'effects','reference','referenceLast' are standard dummy
%               var encoding options (see LinearMixedModel), while 'Rouder'
%               is a variant to implement the zero--sum constraint for
%               fixed effects that gives each level equal marginal prior
%               (and does not favor a specific level). This is the default.
%
% 'sharedPriors' - Which columns (i.e. factors) in the table should share a
%                   their prior on their effect size. The default is that
%                   all levels within a factor share a prior.
%                   To share priors across factors, use
%                      {'{'a','b'},{'c','d'}}   -> share priors for a and b
%                       and, separately, for c and d.
%                      {{'a','b','c','d'}} ->  share priors for all factors
%                      (this is the 'single g' approach in Rouder.
%                      Shortcuts:
%                       'within' - share within a fixed effect factor, not across
%                       'singleG' - share across all fixed effects.
% 'interactions'    - Include interactions ('all') or not ('none'). Default
%                       is  'none'.  Or specify a set {'a:b','c:d'}
% 'nDimsForMC'  - The number of dimensions at which we start to use MC
% integration. [4] (3D can be quite slow for quadrature integration, so you
% could reduce this to 3 to speed things up, at little cost in accuracy).
%
%
% BK -2018

p=inputParser;
p.addParameter('interactions','none',@(x) (ischar(x) || iscell(x)));
p.addParameter('encoding','rouder',@(x) (ischar(x) && ismember(lower(x),{'effects','reference','referenceLast','rouder'})))
p.addParameter('sharedPriors','within',@(x) ischar(x) || (iscell(x) && iscell(x{1}))); % Cell containing cells with factors(columns) that share a prior.
p.addParameter('nDimsForMC',4,@(x)(x<=4)); % Passed to bf.anova
p.parse(varargin{:});


nrFixedEffects =numel(fixedEffects);
nrRandomEffects =numel(randomEffects);

if nrRandomEffects>0
    warning('Random effects modeling does not seem right... see testSuite -rouderFigure5');
end
%% Setup list of interactions
if ischar(p.Results.interactions)
    switch upper(p.Results.interactions)
        case 'FIXED'
            interactions = bf.allInteractions(fixedEffects);
        case 'RANDOM'
            interactions = bf.allInteractions(randomEffects);
        case 'ALL'
            interactions = bf.allInteractions(cat(2,fixedEffects,randomEffects));
        case 'NONE'
            interactions = {};
    end
else
    interactions = p.Results.interactions;
end

nrInteractions = numel(interactions);
allTerms = cat(2,fixedEffects,randomEffects,interactions);
nrAllTerms = nrFixedEffects+nrRandomEffects+nrInteractions;

%% Setup sharing of priors as requested
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

%% Construct the design matrix for the fixed effects
designMatrix = cell(1,nrFixedEffects+nrRandomEffects +nrInteractions);
allMain = cat(2,fixedEffects,randomEffects);
for i=1:nrFixedEffects+nrRandomEffects
    thisX = classreg.regr.modelutils.designmatrix(tbl,'model','linear','intercept',false,'DummyVarCoding','full','PredictorVars',allMain{i},'responseVar','');
    if ismember(allMain{i},randomEffects)
        % Random effect - Keep full dummy X
    else
        switch upper(p.Results.encoding)
            case 'REFERENCE'
                % Treatment contrast: effects defined relative to the first level
                thisX =thisX(:,2:end);
            case 'REFERENCELAST'
                % Treatment contrast: effects defined relative to the first level
                thisX =thisX(:,1:end-1);
            case 'EFFECTS'
                % Relative to mean with -1 for last level.
                last = thisX(:,end)==1;
                thisX =thisX(:,1:end-1);
                thisX(last,:) = -1;
            case 'ROUDER'
                % Sum-to-zero contrasts that equates marginal priors across levels.
                % (Rouder 2019)
                [~,thisX] = bf.fixedEffectConstraint(thisX); 
        end
    end
    designMatrix{i} =thisX; % Store for later use
end

%% Construct the interaction parts of the design matrix for all requested interactions
for i=1:nrInteractions
    aName =extractBefore(interactions{i},':');
    bName = extractAfter(interactions{i},':');
    a =ismember(allMain,aName);
    b =ismember(allMain,bName);
    
    thisA = classreg.regr.modelutils.designmatrix(tbl,'model','linear','intercept',false,'DummyVarCoding','full','PredictorVars',aName,'responseVar','');
    thisB = classreg.regr.modelutils.designmatrix(tbl,'model','linear','intercept',false,'DummyVarCoding','full','PredictorVars',bName,'responseVar','');
    if ismember(aName,fixedEffects)
        [~,thisA] = bf.fixedEffectConstraint(thisA); 
    end
    if ismember(bName,fixedEffects)
        [~,thisB] = bf.fixedEffectConstraint(thisB); 
    end
    designMatrix{nrFixedEffects+nrRandomEffects+i} = bf.interaction(thisA,thisB);          
end


%% Assign groups of effects to use the same prior on effect size.
soFar  =0;
if isempty(sharedPriors)
    sharedPriorIx = {};
else
    for i=1:numel(allTerms)
        match = cellfun(@any,cellfun(@(x)(strcmp(allTerms{i},x)),sharedPriors,'UniformOutput',false));
        if ~any(match)
            error(['Shared priors not defined for ' allTerms{i}]);
        end
        nrInThisTerm  = size(designMatrix{i},2);
        sharedPriorIx{match} = cat(2,sharedPriorIx{match},soFar+(1:nrInThisTerm));
        soFar = soFar+nrInThisTerm;
    end
end
%% Call the anova function for the actual analysis
bf10 = bf.anova(tbl.(response),[designMatrix{:}],'sharedPriors',sharedPriorIx,'nDimsForMC',p.Results.nDimsForMC);


