function [results] = modmed(T,varargin)
% [results] = modmed(T,varargin)
% 
% Moderated Mediation and Mediated Moderation, following the definitions in
% Muller D, Judd CM, Yzerbyt VY (2005) When moderation is mediated and mediation is moderated. J Pers Soc Psychol 89:852–863.
% and the recommended definition of effect size in 
% Preacher KJ, Kelley K (2011) Effect size measures for mediation models:
% Quantitative strategies for communicating indirect effects. Psychol Methods 16:93–115.
%
% Linear models are fit using fitlme's standard options (MLE). 
%
% If a moderator is included, paths are evaluated at +/- 1 stdev of the
% (continuous) moderator, or at the first and last categorical value of the
% categorical moderator.
%
% For models including random effects, the effect size kappa2 is computed
% after removing the fitted random effects.
% 
% SEE ALSO lm.plotModMed for a graphical representation of the results.
%
% INPUT
%
% T - Data table.
% 'outcome' -  the name of the column that stores the outcome (dependent)
% variable.
% 'treatment' - name of the column that stores the treatment (independent) variable.
% mediator - name or cell array of names of columns with candidate
%               mediators.
% moderator - name of the column that is a candidate moderator (can be
%               empty).
% randomEffects - Formula to specify random effects (e.g. (1|subject) to
%                   allow a random intercept per subject.
% alpha - significance level to use (0.05)
% dummyVarCoding - Coding to use for categorical variables (Effects)
% bootstrap - Number of bootstrap samples to use to estimate confidence
%               intervals.
% confLimits - Percentiles of the confidence limits to determine [2.5 97.5]
%
% OUTPUT
% Indirect paths are evaluated for each mediator at two values
% of the moderator [nrMediators 2];

% results.a = a path (Treatment->Mediator)
%        .b = b path (Mediator->Outcome_
%       .ab = indirect path (=product of a and b)
%       .c  = overall effect (treatment->outcome regression without mediators).
%       .cPrime =  direct effect (treatment->outcome  after accounting for
%                    mediators)
%       .kappa2 = Effect size; ranges between 0 and 1: the proportion of the maximum possible indirect effect that could have occurred. 
%       .clim = bootstrapped confidence limits for each of the parameters.
%       .parms = parameters passed to this function.
%       .lm6 =  Linear Mixed Model of the full regression (Eq6 in Muller et al)
% The .style struct interprets the mediation analysis outcomes following
% the rules of Muller et al. This is based on the significance of paths 
%       .style.modMed = Logical [1 nrMediators]. True for each mediator that is determined to be a moderated mediator
%       .style.medMod = Logical [1 nrMediators]. True for each mediator that mediates the moderator.
%       .style.prototypical = Logical [1 nrMediators]. True for a moderator
%          that does not interact significantly with the treatment in
%           moderated mediation. (i.e. the moderation is only significant with
%           the mediator included). Muller et all call this prototypical
%           moderated mediation.
%       .style.med = Logical [1 nrMediators]. True for each mediator (only assessed if no moderator is specified).
%
% EXAMPLE
%
% Mediated Moderation and Moderated Mediation
% Example data from Muller et al are stored in the examples folder.
%
% medmodT = readtable('yzerbyt.medmod.csv'); 
%  Run a small number of bootstraps to save time. 
% results =  lm.modmed(medmodT,'treatment','PRIME','moderator','SVO','mediator','EXP','outcome','BEH','bootstrap',100);
% Show results with kappa2
% lm.plotModMed(results,'effectSize','kappa2');
%
% modMedT = readtable('yzerbyt.modmed.csv');
% results = lm.modmed(modMedT,'treatment','MOOD','moderator','NFC','mediator','POS','outcome','ATT','bootstrap',100);
% lm.plotModMed(results,'effectSize','kappa2');
% 
% BK -    August 2021

p= inputParser;
p.addParameter('outcome','',@ischar);
p.addParameter('treatment','',@ischar);
p.addParameter('mediator','',@(x) ischar(x) || iscellstr(x));
p.addParameter('moderator','',@ischar);
p.addParameter('randomEffects','',@ischar);
p.addParameter('alpha',0.05,@isnumeric);
p.addParameter('dummyVarCoding','effects',@ischar);
p.addParameter('bootstrap',1,@isnumeric);
p.addParameter('confLimits',[2.5 97.5],@isnumeric);
p.parse(varargin{:})

if ischar(p.Results.mediator)
    mediators = {p.Results.mediator};
else
    mediators = p.Results.mediator;
end
nrMediators = numel(mediators);

% Construct the formulas for equations 4,5,6
eq4 = [p.Results.outcome '~' p.Results.treatment];
if isempty(p.Results.moderator)
    hasModerator = false;
else
    hasModerator = true;
    eq4 = [eq4  ' + ' p.Results.moderator '+' p.Results.moderator ':' p.Results.treatment];
end
eq5= cell(1,nrMediators);
for i=1:nrMediators
    eq5{i} = [mediators{i} '~' p.Results.treatment];
    if hasModerator
        eq5{i}= [eq5{i}  ' + ' p.Results.moderator '+' p.Results.moderator ':' p.Results.treatment];
    end
    if ~isempty(p.Results.randomEffects)
        eq5{i}  = [eq5{i} '+' p.Results.randomEffects];
    end
end

eq6 = [p.Results.outcome '~' p.Results.treatment ];
if hasModerator
    eq6 = [eq6  ' + ' p.Results.moderator '+' p.Results.moderator ':' p.Results.treatment];
end
for i=1:nrMediators
    eq6 = [eq6 ' + ' mediators{i}]; %#ok<AGROW>
    if hasModerator
        eq6 = [eq6 ' + ' mediators{i} ':' p.Results.moderator]; %#ok<AGROW>
    end
end


if ~isempty(p.Results.randomEffects)
    eq4  = [eq4 '+' p.Results.randomEffects];
    eq6  = [eq6 '+' p.Results.randomEffects];
end

[results.a,results.b,results.c,results.cPrime,results.ab,results.kappa2,results.moderatorValues, results.style, results.lm6] = locRegression(T,eq4,eq5,eq6,p.Results.treatment,mediators,p.Results.moderator,p.Results.dummyVarCoding,p.Results.alpha);
%% Boostrap confidence limits by resampling the rows in T (trials, presumably)
if hasModerator
    nrModeratorValues= size(results.a,2);
else
    nrModeratorValues= 1;
end
if p.Results.bootstrap>0
    bsA = nan([nrMediators nrModeratorValues p.Results.bootstrap]);
    bsB =  nan([nrMediators nrModeratorValues p.Results.bootstrap]);
    bsC = nan([1 nrModeratorValues p.Results.bootstrap]);
    bsCPrime = nan([1, nrModeratorValues p.Results.bootstrap]);
    bsKappa2 = nan([nrMediators nrModeratorValues p.Results.bootstrap]);
    bsAB = nan([nrMediators nrModeratorValues p.Results.bootstrap]);
    for bs = 1:p.Results.bootstrap
        ix= randi(height(T),height(T),1); % Resample with replacement
        [bsA(:,:,bs),bsB(:,:,bs),bsC(:,:,bs),bsCPrime(:,:,bs),bsAB(:,:,bs),bsKappa2(:,:,bs)] = locRegression(T(ix,:),eq4,eq5,eq6,p.Results.treatment,mediators,p.Results.moderator,p.Results.dummyVarCoding,p.Results.alpha);
    end
    % Determine specified percentiles
    results.clim.a = prctile(bsA,p.Results.confLimits,3);
    results.clim.b = prctile(bsB,p.Results.confLimits,3);
    results.clim.c = prctile(bsC,p.Results.confLimits,3);
    results.clim.cPrime = prctile(bsCPrime,p.Results.confLimits,3);
    results.clim.kappa2 = prctile(bsKappa2,p.Results.confLimits,3);
    results.clim.ab = prctile(bsAB,p.Results.confLimits,3);
    % And store a range of percentiles just in case...
    step = 100/(p.Results.bootstrap/10); % 1% bins for 1000 bs, 0.1% for 10000
    bins = 0:step:100;
    results.bs.a =prctile(bsA,bins,3);
    results.bs.b =prctile(bsB,bins,3);
    results.bs.c =prctile(bsC,bins,3);
    results.bs.cPrime =prctile(bsCPrime,bins,3);
    results.bs.kappa2 =prctile(bsKappa2,bins,3);
    results.bs.ab =prctile(bsAB,bins,3);    
else
    results.clim = [];
end
results.parms = p.Results;

end

%%
function [a,b,c,cPrime,ab,kappa2,moderatorValue,style,lm6] = locRegression(T,eq4,eq5,eq6,treatment,mediators,moderator,dummyVar,alpha)
%% Fit the linear models and extract betas and p-valus using the nomenclature of the paper.
nrMediators = numel(mediators);
hasModerator = ~isempty(moderator);
lm4 = fitlme(T,eq4,'DummyVarCoding',dummyVar);
for i=1:nrMediators
    lm5{i} = fitlme(T,eq5{i},'DummyVarCoding',dummyVar); %#ok<AGROW>
end
lm6 = fitlme(T,eq6,'DummyVarCoding',dummyVar);

% Extract from models
is41 = strcmpi(lm4.anova.Term,treatment);
if hasModerator
    is42 = strcmpi(lm4.anova.Term,moderator);
    is43 = strcmpi(lm4.anova.Term,[treatment ':' moderator])| strcmpi(lm4.anova.Term,[moderator ':' treatment ]);
end
is51 = cell(1,nrMediators);
is52 = cell(1,nrMediators);
is53 = cell(1,nrMediators);
for i=1:nrMediators
    is51{i} = strcmpi(lm5{i}.anova.Term,treatment);
    if hasModerator
        is52{i} = strcmpi(lm5{i}.anova.Term,moderator);
        is53{i} = strcmpi(lm5{i}.anova.Term,[treatment ':' moderator]) | strcmpi(lm5{i}.anova.Term,[moderator ':' treatment]);
    end
end

is61 = strcmpi(lm6.anova.Term,treatment);
is64 = cell(1,nrMediators);
is65 = cell(1,nrMediators);
for i=1:nrMediators
    is64{i} = strcmpi(mediators{i},lm6.anova.Term);
    if hasModerator
        is65{i} = strcmpi(lm6.anova.Term,[mediators{i} ':' moderator]) | strcmpi(lm6.anova.Term,[moderator ':' mediators{i}]);
    end
end
if hasModerator
    is62 = strcmpi(lm6.anova.Term,moderator);
    is63 = strcmpi(lm6.anova.Term,[treatment ':' moderator]) | strcmpi(lm6.anova.Term,[moderator ':' treatment]);
end


b41 = lm4.Coefficients.Estimate(is41);
p41   = lm4.Coefficients.pValue(is41);
if hasModerator
    b42 = lm4.Coefficients.Estimate(is42);    
    b43 = lm4.Coefficients.Estimate(is43);
    p42  = lm4.Coefficients.pValue(is42);
    p43  = lm4.Coefficients.pValue(is43);
else
    b42 = 0;
    b43 = 0;
    p42 = 1;
    p43 =1;
end

for i=1:nrMediators
    b51(i) = lm5{i}.Coefficients.Estimate(is51{i});
    p51(i)   = lm5{i}.Coefficients.pValue(is51{i});
    b64(i) = lm6.Coefficients.Estimate(is64{i});
    p64(i)  = lm6.Coefficients.pValue(is64{i});
    if hasModerator
        b52(i) = lm5{i}.Coefficients.Estimate(is52{i});
        b53(i) = lm5{i}.Coefficients.Estimate(is53{i});           
        p52(i)  = lm5{i}.Coefficients.pValue(is52{i});
        p53(i)  = lm5{i}.Coefficients.pValue(is53{i});    
        b65(i) = lm6.Coefficients.Estimate(is65{i});
        p65(i)  = lm6.Coefficients.pValue(is65{i});
    else 
        b53(i) = 0;
        b65(i) = 0;
        p53(1)  = 1;
        p65(i) = 1;
    end
    
end
b61 = lm6.Coefficients.Estimate(is61);
p61   = lm6.Coefficients.pValue(is61);
if hasModerator
    b62 = lm6.Coefficients.Estimate(is62);
    b63 = lm6.Coefficients.Estimate(is63);
    p62  = lm6.Coefficients.pValue(is62);
    p63  = lm6.Coefficients.pValue(is63);
else
    b62 = 0;
    b63 = 0;
    p62 = 1;
    p63 = 1;
end




isMed =false(1,nrMediators);
if hasModerator
    %% Assess mediated moderation
    overallTreatmentModeration = p43<alpha;
    for i=1:nrMediators
        medModPattern1 = p53(i) < alpha && p64(i) < alpha; %Mod affects treatment effect of mediator  && mediator affects outcome
        medModPattern2 = p51(i) < alpha && p65(i) < alpha; %Treatment affects mediator && moderator affects the mediator's effect on the outcome
        isMedMod(i) = overallTreatmentModeration && ( medModPattern1 || medModPattern2);
        
        %% Asses moderated mediation
        overalTreatmentEffect = p41< alpha;
        modMedPattern1 = p53(i) < alpha && p64(i) < alpha; % Moderaotr affects how mediator depends on treatment && mediator affects outcome
        modMedPattern2 = p65(i)< alpha && p51(i) <alpha; % moderator affects how mediator affects outcome &&  treatment affects mediator.
        isPrototypical = p43 > alpha; % no interaction between treatment and mdoerator.
        isModMed(i)= overalTreatmentEffect && ( modMedPattern1 || modMedPattern2);
        isModMedPrototypical(i) = isModMed(i) && isPrototypical;
    end
else
    %% Asses mediation
    overalTreatmentEffect = p41< alpha;
    isModMedPrototypical = false(1,nrMediators);
    isMedMod = false(1,nrMediators);
    isModMed = false(1,nrMediators);        
    for i=1:nrMediators
        treatmentAffectsMediator = p51(i) < alpha;
        mediatorReducesTreatment  = p64(i) < alpha && abs(b61) < abs(b41); % The latter term needs statistical evaluation...
        isMed(i) = overalTreatmentEffect && treatmentAffectsMediator && mediatorReducesTreatment;
    end
end

style.modMed = isModMed;
style.medMod = isMedMod;
style.med    = isMed;
style.prototypical = isModMedPrototypical;

%% Calculate a,b,c,c'
% In the presence of moderators, we compute these at the lowest and highest
% categorical value of the moderator , or at +/- 1 stdev of the moderator
% for continuous moderators.
if hasModerator
    if strcmpi(lm6.VariableInfo{moderator,'Class'},'double')
        % Continious variable as moderator
        sd = std(T.(moderator));
        m = mean(T.(moderator));
        upperModerator = m +sd;
        lowerModerator = m-sd;
    elseif strcmpi(lm6.VariableInfo{moderator,'Class'},'categorical')
        % need to check the designmatrix as the number depends on the
        % coding of dummyvars
        dm = lm6.designMatrix;
        lowerModerator = min(dm(:,is62));
        upperModerator = max(dm(:,is62));
    end
    moderatorValue = [lowerModerator upperModerator];
else
    moderatorValue  =1;
end

% X ---- c ---> Y
% X --- > a ---> Me
% Me ---> b ----> Y
nrModeratorValues = numel(moderatorValue);
c = b41+b43.*moderatorValue;
a =nan(nrMediators,nrModeratorValues);
b =nan(nrMediators,nrModeratorValues);
for i=1:nrMediators
    a(i,:) =  b51(i)+b53(i).*moderatorValue;
    b(i,:) =  b64(i)+b65(i).*moderatorValue;
end
cPrime = b61+b63.*moderatorValue;
ab = a.*b;
%% Effect size kappa2
% Preacher KJ, Kelley K (2011) Effect size measures for mediation models:
% Quantitative strategies for communicating indirect effects. Psychol Methods 16:93–115.

% Remove the fitted random effects from the response
reDM = designMatrix(lm6,'Random');
if ~isempty(reDM)
    y  = lm6.response - reDM*lm6.randomEffects;
else
    y = lm6.response;
end
S = cov([y lm6.designMatrix],'omitrows');
isX = find(strcmpi(lm6.anova.Term,treatment))+1;
isY = 1;
sxx = S(isX,isX);
syx = S(1,isX);
syy = S(1,1);
kappa2 = nan(nrMediators,nrModeratorValues);
for i=1:nrMediators
    isMe = find(strcmpi(lm6.anova.Term,mediators{i}))+1;
    smx = S(isX,isMe);
    smm = S(isMe,isMe);
    sym = S(isY,isMe);
    tmp = sqrt(smm*syy-sym^2)*sqrt(sxx*syy-syx^2);
    aLim = [sym*syx-tmp, sym*syx+tmp]./(sxx*syy);
    bLim = [-1 1].* sqrt(sxx*syy-syx^2)/sqrt(sxx*smm-smx^2);
    for mo=1:nrModeratorValues
        extremeA = aLim(sign(aLim)==sign(a(i,mo)));
        extremeB = bLim(sign(bLim)==sign(b(i,mo)));
        kappa2(i,mo) = a(i,mo).*b(i,mo)./(extremeA*extremeB);
    end
end
end


