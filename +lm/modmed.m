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
% nuisance - Formula to specify nuisance variables to include in each model
%
% alpha - significance level to use (0.05)
% dummyVarCoding - Coding to use for categorical variables (Effects)
% exclude - Rows of the datatable to exclude from all models.
% bootstrap - Number of bootstrap samples to use to estimate confidence
%               intervals.
%               The default is to resample rows of the data table T.
%               Presumably these are observations (e.g. trials). To
%               resample subjects or some unique combination of column
%               values in the table, use groupResampling.
% groupResampling - Defaults to ''. Used only for the resampling of
% bootstrap sets. Set this to a column name to resample (with replacement)
% unique values of that column in the table. Use a cell array of strings to
% combine multiple columns. For instance, 'subjectNr' to resample all rows
% corresponding to a unique subjectNr, or {'subjectNr','runNr'} to resample
% runs per subject.
%
% matchMissing - [true] If the data set has missing values, the models for
% the overall outcome and mediator could be based on different subdatasets
% (depedning on what is missing). By setting this to false, each model uses
% all the data that are available. By setting this to false, all models use
% the same data (i.e. the minimal set of rows that has no missing data in
% any of the relevant parameters).
% nrWorkers -  how many workers to use in a parpool for the bootstrapping.
% [1].
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
%       .rho   = Correlation ebtween the residuals of the
%       treatment->mediator and treatment->outcom (full model). The
%       assumption of the analysis is that this is zero. ("Sequential
%       Ignorability")
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
p.addParameter('mediator','',@(x) ischar(x) || iscellstr(x) || isstring(x));
p.addParameter('earlyMediatorSelection',false);
p.addParameter('moderator','',@ischar);
p.addParameter('nuisance','',@(x) (isempty(x) || ischar(x)));
p.addParameter('randomEffects','',@ischar);
p.addParameter('alpha',0.05,@isnumeric);
p.addParameter('dummyVarCoding','effects',@ischar);
p.addParameter('exclude',[],@islogical);
p.addParameter('matchMissing',true,@islogical);
p.addParameter('bootstrap',1,@isnumeric);
p.addParameter('groupResampling','',@(x) (iscellstr(x) || ischar(x) ||isstring(x)));
p.addParameter('nrWorkers',1,@isnumeric);
p.parse(varargin{:})

if ischar(p.Results.mediator)
    mediators = {p.Results.mediator};
else
    mediators = p.Results.mediator;
end
nrMediators = numel(mediators);
hasModerator = ~isempty(p.Results.moderator);
hasNuisance = ~isempty(p.Results.nuisance);
hasRandomEffects = isempty(p.Results.randomEffects);

if ~isempty(p.Results.exclude)
    T= T(~p.Results.exclude,:);
end

if p.Results.matchMissing
    % Check for missing observations (rotws in T) so that each model will be fit to the same
    % dataset.
    missing = false(height(T),1);
    missing = missing | ismissing(T.(p.Results.outcome)) | ismissing(T.(p.Results.treatment)) ;
    if hasModerator
        missing = missing | ismissing(T.(p.Results.moderator));
    end

    for i=1:nrMediators
        missing = missing | ismissing(T.(mediators{i}));
    end

    if any(missing)
        fprintf('Removing %d (%3.2f%%) missing values\n',sum(missing),100*mean(missing));
        T(missing,:) =[];
    end
end

% The variable naming below follows the numbering of the Muller et al
% paper, with eq4 the  overall effect, eq5 models how the treatment affects
% the mediators, and eq6 is the joint model of mediators moderators and
% treatment affect outcome.
% To make the code a bit more readable without the paper I also use X to
% stand for the treatment, Me for mediators, Mo for moderators, so
% pX41 is the p-value of the treatment effect in equation 4 (called 41 in
% Muller et al) and
% bMeMo refers to an interaction beta for Mediators and Moderator.

% Construct the formulas for equations 4,5,6
% Eq 4: treatment ->outcome
eq4 = [p.Results.outcome '~' p.Results.treatment];
if hasModerator
    eq4 = [eq4  ' + ' p.Results.moderator '+' p.Results.moderator ':' p.Results.treatment];
end
if hasNuisance
    eq4  = [eq4 '+' p.Results.nuisance];
end
if hasRandomEffects
    eq4  = [eq4 '+' p.Results.randomEffects];
end



% Eq 5: treatment -> mediator
eq5= cell(1,nrMediators);
for i=1:nrMediators
    eq5{i} = [mediators{i} '~' p.Results.treatment];
    if hasModerator
        eq5{i}= [eq5{i}  ' + ' p.Results.moderator '+' p.Results.moderator ':' p.Results.treatment];
    end
    if hasNuisance
        eq5{i}  = [eq5{i} '+' p.Results.nuisance];
    end
    if hasRandomEffects
        eq5{i}  = [eq5{i} '+' p.Results.randomEffects];
    end
end


if p.Results.earlyMediatorSelection
    % Throw out non-significant candidate mediators before bootstrapping
    % and before constructing equation 6. This has the advantage that the
    % eq6 will potentially be less complex, with more power to find effects
    % in the mediators that are significnatly modulated by the treatment.
    % This code basically duplicates the eq5 code in the locRegression
    % function
    keepMediator = nan(1,nrMediators);
    for i=1:nrMediators
        tmpLm5 = fitlme(T,eq5{i},'DummyVarCoding','effects');
        isX51 = strcmpi(tmpLm5.anova.Term,p.Results.treatment);
        keepMediator(i) = tmpLm5.anova.pValue(isX51)<p.Results.alpha;
        if hasModerator && ~keepMediator(i)
            % Also keep if there is a significant interaction with the
            % treatment
            isXMo53 = strcmpi(tmpLm5.anova.Term,[p.Results.treatment ':' p.Results.moderator]) | strcmpi(tmpLm5.anova.Term,[p.Results.moderator ':' p.Results.treatment]);
            keepMediator(i) = keepMediator(i) || tmpLm5.anova.pValue(isXMo53);
        end
    end
else
    % Keep all for now.
    keepMediator = true(1,nrMediators);
end


% Eq6 : (treatment, mediators)->outcome
eq6 = [p.Results.outcome '~' p.Results.treatment ];
if hasModerator
    eq6 = [eq6  ' + ' p.Results.moderator '+' p.Results.moderator ':' p.Results.treatment];
end
for i=1:nrMediators
    if keepMediator(i)
        eq6 = [eq6 ' + ' mediators{i}]; %#ok<AGROW>
        if hasModerator
            eq6 = [eq6 ' + ' mediators{i} ':' p.Results.moderator]; %#ok<AGROW>
        end
    end
end
if hasNuisance
    eq6  = [eq6 '+' p.Results.nuisance];
end
if hasRandomEffects
    eq6  = [eq6 '+' p.Results.randomEffects];
end

% Compute in a sub to simplify bootstratpping
[results.a,results.b,results.c,results.cPrime,results.ab,results.kappa2,results.moderatorValues, results.style, results.rho, results.lm6] = locRegression(T,eq4,eq5,eq6,p.Results.treatment,mediators,p.Results.moderator,p.Results.dummyVarCoding,p.Results.alpha,keepMediator);
%% Boostrap confidence limits by resampling the rows in T (trials, presumably)
if hasModerator
    nrModeratorValues= size(results.a,2);
else
    nrModeratorValues= 1;
end
nrBs = p.Results.bootstrap;
if nrBs>0

    bsA = nan([nrMediators nrModeratorValues nrBs]);
    bsB =  nan([nrMediators nrModeratorValues nrBs]);
    bsC = nan([1 nrModeratorValues nrBs]);
    bsCPrime = nan([1, nrModeratorValues nrBs]);
    bsKappa2 = nan([nrMediators nrModeratorValues nrBs]);
    bsAB = nan([nrMediators nrModeratorValues nrBs]);
    treatment = p.Results.treatment;
    moderator= p.Results.moderator;
    dummyVarCoding=p.Results.dummyVarCoding;
    alpha = p.Results.alpha;
    groupResampling = p.Results.groupResampling;
    if ~isempty(groupResampling)
        uGrpT =unique(T(:,p.Results.groupResampling),'rows');
    else
        uGrpT = NaN; % Not used but parfor wants this defined.
    end

   % Use a dataqueue to show progress updates
   dataQueue = parallel.pool.DataQueue;
   hWaithBar = waitbar(0,'Mediation bootstrapping');
    cntr =0;
    

    afterEach(dataQueue,@updateWaitBar);

    parfor (bs = 1:nrBs, p.Results.nrWorkers)
    %for (bs = 1:nrBs)
        if isempty(groupResampling)
            % Resample rows (trials probably) with replacement
            ix= randi(height(T),height(T),1);
            thisT  = T(ix,:);
        else
            % Resample uniue groups (subjects or some combo of subject and run)
            grpIx = randi(height(uGrpT),height(uGrpT),1);
            thisGrpT = uGrpT(grpIx,:);
            thisT = innerjoin(thisGrpT,T,'Keys',groupResampling);
        end
        [bsA(:,:,bs),bsB(:,:,bs),bsC(:,:,bs),bsCPrime(:,:,bs),bsAB(:,:,bs),bsKappa2(:,:,bs)] = locRegression(thisT,eq4,eq5,eq6,treatment,mediators,moderator,dummyVarCoding,alpha,keepMediator);
        % Update the waitbar      
        send(dataQueue,bs);
    end
    close(hWaithBar);

    % Determine specified percentiles
    confLimits = 100*[p.Results.alpha/2 1-p.Results.alpha/2];
    results.clim.a = prctile(bsA,confLimits,3);
    results.clim.b = prctile(bsB,confLimits,3);
    results.clim.c = prctile(bsC,confLimits,3);
    results.clim.cPrime = prctile(bsCPrime,confLimits,3);
    results.clim.kappa2 = prctile(bsKappa2,confLimits,3);
    results.clim.ab = prctile(bsAB,confLimits,3);
    % And store a range of percentiles just in case...
    step = 100/(nrBs/10); % 1% bins for 1000 bs, 0.1% for 10000
    bins = 0:step:100;
    results.bs.bins = bins;
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

    function updateWaitBar(~)
        cntr=  cntr+1;
        waitbar(cntr/nrBs,hWaithBar);
    end

end
    
%%
function [a,b,c,cPrime,ab,kappa2,moderatorLevels,style,rho,lm6] = locRegression(T,eq4,eq5,eq6,treatment,mediators,moderator,dummyVar,alpha,keepMediator)
%% Fit the linear models and extract betas and p-valus using the nomenclature of the paper.
nrMediators = numel(mediators);
hasModerator = ~isempty(moderator);

%% Extract p-values
% Because we're using the ANOVA output to assess significance of terms, we
% have to run the lme with effects coding.

% Eq 4: Overall treatment-> outcome effect (with or without moderator)
lm4 = fitlme(T,eq4,'DummyVarCoding','effects');
isX41 =  strcmpi(lm4.anova.Term,treatment);
pX41 = lm4.anova.pValue(isX41);
isXMo43 =  strcmpi(lm4.anova.Term,[treatment ':' moderator])| strcmpi(lm4.anova.Term,[moderator ':' treatment ]);
pXMo43  = lm4.anova.pValue(isXMo43);


% Eq 5 treatment -> mediator effect (withor without moderator)
lm5 = cell(1,nrMediators);
pX51 = nan(1,nrMediators);
pXMo53=nan(1,nrMediators*hasModerator);
for i=1:nrMediators
    if keepMediator(i)
        lm5{i} = fitlme(T,eq5{i},'DummyVarCoding','effects');
        isX51 = strcmpi(lm5{i}.anova.Term,treatment);
        pX51(i) = lm5{i}.anova.pValue(isX51);
        if hasModerator
            isXMo53 = strcmpi(lm5{i}.anova.Term,[treatment ':' moderator]) | strcmpi(lm5{i}.anova.Term,[moderator ':' treatment]);
            pXMo53(i) = lm5{i}.anova.pValue(isXMo53);
        end
    end
end


% Eq 6. Full model treatment -> outcome with mediators and optional moderator
lm6 = fitlme(T,eq6,'DummyVarCoding','effects');
pMe64 = nan(1,nrMediators);
pMeMo65 = nan(1,nrMediators*hasModerator);
for i=1:nrMediators
    if keepMediator(i)
        isMe64 = strcmpi(lm6.anova.Term,mediators{i});
        pMe64(i) = lm6.anova.pValue(isMe64);
        if hasModerator
            isMeMo65 = strcmpi(lm6.anova.Term,[mediators{i} ':' moderator]) | strcmpi(lm6.anova.Term,[moderator ':' mediators{i}]);
            pMeMo65(i) = lm6.anova.pValue(isMeMo65);
        end
    end
end


if ~strcmpi(dummyVar,'effects')
    % The request had a different coing style. Refit to get the appropriate
    % coefficients.
    lm4 = fitlme(T,eq4,'DummyVarCoding',dummyVar);
    for i=1:nrMediators
        if keepMediator(i)
            lm5{i} = fitlme(T,eq5{i},'DummyVarCoding',dummyVar);
        end
    end
    lm6 = fitlme(T,eq6,'DummyVarCoding',dummyVar);
end

%% Extract coefficients to determine paths
%% Calculate a,b,c,c'
% In the presence of moderators, we compute these at the lowest and highest
% categorical value of the moderator , or at +/- 1 stdev of the moderator
% for continuous moderators.
% X ---- c ---> Y
% X --- > a ---> Me
% Me ---> b ----> Y


isX41 =  startsWith(lm4.CoefficientNames,treatment) & ~contains(lm4.CoefficientNames,':');
bX41 = lm4.Coefficients.Estimate(isX41);
isX61 =  startsWith(lm6.CoefficientNames,treatment) & ~contains(lm6.CoefficientNames,':');
bX61 = lm6.Coefficients.Estimate(isX61);

if hasModerator
    if isa(T.(moderator),'double')
        sd = std(T.(moderator));
        m = mean(T.(moderator));
        moderatorLevels=m + sd.*[-1 1];
    else
        moderatorLevels= unique(T.(moderator));
    end
else
    moderatorLevels = 1;
end
nrModeratorLevels= numel(moderatorLevels);
a =nan(nrMediators,nrModeratorLevels);
b =nan(nrMediators,nrModeratorLevels);
bX51 = nan(nrMediators,1);
for i=1:nrMediators
    if keepMediator(i)
        isX51 =  startsWith(lm5{i}.CoefficientNames,treatment)  & ~contains(lm5{i}.CoefficientNames,':');
        bX51(i) = lm5{i}.Coefficients.Estimate(isX51);

        isMe64 =  startsWith(lm6.CoefficientNames,mediators{i}) & ~contains(lm6.CoefficientNames,':');
        bMe64 = lm6.Coefficients.Estimate(isMe64);

        if hasModerator
            % Have to pick two (Example) moderator values at which to evaluate the
            % mediation.
            if isa(T.(moderator),'double')
                % Continious variable as moderator: pick +/- 1 std of the moderator
                if i==1
                    isXMo43 =  strcmpi(lm4.CoefficientNames,[lm4.CoefficientNames{isX41} ':' moderator]) | strcmpi(lm4.CoefficientNames,[moderator ':' lm4.CoefficientNames{isX41} ]);
                    bXMo43 = lm4.Coefficients.Estimate(isXMo43)*moderatorLevels;
                    isXMo63 =  strcmpi(lm6.CoefficientNames,[lm6.CoefficientNames{isX61} ':' moderator]) | strcmpi(lm6.CoefficientNames,[moderator ':' lm6.CoefficientNames{isX61} ]);
                    bXMo63 =  lm6.Coefficients.Estimate(isXMo63)*moderatorLevels;
                end
                isXMo53 =  strcmpi(lm5{i}.CoefficientNames,[lm5{i}.CoefficientNames{isX51} ':' moderator]) |  strcmpi(lm5{i}.CoefficientNames,[moderator ':' lm5{i}.CoefficientNames{isX51} ]);
                bXMo53 =  lm5{i}.Coefficients.Estimate(isXMo53)*moderatorLevels;
                isMeMo65 =  strcmpi(lm6.CoefficientNames,[lm6.CoefficientNames{isMe64} ':' moderator]) | strcmpi(lm6.CoefficientNames,[moderator ':' lm6.CoefficientNames{isMe64}]);
                bMeMo65 =  lm6.Coefficients.Estimate(isMeMo65)*moderatorLevels;
            elseif isa(T.(moderator),'categorical') || isa(T.(moderator),'ordinal')
                if i==1
                    bXMo43 =moderatorBeta(lm4,treatment,moderator);
                    bXMo63 = moderatorBeta(lm6,treatment,moderator);
                end
                bXMo53 =  moderatorBeta(lm5{i},treatment,moderator);
                bMeMo65 =  moderatorBeta(lm6,mediators{i},moderator);
            end
        else
            bXMo43 = 0;
            bXMo63 = 0;
            bXMo53 = 0;
            bMeMo65=0;
        end
    else
        bXMo53 = NaN;
    end
    a(i,:) = bX51(i) + bXMo53;
    b(i,:) = bMe64 +bMeMo65;
end
c = bX41+bXMo43;
cPrime = bX61+bXMo63;
ab = a.*b;


%% With these results we can assess whether there is mediation/moderation


isMed =false(1,nrMediators);
isMedMod = false(1,nrMediators);
isModMed =false(1,nrMediators);
isModMedPrototypical =false(1,nrMediators);
if hasModerator
    %% Assess mediated moderation
    overallTreatmentModeration = pXMo43<alpha;
    for i=1:nrMediators
        if keepMediator
            isPattern1 = pXMo53(i) < alpha && pMe64(i) < alpha; %Mod affects treatment effect of mediator  && mediator affects outcome
            isPattern2 = pX51(i) < alpha && pMeMo65(i) < alpha; %Treatment affects mediator && moderator affects the mediator's effect on the outcome
            isMedMod(i) = overallTreatmentModeration && ( isPattern1 || isPattern2);

            %% Asses moderated mediation
            overalTreatmentEffect = pX41< alpha;
            isPrototypical = pXMo43 > alpha; % no interaction between treatment and mdoerator.
            isModMed(i)= overalTreatmentEffect && ( isPattern1 || isPattern2);
            isModMedPrototypical(i) = isModMed(i) && isPrototypical;
        end
    end
else
    %% Asses mediation
    overalTreatmentEffect = pX41< alpha;
    for i=1:nrMediators
        if keepMediator(i)
            treatmentAffectsMediator = pX51(i) < alpha;
            mediatorReducesTreatment  = pMe64(i) < alpha && abs(bX61) < abs(bX41); % The latter term needs statistical evaluation...
            isMed(i) = overalTreatmentEffect && treatmentAffectsMediator && mediatorReducesTreatment;
        end
    end
end

style.modMed = isModMed;
style.medMod = isMedMod;
style.med    = isMed;
style.prototypical = isModMedPrototypical;



%% Sensitivity
%  Imai et al correlation based sensitivity analysis.
% Not clear if this is valid for moderated mediation, multiple mediators,
% and mixed models....so commented out.
%
rho = nan(nrMediators,1);
for i=1:nrMediators
    if keepMediator(i)
        rho(i) = corr(lm5{i}.residuals,lm6.residuals);
    end
end
% uT = unique(T.(treatment));
% nrTreatment = numel(uT);
% e1 = lm4.residuals;
% for t=1:nrTreatment
%     stayT =T.(treatment)==uT(t);
%     sigma1t(t) =sqrt(var(e1(stayT)));
%     for i=1:nrMediators
%         e2 = lm5{i}.residuals;
%         sigma2t(i,t) = sqrt(var(e2(stayT)));
%         thisR = corrcoef(e1(stayT),e2(stayT));
%         rhot(i,t) = thisR(1,2);
%     end
% end
% rho = repmat(rho,[1 nrTreatment]);
% delta  = (repmat(bX51,[1 nrTreatment]).*repmat(sigma1t,[nrMediators 1])./sigma2t).* (rhot-rho.*sqrt((1-rhot.^2)./(1-rho.^2)))

%% Effect sizes
% Preacher KJ, Kelley K (2011) Effect size measures for mediation models:
% Quantitative strategies for communicating indirect effects. Psychol Methods 16:93–115.

% kappa2 - ab/max-possible(ab)
% This is really only defiend for 1 mediator and 0 moderators.
% This comutation simply ignores the other mediators (and the moderator) but it is not
% clear that this preserves the meanin of  'fraction of maximum possible mediation'
if hasModerator || nrMediators > 1
    kappa2 = nan(nrMediators,nrModeratorLevels);
else
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
    kappa2 = nan(nrMediators,hasModerator+1);
    for i=1:nrMediators
        if keepMediator(i)
            isMe = find(strcmpi(lm6.anova.Term,mediators{i}))+1;
            smx = S(isX,isMe);
            smm = S(isMe,isMe);
            sym = S(isY,isMe);
            tmp = sqrt(smm*syy-sym^2)*sqrt(sxx*syy-syx^2);
            aLim = [sym*syx-tmp, sym*syx+tmp]./(sxx*syy);
            bLim = [-1 1].* sqrt(sxx*syy-syx^2)/sqrt(sxx*smm-smx^2);
            extremeA = aLim(sign(aLim)==sign(a(i)));
            extremeB = bLim(sign(bLim)==sign(b(i)));
            kappa2(i) = a(i).*b(i)./(extremeA*extremeB);
        end
    end
end
%% Explained variance

%% Residuals



end
function [beta] = moderatorBeta(m, treatment,moderator)
isX =  startsWith(m.CoefficientNames,treatment) & ~contains(m.CoefficientNames,':');

% Find the XMo terms
if sum(isX)>1
    error('Not implemented yet; categorical treatment variable with more than 2 levels');
end

moderatorLevels = unique(m.Variables.(moderator));
nrLevels = numel(moderatorLevels);
beta= nan(1,nrLevels);
cntr=1;
coding = lm.dummyVarCoding(m);
for mo=moderatorLevels'
    % Find the first of the categorical/ordinal moderator levels that is in the model
    isXMo =  strcmpi(m.CoefficientNames,[m.CoefficientNames{isX} ':' moderator '_' char(mo)]) | strcmpi(m.CoefficientNames,[moderator '_' char(mo) ':' m.CoefficientNames{isX}]);
    if any(isXMo)
        beta(cntr) = m.Coefficients.Estimate(isXMo);
    else
        % term not in in the model
        switch upper(coding)
            case 'EFFECTS'
                isOtherXMo =  startsWith(m.CoefficientNames,[m.CoefficientNames{isX} ':' moderator '_' ]) | (startsWith(m.CoefficientNames,[moderator '_' ]) & endsWith(m.CoefficientNames, [':' m.CoefficientNames{isX}]));
                beta(cntr) = -sum(m.Coefficients.Estimate(isOtherXMo));
            case 'REFERENCE'
                beta(cntr) = 0;
            otherwise
                error(['Not implemented yet : ' coding])
        end
    end
    cntr= cntr+1;
end



end

