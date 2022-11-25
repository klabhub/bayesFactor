function [anovaPower,contrastPower,equivalencePower,deltaContrast,h,fixedEffects,randomEffects] = powerAnalysis(m,varargin)
% Simulate a linear mixed model to generate a power analysis for factors in
% the model, specific posthoc contrasts, or tests of equivalence.
%
% INPUT
% lm  =  A linear model with pilot data
% Parm/Value pairs:
% subjectVariable  - the name of the variable in the lm that contains the
% subject  ID. [subject].
% betweenSubjectVariable - The name of the variable that identifies the
% group for a between subject design. ['']. If this is specified,
% simulations will pick the same nrSubjects from each group identified by
% this variable. 
% nrSubjects: speciy the vector of subject sample sizes to investigate
% nrMonteCarlo -  MC simulations per sample size. [100]
% alpha - significance level
% qcFunction - This should be a function that takes a data table (of the
% same kind as used in the original fit(g)lme and returns a pruned data
% table. If specified, this function is applied to each simulated data set
% to ensure the same quality control as for the real experiment.
%
% fixedEffectsScale = Row vector (matching the size of fixedEffects) setting the factor with which each fixed
% effect is multiplied (e..g by 0 to simulate a null hypothesis).  Defaults
% to 1  (i.e. no scaling)
% randomEffectsScale = Row vector to scale random effects. Defaults to 1
% (no scaling)
% multipleComparisonsProcedure = 'None','FDR','Bonferroni'  - applied only to the
%                                                           (multiple) contrasts and equivalences
% Power Analysis will be done for the following model terms:
% anovaTerm - Cell array with anova terms
% contrast  - Cell array in which each row specifies a specific contrast to probe. (see lm.posthoc or
% coefTest how to specify contrasts). 
% equivalence - Cell array where each row defines an equivalence test
%               {A,B,bounds}; see lm.tost
% nrWorkers = The number of workers for parfor [0]
% graph - set to true to show the results in graphical form [false]
% names - Name each of the model terms (anova terms, contrasts,
% equivalences, in that order) - used only by the legend of the graph.
% 
% OUTPUT
% Each row corresponds to a number of simulated subjects, each column to a
% term/contrast/equivalence test.
%
% anovaPower - Power for each factor specified in anovaTerm
% contrastPower - Power for each row specified in contrast
% equivalencePower - Power for reach row in equivalence
% deltaContrast - Actual delta computed for each contrast .
% fixedEffects - [nrFixedEffects nrMonteCarlo] array of fixed effect
%                   estimates from each of the MC simulations
% ranodomEffects - [nrRandomEffects nrMonteCarlo] array of random effect
%                   estimates from each of the MC simulations
% BK -Dec 2020

p=inputParser;
p.addRequired('lm',@(x) (isa(x,'GeneralizedLinearMixedModel') || isa(x,'LinearMixedModel')))
p.addParameter('subjectVariable','subject',@ischar);
p.addParameter('betweenSubjectsVariable','',@ischar);
p.addParameter('nrSubjects',10,@isnumeric);
p.addParameter('nrMonteCarlo',100,@isnumeric);
p.addParameter('alpha',0.05,@isnumeric); % Significance level
p.addParameter('anovaTerm',m.anova.Term,@iscell);
p.addParameter('contrast',{},@iscell);
p.addParameter('equivalence',{},@iscell);
p.addParameter('nrWorkers',0,@isnumeric); % By default no parfor
p.addParameter('sesoi',[],@isnumeric);
p.addParameter('fixedEffectsScale',[],@isnumeric);
p.addParameter('randomEffectsScale',[],@isnumeric);
p.addParameter('qcFunction',@(x)(x),@(x)isa(x,'function_handle')); % Apply QC to each simulated set.
p.addParameter('maxLossToQc',1,@(x) (isnumeric(x) && x>0 && x <1)); % If QC removes more than this fraction [0 1], an error is gnerated.
p.addParameter('graph',false,@islogical);
p.addParameter('names',{},@iscell);
p.addParameter('multipleComparisonsProcedure','none',@(x)(ischar(x) && ismember(upper(x),{'FDR','BONFERRONI','NONE'})));
p.parse(m,varargin{:});

if any(m.ObservationInfo.Excluded)
    error('This model has excluded some observations by using the ''Exclude'' argument to fitglme/fitlme. Please remove the data from the data table instead, then call fit without ''Exclude'' and then pass to this funciton');
end

dummyVarCoding = lm.dummyVarCoding(m);
nrFixedEffects = size(m.fixedEffects,1);
nrRandomEffects =  size(m.randomEffects,1);
% Extract from p to avoid broadcasting and initialize outputs
alpha = p.Results.alpha;
maxLossToQc = p.Results.maxLossToQc;
nrSubjectsToSimulate= p.Results.nrSubjects(:)';
nrWorkers = p.Results.nrWorkers;
nrMonteCarlo = p.Results.nrMonteCarlo;
subjectVariable = p.Results.subjectVariable;
contrast = p.Results.contrast;
equivalence  =  p.Results.equivalence;
nrEquivalenceTests= size(equivalence,1);
nrEquivalenceTestsHack = max(1,nrEquivalenceTests);
anovaTerm = p.Results.anovaTerm;
mcp = p.Results.multipleComparisonsProcedure;
nrAnovaTerms = numel(anovaTerm);
nrAnovaTermsHack = max(nrAnovaTerms,1);
uSubjectID = unique(m.Variables.(subjectVariable));
nrSubjectsAvailable  = numel(uSubjectID);
nrN = numel(nrSubjectsToSimulate);
nrContrasts = size(contrast,1);
nrContrastsHack = max(nrContrasts,1);
anovaIsSig= nan(nrAnovaTerms,nrN,p.Results.nrMonteCarlo);
contrastIsSig= nan(nrContrasts,nrN,p.Results.nrMonteCarlo);
deltaContrast = nan(nrContrasts,nrN,p.Results.nrMonteCarlo);
equivalenceIsSig = nan(nrEquivalenceTests,nrN,p.Results.nrMonteCarlo);
fixedEffects = nan(nrFixedEffects,nrN,p.Results.nrMonteCarlo);
randomEffects = nan(nrRandomEffects,nrN,p.Results.nrMonteCarlo);
qcFunction  =p.Results.qcFunction;
nrRandomEffects = numel(m.randomEffects);
if ~isempty(p.Results.sesoi)
    % Use the sesoi to scale fixed effects (only)
    assert(size(contrast,1)==1 && numel(p.Results.sesoi) ==1,'Sesoi can only be computed for a single contrast at a time.');
    assert(isempty(p.Results.fixedEffectsScale) && isempty(p.Results.randomEffectsScale),'Choose either sesoiScale or fixed/random effects scale, but not both');
    fe = m.fixedEffects;    
    con = lm.contrast(m,contrast{1,1:2});
    modelContrast = con*fe;
    scale  =p.Results.sesoi./modelContrast;
    fixedEffectsScale  =ones(size(con));
    changeFe = con~=0;
    fixedEffectsScale(changeFe) = scale.*con(changeFe);
    fixedEffectsScale  = repmat(fixedEffectsScale,[m.NumObservations 1]);
    randomEffectsScale =  ones(m.NumObservations ,nrRandomEffects);
else
    % Setup scaling for fixed and random effects based on the user supplied
    % fixed/random effects scale
    fixedEffectsScale = p.Results.fixedEffectsScale;
    if isempty(fixedEffectsScale)
        fixedEffectsScale = ones(1,m.NumCoefficients);
    end
    assert(size(fixedEffectsScale,2) == m.NumCoefficients,'The numbers of elements in fixedEffectsScale should match the number of fixed effects coefficients in the model');
    fixedEffectsScale  = repmat(fixedEffectsScale,[m.NumObservations 1]);
    randomEffectsScale = p.Results.randomEffectsScale;    
    if isempty(randomEffectsScale)
        randomEffectsScale = ones(1,nrRandomEffects);
    end
    assert( size(randomEffectsScale,2) == nrRandomEffects,'The numbers of elements in randomEffectsScale should match the number of random effects coefficients in the model');
    randomEffectsScale  = repmat(randomEffectsScale,[m.NumObservations 1]);
end
% Use a dataqueue to show progress updates
dataQueue = parallel.pool.DataQueue;
hWaithBar = waitbar(0,['Power analysis for ' m.Formula.char ]);
cntr =0;
nrTotal = nrN*nrMonteCarlo;
    function updateWaitBar(~)
        cntr=  cntr+1;
        waitbar(cntr/nrTotal,hWaithBar);
    end
afterEach(dataQueue,@updateWaitBar);


% Hack to get access to the private slme member of m
st= warning('query');
warning('off', 'MATLAB:structOnObject'); % Avoid the warning
modelStruct = struct(m);
slme = modelStruct.slme;
warning(st);

%Modify the design matrix with the fixed effect scaling  (defaults to no
%change with fixedEffectsScale /randomEffectsScale empty
X = slme.X .* fixedEffectsScale; % Fixed effects
subset       = modelStruct.ObservationInfo.Subset;
Z = slme.Z.*randomEffectsScale; % Random effects
if isa(m,'LinearMixedModel')
    weights = modelStruct.ObservationInfo.Weights;
    weights = weights(subset);
    offset = NaN; % Not used but have to define parfor
    ntrials = NaN;
elseif isa(m,'GeneralizedLinearMixedModel')
    weights     = slme.PriorWeights;
    offset    = slme.Offset;
    ntrials = slme.BinomialSize;
else
    error('Unknown model type??');
end

mcpfun = @multicomp; % Hack to use the nested function in the parfor using feval

parfor (i=1:nrMonteCarlo ,nrWorkers ) % Parfor for the largest number
    %for i=1:nrMonteCarlo  % For debugggin without parfor
    for n=1:nrN        
        % Generate surrogate data based on the model
        if isempty(p.Results.betweenSubjectsVariable)
            % Within subjects design, just pick subjects from the pool
            subjectsToKeep = randi(nrSubjectsAvailable,[nrSubjectsToSimulate(n) 1]); %#ok<PFBNS> % Sample subjects with replacement                        
        else
            % A design with one between subjects variable; make sure to sample equally from each group (
            % nrSubjectsToSimulate is now interpreted as N per group.
            betweenSubjects = m.Variables.(p.Results.betweenSubjectsVariable);
            subjectsToKeep =[];
            for u=unique(betweenSubjects)'
                thisSubjects = unique(m.Variables.(subjectVariable)(ismember(betweenSubjects,u)));
                thisSubjectsToKeep = thisSubjects(randi(numel(thisSubjects),[nrSubjectsToSimulate(n) 1]));
                subjectsToKeep = [subjectsToKeep,thisSubjectsToKeep];                % Each column a group
            end
        end
        nrSubjectsSoFar = 0;
        simT= table; % Start with empty table
        subjectNr = 0;
        while nrSubjectsSoFar < nrSubjectsToSimulate(n)
            % For each s we generates the same number of responses as in the original data.
            % The builtin random function generates random effects as well
            % as error on each call. So these are "new" subjects.
            
            % Use code copied from
            % (Generalized)LinearMixedModel/random
            % To call random on the slme
            if isa(m,'LinearMixedModel')
                ysim = random(slme,[],X,Z);
                ysim = ysim ./ sqrt(weights);
            elseif isa(m,'GeneralizedLinearMixedModel')
                ysim    = random(slme,[],X,Z,offset,weights,ntrials);
            else 
                ysim = []; %#ok<NASGU> %Needed to avoid parfor warning
                error('Unknown model type');
            end
            simResponse= NaN(length(subset),1);
            simResponse(subset) = ysim;
            
            from = nrSubjectsSoFar+1;
            to = min(size(subjectsToKeep,1),nrSubjectsToSimulate(n));

            thisSubjects = subjectsToKeep(from:to,:);
            nrSubjectsSoFar= nrSubjectsSoFar+size(thisSubjects,1);
            for sub = uSubjectID(thisSubjects(:))'  %#ok<PFBNS>
                keep = ismember(m.Variables.(subjectVariable),sub);
                thisSim = m.Variables(keep,:);
                thisSim.(m.ResponseName) = simResponse(keep);
                % Assign a new ID to each subject
                thisSim.(subjectVariable) = repmat(categorical(subjectNr),[sum(keep) 1]);
                subjectNr = subjectNr +1;
                simT = [simT; thisSim]; %#ok<AGROW> 
            end
        end
        if ~isempty(qcFunction)
            % Apply quality control function to the simulated data set.
            before  = height(simT);
            simT = qcFunction(simT);
            after  = height(simT);
            loss = abs(after-before)/before;
            if loss > maxLossToQc
                fprintf('Losing too many observations to QC (%2.2f %%). Simulation skipped.\n', 100*loss);
                continue; % Go to the next
            end
        end
        
        % Refit the model
        if isa(m,'GeneralizedLinearMixedModel')
            lmSim =fitglme(simT,char(m.Formula),'Distribution',m.Distribution,'link',m.Link,'DummyVarCoding',dummyVarCoding);
        else
            lmSim =fitlme(simT,char(m.Formula),'DummyVarCoding',dummyVarCoding);
        end
        
        %Evaluate standard anova and store pValues
        for a = 1:nrAnovaTermsHack  % Hack is needed to trick the Matlab parfor parser in case nrAnovaTerms ==0
            if nrAnovaTermsHack >nrAnovaTerms
                break;
            else
                thisAnova = anova(lmSim);
                stay = strcmp(anovaTerm{a},thisAnova.Term); %#ok<PFBNS>
                anovaIsSig(a,n,i) = thisAnova.pValue(stay)<alpha;
            end
        end
        
        
        % Evaluate specific contrasts, if requested
        thisContrastP =nan(nrContrasts,1);
        for c = 1:nrContrastsHack
            if nrContrastsHack >nrContrasts
                break;
            else                
                [thisContrastP(c),~,~,deltaContrast(c,n,i)] = lm.posthoc(lmSim,contrast{c,:}); %#ok<PFBNS>
            end
        end
        contrastIsSig(:,n,i) = feval(mcpfun,thisContrastP,alpha,mcp); %#ok<FVAL>
        
        % Evaluate equiavlence tests, if requested
        thisEqP =nan(nrEquivalenceTests,1);
        for e = 1:nrEquivalenceTestsHack
            if nrEquivalenceTestsHack >nrEquivalenceTests
                break;
            else
                thisEqP(e) = lm.tost(lmSim,equivalence{e,:}); %#ok<PFBNS>                
            end
        end        
        equivalenceIsSig(:,n,i) = feval(mcpfun,thisEqP,alpha,mcp); %#ok<FVAL>
        
        % Store the fixed and random effects
        fixedEffects(:,n,i) = lmSim.fixedEffects;
        randomEffects(:,n,i) = lmSim.randomEffects;

        % Update the data queue
        send(dataQueue,(n-1)+i);
    end
end

close(hWaithBar);

anovaPower = mean(anovaIsSig,3,'omitnan')';
contrastPower = mean(contrastIsSig,3,'omitnan')';
equivalencePower = mean(equivalenceIsSig,3,'omitnan')';


h=[];

%% Show graphical results
if p.Results.graph
    %Interpolate subjects for the graph
    if nrN>1
        iSubjects= (min(nrSubjectsToSimulate):1:max(nrSubjectsToSimulate))';
    else
        iSubjects = nrSubjectsToSimulate;
    end
    powerValues = {anovaPower,contrastPower,equivalencePower};
    for i=1:3
        if ~isempty(powerValues{i}) && ~any(isnan(powerValues{i}),'all')
            if nrN>1
                iPower = interp1(nrSubjectsToSimulate,powerValues{i},iSubjects,'pchip');
            else
                iPower = powerValues{i};
            end
            
            % Binomial confidence intervals
            [x,ci] = binofit(round(powerValues{i}(:)*nrMonteCarlo),nrMonteCarlo);
            neg = reshape(x-ci(:,1),size(powerValues{i}));
            pos = reshape(ci(:,2)-x,size(powerValues{i}));
            
            hE = errorbar(repmat(nrSubjectsToSimulate',[1 size(powerValues{i},2)]),powerValues{i},neg,pos,'o');
            hold on
            for j=1:numel(hE)
                h = [h  plot(iSubjects,iPower(:,j),'LineWidth',2,'Color',hE(j).Color)];             %#ok<AGROW>
            end
        end
        xlabel (['# ' subjectVariable 's'])
        ylabel 'Power'
        set(gca,'YLim',[0 1],'XTick',nrSubjectsToSimulate)
    end
    if ~isempty(p.Results.names )
        legend(h,p.Results.names)
    end
    drawnow;
end


function isSig = multicomp(p,alpha,mcp)
% Multiple comparison procedure
nrComps = numel(p);
switch (upper(mcp))
    case 'NONE'
        isSig = p<alpha;
    case 'BONFERRONI'
        isSig = p<alpha./nrComps;
    case 'FDR'
        %Benjamini - Hochberg FDR.
        [sortedP,ix] = sort(p,'ascend');
        isSig = sortedP< (1:nrComps)*alpha/nrComps;
        [~,ix ] =sort(ix,'ascend');
        isSig = isSig(ix);
    otherwise
        error('Unknown MCP %s',mcp);
end
end
end

