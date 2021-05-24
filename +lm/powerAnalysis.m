function [anovaPower,contrastPower,equivalencePower,hWaithBar] = powerAnalysis(m,varargin)
% Simulate a linear mixed model to generate a power analysis for factors in
% the model, specific posthoc contrasts, or tests of equivalence.
%
% INPUT
% lm  =  A linear model with pilot data
% Parm/Value pairs:
% subjectVariable  - the name of the variable in the lm that contains the
% subject  ID. [subject].
% nrSubjects: speciy the vector of subject sample sizes to investigate
% nrMonteCarlo -  MC simulations per sample size. [100]
% alpha - significance level
%
% Power Analysis will be done for the following model terms:
% anovaTerm - Cell array with anova terms
% contrast  - Matrix with specific model contrasts (see lm.posthoc or coefTest how to specify contrasts)
% equivalence - Cell array where each row defines an equivalence test
%               {A,B,bounds}; see lm.tost
% nrWorkers = The number of workers for parfor [0]
% OUTPUT
% Each row corresponds to a number of simulated subjects, each column to a
% term/contrast/equivalence test.
%
% anovaPower - Power for each factor specified in anovaTerm
% contrastPower - Power for each row specified in contrast
% equivalencePower - Power for reach row in equivalence
%
%
% BK -Dec 2020

p=inputParser;
p.addRequired('lm',@(x) (isa(x,'GeneralizedLinearMixedModel') || isa(x,'LinearMixedModel')))
p.addParameter('subjectVariable','subject',@ischar);
p.addParameter('nrSubjects',10,@isnumeric);
p.addParameter('nrMonteCarlo',100,@isnumeric);
p.addParameter('alpha',0.05,@isnumeric); % Significance level
p.addParameter('anovaTerm',m.anova.Term,@iscell);
p.addParameter('contrast',{},@iscell);
p.addParameter('equivalence',{},@iscell);
p.addParameter('nrWorkers',0,@isnumeric); % By default no parfor
p.addParameter('fixedEffectsScale',[],@isnumeric);
p.addParameter('graph',false,@islogical);
p.addParameter('names',{},@iscell);
p.parse(m,varargin{:});

if any(m.ObservationInfo.Excluded)
    error('This model has excluded some observations by using the ''Exclude'' argument to fitglme/fitlme. Please remove the data from the data table instead, then call fit without ''Exclude'' and then pass to this funciton');
end

dummyVarCoding = lm.dummyVarCoding(m);
% Extract from p to avoid broadcasting and initialize outputs
nrSubjectsToSimulate= p.Results.nrSubjects(:)';
nrWorkers = p.Results.nrWorkers;
nrMonteCarlo = p.Results.nrMonteCarlo;
subjectVariable = p.Results.subjectVariable;
contrast = p.Results.contrast;
equivalence  =  p.Results.equivalence;
nrEquivalenceTests= size(equivalence,1);
nrEquivalenceTestsHack = max(1,nrEquivalenceTests);
anovaTerm = p.Results.anovaTerm;
nrAnovaTerms = size(anovaTerm ,1);
nrAnovaTermsHack = max(nrAnovaTerms,1);
uSubjectID = unique(m.Variables.(subjectVariable));
nrSubjectsAvailable  = numel(uSubjectID);
nrN = numel(nrSubjectsToSimulate);
nrContrasts = size(contrast,1);
nrContrastsHack = max(nrContrasts,1);
anovaPValue= nan(nrAnovaTerms,nrN,p.Results.nrMonteCarlo);
contrastPValue= nan(nrContrasts,nrN,p.Results.nrMonteCarlo);
equivalencePValue = nan(nrEquivalenceTests,nrN,p.Results.nrMonteCarlo);
if isempty(p.Results.fixedEffectsScale)
    fixedEffectsScale = [];
elseif numel(p.Results.fixedEffectsScale)==numel(m.fixedEffects)
    % Same scaling for each observatin
    fixedEffectsScale  = repmat(p.Results.fixedEffectsScale,[m.NumObservations 1]);
elseif all(size(p.Results.fixedEffectsScale)==[m.NumObservations numel(m.fixedEffects)])
    % Each observation with its own scale
    fixedEffectsScale  =p.Results.fixedEffectsScale;
else
    error(['The size of the fixed effect scaling [' num2str(size(p.Results.fixedEffectsScale)) '] does not match the data [' num2str([m.NumObservations numel(m.fixedEffects)]) ']'])
end

if ~isempty(fixedEffectsScale)
    assert(numel(m.covarianceParameters)==1 && numel(m.covarianceParameters{1})==1,'Effect scaling only works for a single random effect/grouping variable');
end


% Use a dataqueue to show progress updates
dataQueue = parallel.pool.DataQueue;
hWaithBar = waitbar(0,['Power analysis for ' m.Formula.char ]);
cntr =0;
    function updateWaitBar(~)
        cntr=  cntr+1;
        waitbar(cntr/(nrN*nrMonteCarlo),hWaithBar);
    end
afterEach(dataQueue,@updateWaitBar);

parfor (n=1:nrN,nrWorkers)
    
    for i=1:nrMonteCarlo
        
        % Generate surrogate data based on the model
        subjectsToKeep = randi(nrSubjectsAvailable,[nrSubjectsToSimulate(n) 1]); % Sample subjects with replacement
        nrSims = ceil(nrSubjectsToSimulate(n)/nrSubjectsAvailable);
        nrSubjectsSoFar = 0;
        simT= [];
        subjectNr = 0;
        for s = 1:nrSims
            % For each s we generates the same number of responses as in the original data.
            % The builtin random function generates random effects as well
            % as error on each call. So these are "new" subjects. In the
            % manual calculation we do the same (while allowing some
            % scaling of fixed effects)
            if isempty(fixedEffectsScale)
                simResponse = random(m); % Use built-in - it handles all random effecs, including multiple grouping.
            else
                % We are scaling fixed effects
                if isa(m,'LinearMixedModel')
                    % For lmm random takes a modified deisng matrix
                    simResponse = random(m,designMatrix(m,'Fixed').*fixedEffectsScale,designMatrix(m,'Random'));
                else %Generalized LMM don't have this option.
                    % Hack to get access to the private slme member of m
                    %  if s==1
                    % Only have to do this once.
                    st= warning('query');
                    warning('off', 'MATLAB:structOnObject'); % Avoid the warning
                    modelStruct = struct(m);
                    warning(st);
                    % end
                    %Modify the design matrix with the fixed effect scaling
                    X = modelStruct.slme.X .* fixedEffectsScale;
                    % Then use code copied from
                    % GeneralizedLinearMixedModel/random
                    % To call random on the slme
                    wp      = modelStruct.slme.PriorWeights;
                    delta   = modelStruct.slme.Offset;
                    ntrials = modelStruct.slme.BinomialSize;
                    ysim    = random(modelStruct.slme,[],X,modelStruct.slme.Z,delta,wp,ntrials);
                    subset       = modelStruct.ObservationInfo.Subset;
                    simResponse= NaN(length(subset),1);
                    simResponse(subset) = ysim;
                end
            end
            % Use what we need, while making sure to sample a complete "data set" for each subject.
            if s==nrSims
                nrSubjectsNow = nrSubjectsToSimulate(n)- nrSubjectsSoFar;
            else
                % All
                nrSubjectsNow = nrSubjectsAvailable;
            end
            thisSubjects = subjectsToKeep(nrSubjectsSoFar+(1:nrSubjectsNow));
            nrSubjectsSoFar= nrSubjectsSoFar+nrSubjectsNow;
            for sub = uSubjectID(thisSubjects)'  %#ok<PFBNS>
                keep = ismember(m.Variables.(subjectVariable),sub);
                thisSim = m.Variables(keep,:);
                thisSim.(m.ResponseName) = simResponse(keep);
                % Assign a new ID to each subject
                thisSim.(subjectVariable) = repmat(categorical(subjectNr),[sum(keep) 1]);
                subjectNr = subjectNr +1;
                simT = [simT;thisSim];
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
                stay = strcmp(anovaTerm{a,:},thisAnova.Term); %#ok<PFBNS>
                anovaPValue(a,n,i) = thisAnova.pValue(stay);
            end
        end
        
        
        % Evaluate specific contrasts, if requested
        for c = 1:nrContrastsHack
            if nrContrastsHack >nrContrasts
                break;
            else
                contrastPValue(c,n,i)= lm.posthoc(lmSim,contrast{c,:}); %#ok<PFBNS>
            end
        end
        
        
        % Evaluate equiavlence tests, if requested
        for e = 1:nrEquivalenceTestsHack
            if nrEquivalenceTestsHack >nrEquivalenceTests
                break;
            else
                equivalencePValue(e,n,i) = lm.tost(lmSim,equivalence{e,:}); %#ok<PFBNS>
            end
        end
        
        % Update the data queue
        send(dataQueue,(n-1)+i);
    end
end

close(hWaithBar);

anovaPower = nanmean(anovaPValue<p.Results.alpha,3)';
contrastPower = nanmean(contrastPValue<p.Results.alpha,3)';
equivalencePower = nanmean(equivalencePValue<p.Results.alpha,3)';


h=[];

%% Show graphical results
if p.Results.graph    
    %Interpolate subjects for the graph
    iSubjects= (min(nrSubjectsToSimulate):1:max(nrSubjectsToSimulate))';
    powerValues = {anovaPower,contrastPower,equivalencePower};
    for i=1:3
        if ~isempty(powerValues{i})
            iPower = interp1(nrSubjectsToSimulate,powerValues{i},iSubjects,'pchip');
            
            % Binomial confidence intervals
            [x,ci] = binofit(round(powerValues{i}(:)*nrMonteCarlo),nrMonteCarlo);
            neg = reshape(x-ci(:,1),size(powerValues{i}));
            pos = reshape(ci(:,2)-x,size(powerValues{i}));
            
            hE = errorbar(repmat(nrSubjectsToSimulate',[1 size(powerValues{i},2)]),powerValues{i},neg,pos,[ 'o']);
            hold on
            for j=1:numel(hE)
                h = [h  plot(iSubjects,iPower(:,j),'LineWidth',2,'Color',hE(j).Color)];             %#ok<AGROW>
            end
            if ~isempty(p.Results.names )
                legend(h,p.Results.names)
            end
        end
        xlabel '#Subjects'
        ylabel 'Power'
        set(gca,'YLim',[0 1],'XTick',nrSubjectsToSimulate)
    end
    drawnow;
end
end
