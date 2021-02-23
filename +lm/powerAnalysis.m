function [anovaPower,contrastPower,equivalencePower] = powerAnalysis(m,dummyVarCoding,varargin)
% Simulate a linear mixed model to generate a power analysis for factors in
% the model, specific posthoc contrasts, or tests of equivalence.
% 
% INPUT
% lm  =  A linear model with pilot data
% dummyVarCoding = The coding used for categorical variables in the model
%                   (Effects or Reference)
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
% NOTE
% With long/parfor evaluation, getting some updates on progress is useful. This
% script uses parfor_progress (See Matlab File exchange) for this if it is
% found on the path. If this funciton is not found, progress is indicated
% by dots printing to the screen.
%
% BK -Dec 2020

p=inputParser;
p.addRequired('lm',@(x) isa(x,'LinearMixedModel'))
p.addRequired('dummyVarCoding',@(x)(ischar(x) && ismember(upper(x),upper({'Effects','Reference'}))));
p.addParameter('subjectVariable','subject',@ischar);
p.addParameter('nrSubjects',10,@isnumeric);
p.addParameter('nrMonteCarlo',100,@isnumeric);
p.addParameter('alpha',0.05,@isnumeric); % Significance level
p.addParameter('anovaTerm',m.anova.Term,@iscell);
p.addParameter('contrast',[],@isnumeric); 
p.addParameter('equivalence',{},@iscell);
p.addParameter('nrWorkers',0,@isnumeric); % By default no parfor
p.addParameter('fixedEffectsScale',[],@isnumeric);
p.parse(m,dummyVarCoding,varargin{:});

if any(m.ObservationInfo.Excluded)
    error('This model has excluded some observations by using the ''Exclude'' argument to fitglme/fitlme. Please remove the data from the data table instead, then call fit without ''Exclude'' and then pass to this funciton');
end

if ~isempty(which('parfor_progress'))
    useParforProgress  = true; % Use parfor_pross
else
    useParforProgress  =false; % Show dots as signs of life.
end

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

clc;
startTime = tic;
if useParforProgress
    parfor_progress(nrN*nrMonteCarlo);
end
designMatrixRE = designMatrix(m,'Random');

parfor (n=1:nrN,nrWorkers)
    for i=1:nrMonteCarlo
         
        % Generate surrogate data based on the model
        subjectsToKeep = randi(nrSubjectsAvailable,[nrSubjectsToSimulate(n) 1]); % Sample subjects with replacement 
        nrSims = ceil(nrSubjectsToSimulate(n)/nrSubjectsAvailable);
        nrSubjectsSoFar = 0;
        simT= [];
        for s = 1:nrSims
            % For each s we generates the same number of responses as in the original data.
            % The builtin random function generates random effects as well
            % as error on each call. So these are "new" subjects. In the
            % manual calculation we do the same (while allowing some
            % scaling of fixed effects)
            if isempty(fixedEffectsScale)
                simResponse = random(m); % Use built-in - it handles all random effecs, including multiple grouping.
            else
                % We are scaling fixed effects; only implemented for a
                % single (subject) grouping va/random effect.
                simRandomEffects = m.covarianceParameters{1}*randn([size(designMatrixRE,2) 1]);
                simResponse = (designMatrix(m,'Fixed').*fixedEffectsScale)*m.fixedEffects + designMatrixRE*simRandomEffects + sqrt(m.MSE)*randn(height(m.Variables),1);             
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
            for sub = uSubjectID(thisSubjects)' %#ok<PFBNS>
                keep = ismember(m.Variables.(subjectVariable),sub);
                thisSim = m.Variables(keep,:);
                thisSim.(m.ResponseName) = simResponse(keep);
                simT = [simT;thisSim];  %#ok<AGROW>
            end
        end
        % Refit the model
        if isa(m,'GeneralizedMixedModel')
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
                contrastPValue(c,n,i)= coefTest(lmSim,contrast(c,:)); %#ok<PFBNS>
            end
        end
        
        
        % Evaluate equiavlence tests, if requested 
        for e = 1:nrEquivalenceTestsHack
            if nrEquivalenceTestsHack >nrEquivalenceTests
                break;
            else
                equivalencePValue(e,n,i) = lm.tost(lmSim,'Effects',equivalence{e,:}); %#ok<PFBNS>
            end
        end
                             
        if useParforProgress
            it = parfor_progress;
            showTimeToCompletion( it/100, [], [], startTime );
        else
            fprintf('.');
        end
    end
end
 
if ~useParforProgress
    fprintf('\n');
end

anovaPower = nanmean(anovaPValue<p.Results.alpha,3)';
contrastPower = nanmean(contrastPValue<p.Results.alpha,3)';
equivalencePower = nanmean(equivalencePValue<p.Results.alpha,3)';


end
