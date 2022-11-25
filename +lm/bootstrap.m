function [fixedEffects,randomEffects,loglike] = bootstrap(m,varargin)
% Bootstrap the fixed and random effect estimates of a linear mixed model
% by resampling (with replacement) at the subject level.
% In each bootstrap set, we take the data from a random subset of the subjects 
% to create a resampled data set with the same number of data rows and
% subjects (i.e. with replacement). Each of these subjects gets a new
% unique name, then the LME is fit using the same settings as the original
% model (m). For each boostrap set, the fixed effects, random effects, and
% the loglikelihood of the model is kept and returned by this function. 
% With 'graph' set to true, the function shows the bootstrapped fixed
% effect distributions.
%
% INPUT
% m - The  (generalized) linear mixed model
% 'subjectVariable' - Specify which variable should be resampled. ['subject']
% 'nrBootstrap' - Number of sets to compute [100]
% 'nrWorkers' - How many parallel workers to use [0]
% 'graph' - Show an graphical output [false]
%
% OUTPUT
% fixedEffects  - [nrFixedEffects nrBoostrap] Matrix with the fixed effects for each of the resampled sets.
% randomEffects - [nrRandomEffects nrBoostrap] Matrix with the  random effects for each of the resampled sets
% loglike  - The log likelihood for each boostrap set.
% BK -  Nov 2022.

p=inputParser;
p.addRequired('lm',@(x) (isa(x,'GeneralizedLinearMixedModel') || isa(x,'LinearMixedModel')))
p.addParameter('subjectVariable','subject',@ischar);
p.addParameter('nrMonteCarlo',100,@isnumeric);
p.addParameter('nrWorkers',0,@isnumeric); % By default no parfor
p.addParameter('graph',false,@islogical); % Show graphs
p.parse(m,varargin{:});

if any(m.ObservationInfo.Excluded)
    error('This model has excluded some observations by using the ''Exclude'' argument to fitglme/fitlme. Please remove the data from the data table instead, then call fit without ''Exclude'' and then pass to this function');
end

% Create local variables to reduce parfor broadcasting
nrMonteCarlo    = p.Results.nrMonteCarlo;
nrWorkers       = p.Results.nrWorkers;
subjectVariable = p.Results.subjectVariable;
T               = m.Variables;
subjects        = T.(subjectVariable);
uSubjects       = unique(subjects);
nrSubjects      = numel(uSubjects);
dummyVarCoding = lm.dummyVarCoding(m);
%#ok<*PFBNS>   % Suppress broadcasting message.


% Preallocate
nrFixedEffects  = size(m.fixedEffects,1);
nrRandomEffects =  size(m.randomEffects,1);
fixedEffects    = nan(nrFixedEffects,nrMonteCarlo);
randomEffects   = nan(nrRandomEffects,nrMonteCarlo);
loglike         = nan(1,nrMonteCarlo); % Log likelihood.
% Use a dataqueue to show progress updates
dataQueue = parallel.pool.DataQueue;
hWaithBar = waitbar(0,['Bootstrapping ' m.Formula.char ]);
cntr =0;
nrTotal = nrMonteCarlo;
    function updateWaitBar(~)
        cntr=  cntr+1;
        waitbar(cntr/nrTotal,hWaithBar);
    end
afterEach(dataQueue,@updateWaitBar);

parfor (i=1:nrMonteCarlo ,nrWorkers )
    %for i=1:nrMonteCarlo % Use for debugging
    % Select random subset of subjects with resampling and same total
    % number.
    subjectsToUse = uSubjects(randi(nrSubjects,[nrSubjects 1]));
    subjectsSoFar = 0;
    % Now create a new data table from the data for this subject, but
    % assign a new ID to each of the resampled subjects
    resampledT = table;
    for sub = subjectsToUse
        subjectsSoFar = subjectsSoFar +1; % Used to generate new IDs
        keep = ismember(subjects,sub);
        resampledSubject= m.Variables(keep,:);
        % Assign a new ID to each subject
        resampledSubject.(subjectVariable) = repmat(categorical(subjectsSoFar),[sum(keep) 1]);
        resampledT = [resampledT; resampledSubject];
    end
    % Estimate model parameters/
    if isa(m,'GeneralizedLinearMixedModel')
        thisM =fitglme(resampledT,char(m.Formula),'Distribution',m.Distribution,'link',m.Link,'DummyVarCoding',dummyVarCoding);
    else
        thisM =fitlme(resampledT,char(m.Formula),'DummyVarCoding',dummyVarCoding);
    end

    % Store the fixed and random effects
    fixedEffects(:,i) = thisM.fixedEffects;
    randomEffects(:,i) = thisM.randomEffects;
    loglike(i) = thisM.ModelCriterion.LogLikelihood;
    % Update the data queue
    send(dataQueue,i);
end
close(hWaithBar);


%% Graphical output
% Plot a histogram of each of the fixed effects plus the log likelihood
FE = m.fixedEffects;
if p.Results.graph
    clf;
    tiledlayout("flow");
    for f=1:nrFixedEffects
        nexttile;
        histogram(fixedEffects(f,:),'Normalization','pdf');
        hold on
        plot(FE(f)*[1 1],ylim,'k','LineWidth',2);
        title (sprintf('%s: %.2f CI [%.2f %.2f]',m.CoefficientNames{f},mean(fixedEffects(f,:),"omitnan"),prctile(fixedEffects(f,:),2.5),prctile(fixedEffects(f,:),97.5)) ,"Interpreter","none");
        xlabel 'Coefficient'
        ylabel 'PDF'
    end
    nexttile   
    histogram(loglike,'Normalization','pdf')
    title (sprintf('%s: %.2f CI [%.2f %.2f]','Log Likelihood:',mean(loglike,"omitnan"),prctile(loglike,2.5),prctile(loglike,97.5)) ,"Interpreter","none");
    xlabel 'Random Effect'
    ylabel 'PDF'    
end

end

