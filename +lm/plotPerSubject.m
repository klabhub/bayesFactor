function [effects,ci] = plotPerSubject(m,dummyVarCoding)
% For a given linear model based on one or more subjects, refit the model
% per subject and show the results to get an idea of variability in the
% sample.
%
% INPUT
% glm -  a (generalized) linear mixed model
% dummyVarCoding - 'effects', 'reference', or 'full', use the same setting
% as one used for the group fitglme
%
% OUTPUT
% effects  - table with fixed effects. First row is the group, the other
%               rows are the individual subjects
% ci - table with confidence intervals.
%
% BK - Feb 2020.

NUMPRECISION =3; % num2str for effects


%% Createa a table of fixed effects and confidence intervals
fe = m.fixedEffects;
fe(1) =[]; % Remove intercept
feNames = m.CoefficientNames;
feNames(1) = [];
effects = table(fe,'RowNames',feNames,'VariableNames',{'Group'});
ci = [m.Coefficients.Lower m.Coefficients.Upper];
ci(1,:) = [];
cis  = table(ci,'RowNames',feNames,'VariableNames',{'Group'});

%% Now fit each subject separately and add to the tables
T = m.Variables;
subjects = unique(cellstr(T.subject(~m.ObservationInfo.Excluded))); % only keep subjects who're relevant to this condition
formula = m.Formula.char;
for s=subjects'
    try
        thisGlm = fitglme(T,formula,'FitMethod',m.FitMethod,'Distribution',...
            m.Distribution,'Link',m.Link,'DummyVarCoding',dummyVarCoding,...
            'Exclude',m.ObservationInfo.Excluded | ~strcmpi(cellstr(T.subject),s{1}));
        fe = thisGlm.fixedEffects;
        fe(1) =[]; % Remove intercept
        thisFeNames = thisGlm.CoefficientNames;
        thisFeNames(1) = [];
        
        effects = [effects   table(fe,'VariableNames',{['s' s{1}]})]; %#ok<AGROW>
        ci = [thisGlm.Coefficients.Lower thisGlm.Coefficients.Upper];
        ci(1,:) = [];
        cis  = [cis table(ci,'VariableNames',{['s' s{1}]})];%#ok<AGROW>
    catch
        disp(['perSubject lmm for ' formula ' failed on ' s{1}])
        continue
    end
    if ~all(strcmpi(feNames,thisFeNames)) % not all levels are matched with the group fit
        error('coefficientNames of individual fit NOT matched with group fit! (check ''dummyVarCoding'' settings)')
    end
end

%% Visualize the main group and individual results
nrSubjects = width(effects)-1;
nrEffects = height(effects);
for e =1:nrEffects
    
    subplot(2,nrEffects,e);
    % Top row shows histogram of effects
    histogram(effects{e,2:end});
    
    if prod(sign(cis{e,1})) >0 % CI both <0 or both >0; significant at alpha level.
        sigStr = '(*)';
    else
        sigStr = '';
    end
    t=  title([feNames{e} ': ' num2str(effects{e,1},NUMPRECISION) ' ' sigStr]);
    t.Interpreter = 'None';
    
    
    subplot(2,nrEffects,e+nrEffects);
    % Bottom row shows line plots with CI per subject.
    line(reshape(cis{e,2:end},[2 nrSubjects]),repmat(1:nrSubjects,[2 1]),'Color','k')
    hold on
    plot(effects{e,2:end},1:nrSubjects,'k.')
    line([0 0],[0 nrSubjects])
    % Show the group effect as a red line
    line(cis{e,1},0.5+nrSubjects/2*[1 1],'Color','r','LineWidth',2);
    xlabel 'Effect'
    ylim([0 nrSubjects]);
    set(gca,'yTickLabels',{})
    if e==1
        ylabel 'Subject #'
        set(gca,'yTick',1:nrSubjects,'ytickLabel',subjects)
    end
end

