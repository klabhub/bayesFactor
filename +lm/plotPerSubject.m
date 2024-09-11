function [effects,ci] = plotPerSubject(m,pv)
% For a given linear model based on one or more subjects, refit the model
% per subject and show the results to get an idea of variability in the
% sample.
%
% INPUT
% glm -  a (generalized) linear mixed model
%
% OUTPUT
% effects  - table with fixed effects. First row is the group, the other
%               rows are the individual subjects
% ci - table with confidence intervals.
%
% BK - Feb 2020.
arguments
    m (1,1)  % Linear model
    pv.showHistogram (1,1) logical = false  % Histogram of effects across subjects
    pv.showLine (1,1)logical = true         % Line graph of effects for each subject
    pv.NUMPRECISION (1,1) double {mustBeInteger,mustBeNonnegative} =3; % num2str for effects
end


dummyVarCoding = lm.dummyVarCoding(m);

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
%
T = m.Variables(:,m.VariableInfo.InModel | ismember(m.VariableNames,m.Formula.ResponseName));

subjects = unique(T{:,m.Formula.GroupingVariableNames{:}}); % only keep subjects who're relevant to this condition

formula = m.Formula.char;
for s=subjects'
    try
        % Extract the relevant subset of data for this subjects
        thisT = T(~m.ObservationInfo.Excluded & T.subject==s,:);
        % Refit for this subject
        if isa(m,'LinearMixedModel')
            thisGlm = fitlme(thisT,formula,'FitMethod',m.FitMethod,'DummyVarCoding',dummyVarCoding) ;
        else
            lastwarn('')
            thisGlm = fitglme(thisT,formula,'FitMethod',m.FitMethod,'Distribution',...
                m.Distribution,'Link',m.Link,'DummyVarCoding',dummyVarCoding) ;
            [msg,id] = lastwarn;
            if strcmpi(id,'stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:Message_PLUnableToConverge')
                lastwarn('')
                error(msg);
            end
        end
        
        fe = thisGlm.fixedEffects;
        fe(1) =[]; % Remove intercept
        thisFeNames = thisGlm.CoefficientNames;
        thisFeNames(1) = [];

        effects = [effects   table(fe,'VariableNames',string(s))]; %#ok<AGROW>
        ci = [thisGlm.Coefficients.Lower thisGlm.Coefficients.Upper];
        ci(1,:) = [];
        cis  = [cis table(ci,'VariableNames',string(s))];%#ok<AGROW>

        %Sanity check that the order of FE is the same in the per subject and
        %group model. This can fail if one of the subjects does not
        % have a complete set of conditions
        assert(all(strcmpi(feNames,thisFeNames)),'coefficientNames of individual fit NOT matched with group fit! (check ''dummyVarCoding'' settings)');
    catch me
        fprintf('perSubject lmm for %s failed on %s (%s)\n',formula,s,me.message)
        continue
    end

end

%% Visualize the main group and individual results
nrSubjects = width(effects)-1;
nrEffects = height(effects);
nrRows = sum(pv.showHistogram+pv.showLine);

for e =1:nrEffects
    if pv.showHistogram
        subplot(nrRows,nrEffects,e);
        % Top row shows histogram of effects
        histogram(effects{e,2:end});


    end

    if pv.showLine

        subplot(nrRows,nrEffects,e+nrEffects*pv.showHistogram);
        % Bottom row shows line plots with CI per subject.
        line(reshape(cis{e,2:end},[2 nrSubjects]),repmat(1:nrSubjects,[2 1]),'Color','k')
        hold on
        plot(effects{e,2:end},1:nrSubjects,'k.')
        line([0 0],[0 nrSubjects])
        % Show the group effect as a red line
        line(cis{e,1},0.5+nrSubjects/2*[1 1],'Color','r','LineWidth',2);
        ylim([0 nrSubjects]);
        set(gca,'yTickLabels',{})
        if e==1
            ylabel 'Subject #'
            set(gca,'yTick',1:nrSubjects,'ytickLabel',subjects)
        end
    end

    if prod(sign(cis{e,1})) >0 % CI both <0 or both >0; significant at alpha level.
        sigStr = '(*)';
    else
        sigStr = '';
    end
    t=  title([feNames{e} ': ' num2str(effects{e,1},pv.NUMPRECISION) ' ' sigStr]);
    t.Interpreter = 'None';

    xlabel 'Effect'

end

