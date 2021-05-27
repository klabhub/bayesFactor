function [effects,ci] = plotPerSubject(m)
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

NUMPRECISION =3; % num2str for effects

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
T.subject = categorical(T.subject);
subjects = unique(T.subject); % only keep subjects who're relevant to this condition

formula = m.Formula.char;
stay  =m.VariableInfo.InModel & ~ismember(m.VariableNames,[m.Formula.GroupingVariableNames{:}]);
varTypes = m.VariableInfo.Class(stay);
varNames = m.VariableNames(stay);
switch upper(dummyVarCoding)
    case {'REFERENCE','EFFECTSFIRST'}
        matchIx = height(T);
    case {'REFERENCELAST','EFFECTS'}
        matchIx = 1;
    otherwise 
        warning('Not tested for  %s coding',dummyVarCoding);
end
rowToMatch = T(matchIx,:);
for s=subjects'
    try
        % Extract the relevant subset of data for this subjects
        thisT = T(~m.ObservationInfo.Excluded & T.subject==s,:);        
        % We have to order this such that the same variable serves as the
        % reference/left-out parameter.
        for r=1:numel(varNames)
            type = varTypes{r};
            switch(type)
                case 'cell'
                    comp = @ismember;
                    otherwise
                    comp = @eq;
            end
            thisStay = comp(thisT.(varNames{r}),rowToMatch.(varNames{r}));
            if r==1 
                stay = thisStay;
            else
                stay = stay & thisStay;
            end
        end        
        ix  = find(stay,1,'first');
        if isempty(ix)
            error('Subject %s does not have the condition that is left out of the model in the population fit.Skipped',s);
        end
        if matchIx==1
            %First value serves as reference : Prepend
            thisT = cat(1,thisT(ix,:),thisT(setdiff(1:height(thisT),ix),:));
        else
            % Last value serves as referrence : Append
            thisT = cat(1,thisT(setdiff(1:height(thisT),ix),:),thisT(ix,:));
        end
        % Refit for this subject
        thisGlm = fitglme(thisT,formula,'FitMethod',m.FitMethod,'Distribution',...
            m.Distribution,'Link',m.Link,'DummyVarCoding',dummyVarCoding) ;
        fe = thisGlm.fixedEffects;
        fe(1) =[]; % Remove intercept
        thisFeNames = thisGlm.CoefficientNames;
        thisFeNames(1) = [];
        
        effects = [effects   table(fe,'VariableNames',{['s' char(s)]})]; %#ok<AGROW>
        ci = [thisGlm.Coefficients.Lower thisGlm.Coefficients.Upper];
        ci(1,:) = [];
        cis  = [cis table(ci,'VariableNames',{['s' char(s)]})];%#ok<AGROW>
    catch
        disp(['perSubject lmm for ' formula ' failed on :'])
        s
        continue
    end
    %Sanity check that the order of FE is the same in the per subject and
    %group model
    assert(all(strcmpi(feNames,thisFeNames)),'coefficientNames of individual fit NOT matched with group fit! (check ''dummyVarCoding'' settings)');
     
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

