function [p,stat,df,delta,contrast,debugStr] = posthoc(m,dummyVarCoding,A,B,predictedDelta,tail)
% Perform a posthoc comparison on a Linear Mixed Model, using coefTest, but
% using words to specify the contrast.
%
% Not tested extensively and this is tricky code, so make sure to check the
% contrast and debugStr amd confirm that the contrast does what you think it
% does.
%
% INPUT
% lm                = The linear model
% dummyVarCoding    = The coding used for categorical variables.
% A                 =  Cell array specifying condition A
% B                 = Cell array specifying condition B
% predictedDelta    = Predicted difference between A and B. Defaults to 0.
% tail              = Specify the tail of the distribution. 'left','right' (use one sided T-tests)
%                       or 'both' (use F-test).
%
% OUTPUT
% p                 = The p-value associated with the test. (see coefTest)
% stat              =  F for two-tailed tests,  T for one sided
% df                = Error degrees of freedom (see coefTest)
% delta             = The difference.
% contrast          = The contrast used for this test.
% debugStr          = Contrast shown together with coefficient names to
%                       help understand why the contrast is the way it is...
% EXAMPLE
% Fit a model with two fixed-effect predictors and a random
%             effect. Test for the significance of the Cylinders term. The
%             p-value is the same as shown in the anova table.
% load carsmall
% T = table(MPG,Weight,Model_Year,Cylinders);
% T.Cylinders = nominal(T.Cylinders);
% glme = fitglme(T,'MPG ~ Weight + Cylinders + (1|Model_Year)','Distribution','Normal','DummyVarCoding','Effects')
% H0:  Cylinders have no effect on MPG (should match anova result)
% p = lm.posthoc(glme,'Effects','Cylinders')
% Or a pairwise: (H0 8 and 6 cylinders are the same.
% p = lm.posthoc(glme,'Effects','Cylinders_8','Cylinders_6')
%
% EXAMPLE
% Use cell arrays to select specific cells in a multifactorial design. E.g. in a design
% with stim and valid factors, each with multiple levels
% {'stim_VIS','valid_1'} represents the cell in which the stim factor has level VIS and the
% valid factor has level 1.
%
% {'stim'} -> Selects all levels of the stim factor. (With this, the test
% results is equivalent to the main effect in an ANOVA).
%
% BK -  Jan 2021

if nargin <6
    tail = 'both';
end
if nargin <5
    predictedDelta = 0;
end
if nargin<2 || ~ismember(upper(dummyVarCoding),upper({'Effects','Reference'}))
    error('DummyVarcoding must be speciied as either Effects or Reference');
end
if ischar(A) ;A ={A};end
vectorA = name2vector(m,A,dummyVarCoding);

if nargin>3 && ~isempty(B)
    if ischar(B) ;B ={B};end
    vectorB = name2vector(m,B,dummyVarCoding);
else
    vectorB = zeros(size(vectorA));
end

contrast = vectorA - vectorB;
delta = contrast*m.fixedEffects;
debugStr = strcat(m.CoefficientNames', ' : ' , cellstr(num2str(contrast')));

predictedDelta = predictedDelta(:).*ones(size(delta));
if ~all(size(delta)==size(predictedDelta))
    error('Predicted delta [nrRows nrCols] has to match the expected delta');
end

switch upper(tail)
    case 'BOTH'
        % Two-sided tests - use the built-in F test.
        [p,stat,~,df] = coefTest(m,contrast,predictedDelta);
    case {'LEFT','RIGHT'}
        cov  = contrast*m.CoefficientCovariance*contrast';
        stat = (delta-predictedDelta)/sqrt(cov);
        df  =  m.DFE;
        if strcmpi(tail,'LEFT')
            p = tcdf(stat,m.DFE);
        else
            p = 1-tcdf(stat,m.DFE);
        end
    otherwise
        error('%s is not a valid tail specification',tail);
end


end

function conditionVector = name2vector(m,X,dummyVarCoding)
%% determine the vector that defines coniditon X
% X = {'A_1','B_2'} means level 1 from factor A *and* level 2 from factor B
% so it will identify terms in the model that have A_1,  B_1, or
% A_1:B_2:xxx.
% For levels that are explicitly in the model this is relatively easy to
% do, but for the left-out level we need to find all other levels and then
% use the zero sum constraint to determine the vector.
%

conditionVector =zeros(1,m.NumCoefficients);
interactionTerms = {};
if numel(X)>2;error('I doubt this code is correct for more than 2 conditions');end
isMissingLevel = false(1,numel(X));
for i=1:numel(X)
    [mainTerms,interactionTerms, isMissingLevel(i)] = getTerms(m,X{i},interactionTerms);
    conditionVector = conditionVector+ terms2conditionVector(m,mainTerms,isMissingLevel(i),dummyVarCoding);    
    if i>1
        % Interactions with the earlier terms. OK for terms in the model,
        % but the "missing" terms are tricky. This wrks for 2 lvels, but
        % for higher order interactions something more complicated will be
        % needed?
        conditionVector = conditionVector+ terms2conditionVector(m,interactionTerms,any(isMissingLevel),dummyVarCoding);   %
    end
end
end
%%
function conditionVector = terms2conditionVector(m,terms,isMissingLevel,dummyVarCoding)
if isMissingLevel
    conditionVector  = zeros(1,m.NumCoefficients);
    if strcmpi(dummyVarCoding,'Effects')
        for j=1:numel(terms)
            conditionVector(1,:) =  conditionVector(1,:) - strcmp(m.CoefficientNames,terms{j});
        end
    else
        % It was the reference - keep zeros
    end
else
    conditionVector  = zeros(numel(terms),m.NumCoefficients);
    for j=1:numel(terms)
        conditionVector(j,:) =  conditionVector(j,:) + strcmp(m.CoefficientNames,terms{j});
    end
end


end
%%
function [mainTerms,interactionTerms, zerosum] = getTerms(m,x,allowedInteractionTerms)
if ~contains(x,'_')
    % Not a condition, possibly a main effect
    if ismember(x,m.VariableNames)
        zerosum = false;
        ix  =contains(m.CoefficientNames,[x '_']) & ~contains(m.CoefficientNames,':');
    else
        error('Posthoc comparison requires conditions |(with a _ in the name), or main effects, not interactions (%s)',x);
    end
else
    % A specific condition (factor_level)
    ix = find(~cellfun(@isempty,regexp(m.CoefficientNames,[x '\>'])));    
    if ~any(ix)
        % Level not in the model.
        % Make sure this is not a typo.
        factorName =extractBefore(x,'_');
        levelName = extractAfter(x,'_');
        if ~( ismember(factorName,m.VariableNames) && ismember(levelName,m.Variables.(factorName)))
            error('%s : not found in the variables of this model',x);
        end
        % It is the missing level; compute using zerosm
        zerosum = true;
        ix  =contains(m.CoefficientNames,[extractBefore(x,'_') '_']);
    else
        zerosum = false;
    end
end
allTerms  = m.CoefficientNames(ix);
interactionTerms = allTerms(contains(allTerms,':'));
mainTerms  = setdiff(allTerms, interactionTerms);
if ~isempty(allowedInteractionTerms)
    interactionTerms = intersect(allowedInteractionTerms,interactionTerms);
end
end