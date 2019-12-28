function [X,y,isContinuous] = designMatrix(lm,allTerms,varargin)
% Extract a design matix from a linear model.
%  Dummy coded designs for categorical variables and zero-sum covariates for
%  contiunuous variables are returned as X
% INPUT
% lm =  A linear mixed efffects model
% allTerms = A cell array of variable names to include in the
% design matrix ( e.g. {'ori','freq','ori:freq'} for two mains and an interaction)
% Parm/Value
% ZeroSumConstraint -  Toggle to apply the zero sum constraint to
%                   each categorical factor (This equates the marginal prior across terms in the
%                           factor; see Rouder et al) [true]
% treatAsRandom  - A char or a cell array of chars with factors that should be
%               treated as random (i.e. not fixed) effects.
%               [{}].
% forceCategorical - Logical to force each term to be treated as a
% categorical variable (useful for grouping variables like subject ID which
% may look continuous...)
% OUTPUT
% X = The complete design matrix. Cell array with one element per term.
% y = The response data
% isContinuous = Logical indicating which columns are continuous co-variates.
%
% BK - 2019.

p = inputParser;
p.addParameter('zeroSumConstraint',true,@islogical);
p.addParameter('treatAsRandom',{},@(x) ischar(x) || iscell(x));
p.addParameter('forceCategorical',false,@islogical);
p.parse(varargin{:});
treatAsRandom =p.Results.treatAsRandom;
if ischar(treatAsRandom);treatAsRandom = {treatAsRandom};end

nrAllTerms =numel(allTerms);
X = cell(1,nrAllTerms);
N = height(lm.Variables);

%Options for modelutils.designmatrix . Its internal alg determines which
%vars are categorical quite well. But use categorical() in the data table
%to be sure (or chars/strings).
isCategorical = bf.internal.isCategorical(lm.Variables);
if p.Results.forceCategorical
    isCategorical = true(size(isCategorical));
end
cateoricalOpts = {'model','linear','intercept',false,'DummyVarCoding','full','responseVar',lm.ResponseName,'CategoricalVars',isCategorical};
isContinuous = false(1,nrAllTerms);
for i=1:nrAllTerms
    if any(allTerms{i}==':')
        % An interaction term.
        aName =extractBefore(allTerms{i},':');
        bName = extractAfter(allTerms{i},':');
        aCategorical  = bf.internal.isCategorical(lm.Variables,aName);
        bCategorical = bf.internal.isCategorical(lm.Variables,bName);  
        bothCategorical = aCategorical &&  bCategorical;
        bothContinuous = ~aCategorical && ~bCategorical;
        if bothCategorical || p.Results.forceCategorical
            thisA = classreg.regr.modelutils.designmatrix(lm.Variables,'PredictorVars',aName,cateoricalOpts{:});
            if ~ismember(aName,treatAsRandom) && p.Results.zeroSumConstraint 
                thisA = bf.internal.zeroSumConstraint(thisA);
            end
            thisB = classreg.regr.modelutils.designmatrix(lm.Variables,'PredictorVars',bName,cateoricalOpts{:});
            if ~ismember(bName,treatAsRandom) && p.Results.zeroSumConstraint
                thisB = bf.internal.zeroSumConstraint(thisB);
            end
            thisX = bf.internal.interaction(thisA,thisB); 
        elseif bothContinuous
            thisX =  lm.Variables.(aName).*lm.Variables.(bName);
            thisX = thisX - mean(thisX);
            isContinuous(i) = true;
        else
            % Interaction between a categorical and a continuous covariate
            if aCategorical
                thisA = classreg.regr.modelutils.designmatrix(lm.Variables,'PredictorVars',aName,'intercept',false,'model','linear','responseVar',lm.ResponseName,'DummyVarCoding','full');
                thisA = thisA(:,2:end); % Remove first category 
            else
                thisA = lm.Variables.(aName);
            end
            if bCategorical
                thisB = classreg.regr.modelutils.designmatrix(lm.Variables,'PredictorVars',bName,'intercept',false,'model','linear','responseVar',lm.ResponseName,'DummyVarCoding','full');
                thisB = thisB(:,2:end); % Remove first category
            else
                thisB = lm.Variables.(bName);
            end            
            thisX  =  thisA.*thisB;
            thisX  = thisX-mean(thisX);
            isContinuous(i) = true;
        end        
    else
        % A main term
        if  bf.internal.isCategorical(lm.Variables,allTerms{i}) || p.Results.forceCategorical        
            thisX = classreg.regr.modelutils.designmatrix(lm.Variables,'PredictorVars',allTerms{i},cateoricalOpts{:});
            % Sum-to-zero contrasts that equates marginal priors across levels.
            if ~ismember(allTerms{i},treatAsRandom) && p.Results.zeroSumConstraint 
                thisX = bf.internal.zeroSumConstraint(thisX);
            end
        else %Continuous
             thisX =  lm.Variables.(allTerms{i});
             isContinuous(i)= true;
        end
        
    end
    X{i} =thisX; % Store for later use
end

if nargout>1
    y= lm.Variables.(lm.ResponseName);
end
end
