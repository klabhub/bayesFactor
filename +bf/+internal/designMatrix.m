function [X,y,contScaleFactor] = designMatrix(lm,allTerms,varargin)
% Extract a design matix from a linear model.
%  Dummy coded designs for categorical variables and zero-sum designs for
%  contiunuous variables are returned as X
% INPUT
% lm =  A linear mixed efffects model
% allTerms = A cell array of variable names to include int he
% design arrrag ( e.g. {'ori','freq','ori:freq'} for two mains and an interaction)
% Parm/Value
% ZeroSumConstraint -  Toggle to apply the zero sum constraint to
%                   each factor (This equates the marginal prior across terms in the
%                           factor; see Rouder et al) [true]
% treatAsRandom  - A char or a cell array of chars with factors that should be
%               treated as random (i.e. not fixed) effects.
%               [{}].
% OUTPUT
% X = The design matrix for categorical variables. Cell array with one element per term.
% y = The response data
% contScaleFactor = Zellner-Siow Prior Scale Factors for continuous
% variables, NaN for Categorical vars
%
% BK - 2019.

p = inputParser;
p.addParameter('zeroSumConstraint',true,@islogical);
p.addParameter('treatAsRandom',{},@(x) ischar(x) || iscell(x));
p.parse(varargin{:});
treatAsRandom =p.Results.treatAsRandom;
if ischar(treatAsRandom);treatAsRandom = {treatAsRandom};end

nrAllTerms =numel(allTerms);
X = cell(1,nrAllTerms);
%Options for modelutils.designmatrix . Its internal alg determines which
%vars are categorical quite well. But use categorical() in the data table
%to be sure (or chars/strings).
cateoricalOpts = {'model','linear','intercept',false,'DummyVarCoding','full','responseVar',lm.ResponseName};
isCategoricalColumn  = @(x) (isa(x,'categorical') || iscellstr(x) || isstring(x) || ischar(x) || islogical(x));
isCategoricalName = @(x) isCategoricalColumn(lm.Variables.(x));
contScaleFactor = nan(1,nrAllTerms);
for i=1:nrAllTerms
    if any(allTerms{i}==':')
        % An interaction term.
        aName =extractBefore(allTerms{i},':');
        bName = extractAfter(allTerms{i},':');
        bothCategorical = isCategoricalName(aName)  && isCategoricalName(bName);
        bothContinuous = ~isCategoricalName(aName) && ~isCategoricalName(bName);
        if bothCategorical 
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
            thisX =  classreg.regr.modelutils.designmatrix(lm.Variables,'PredictorVars',{aName,bName},'intercept',false,'model','interactions','responseVar',lm.ResponseName);
            N=size(thisX,1);
            thisX = thisX-sum(thisX)/N; % Sum =0;Rouder et al page 269.             
            thisScale = inv(thiX'*thisX/N);
            contScaleFactor(i) =thisScale(1,2); % Interaction term is the 1,2 element of the 2x2 matrix 
        else
            error('An interaction betwen categorical and continuous variables has not been implemented yet');            
        end        
    else
        % A main term
        if isCategoricalName(allTerms{i})
            thisX = classreg.regr.modelutils.designmatrix(lm.Variables,'PredictorVars',allTerms{i},cateoricalOpts{:});
            % Sum-to-zero contrasts that equates marginal priors across levels.
            if ~ismember(allTerms{i},treatAsRandom) && p.Results.zeroSumConstraint 
                thisX = bf.internal.zeroSumConstraint(thisX);
            end
        else %Continuous
             thisX =  classreg.regr.modelutils.designmatrix(lm.Variables,'PredictorVars',allTerms{i},'intercept',false,'model','linear','responseVar',lm.ResponseName);
             N=size(thisX,1);
             thisX = thisX-sum(thisX)/N; % Sum =0;Rouder et al page 269.
             contScaleFactor(i) = inv(thisX'*thisX/N); % For a single main effect, this is a scalar.
        end
    end
    X{i} =thisX; % Store for later use
end
if nargout>1
    y= lm.Variables.(lm.ResponseName);
end
end
