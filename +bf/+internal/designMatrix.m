function [X,y] = designMatrix(lm,allTerms,varargin)
% Extract a dummy encoded design matix from a linear model .
% All variables are assumed to be categorical.
% INPUT
% lm =  A linear mixed efffects model
% allTerms = A cell array of variable names to include int he
% design arrrag ( e.g. {'ori','freq','ori:freq'} for two mains and an interaction)
% Parm/Value
% ZeroSumConstraint -  Toggle to apply the zero sum constraint to
%                   each factor (This equates the marginal prior across terms in the
%                           factor; see Rouder et al) [true]
% treatAsRandom  - A cell array of factors which should be
%               treated as random (i.e. not fixed) effects.
%               [{}].
% OUTPUT
% X = The design matrix [nrObservations nrLevels]
% y = The response data
% BK - 2019.

p = inputParser;
p.addParameter('zeroSumConstraint',true,@islogical);
p.addParameter('treatAsRandom',{},@iscell);
p.parse(varargin{:});

nrAllTerms =numel(allTerms);
X = cell(1,nrAllTerms);
isCategorical=  ~ismember(lm.Variables.Properties.VariableNames,lm.ResponseName);
opts = {'model','linear','categoricalvars',isCategorical,'intercept',false,'DummyVarCoding','full','responseVar',lm.ResponseName};
for i=1:nrAllTerms
    if any(allTerms{i}==':')
        % An interaction term.
        aName =extractBefore(allTerms{i},':');
        bName = extractAfter(allTerms{i},':');
        thisA = classreg.regr.modelutils.designmatrix(lm.Variables,'PredictorVars',aName,opts{:});
        thisB = classreg.regr.modelutils.designmatrix(lm.Variables,'PredictorVars',bName,opts{:});
        if ~ismember(aName,p.Results.treatAsRandom) && p.Results.zeroSumConstraint
            thisA = bf.internal.zeroSumConstraint(thisA);
        end
        if ~ismember(bName,p.Results.treatAsRandom) && p.Results.zeroSumConstraint
            thisB = bf.internal.zeroSumConstraint(thisB);
        end
        thisX = bf.internal.interaction(thisA,thisB);
    else
        % A main term
        thisX = classreg.regr.modelutils.designmatrix(lm.Variables,'PredictorVars',allTerms{i},opts{:});
        % Sum-to-zero contrasts that equates marginal priors across levels.
         if ~ismember(allTerms{i},p.Results.treatAsRandom) && p.Results.zeroSumConstraint           
            thisX = bf.internal.zeroSumConstraint(thisX);
        end
    end
    X{i} =thisX; % Store for later use
end
if nargout>1
    y= lm.Variables.(lm.ResponseName);
end
end
