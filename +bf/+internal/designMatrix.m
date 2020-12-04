function [X,y,isContinuous] = designMatrix(tbl,allTerms,responseVar,varargin)
% Extract a design matix from a linear model.
%  Dummy coded designs for categorical variables and zero-sum covariates for
%  contiunuous variables are returned as X
% INPUT
% tbl =  A table
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
%N = height(tbl);

%Options for modelutils.designmatrix . Its internal alg determines which
%vars are categorical quite well. But use categorical() in the data table
%to be sure (or chars/strings).
isCategorical = bf.internal.isCategorical(tbl);
if p.Results.forceCategorical
    isCategorical = true(size(isCategorical));
end
cateoricalOpts = {'model','linear','intercept',false,'DummyVarCoding','full','responseVar',responseVar,'CategoricalVars',isCategorical};
isContinuous = false(1,nrAllTerms);
for i=1:nrAllTerms
    if any(allTerms{i}==':')
        % An interaction term.
        thisTerms =strsplit(allTerms{i},':');
        thisNrTerms = numel(thisTerms);
        thisDm = cell(1,thisNrTerms);
        thisIsCategorical = cellfun(@(name) bf.internal.isCategorical(tbl,name),thisTerms);
        
        if all(thisIsCategorical)
            % Interaction between categorical variables
            for j=1:thisNrTerms
                % Extract design matrices
                thisDm{j} =  classreg.regr.modelutils.designmatrix(tbl,'PredictorVars',thisTerms{j},cateoricalOpts{:});
                if ~ismember(thisTerms{j},treatAsRandom) && p.Results.zeroSumConstraint
                    % Categorical variable with a zero-sum constraint
                    thisDm{j} = bf.internal.zeroSumConstraint(thisDm{j});
                end
            end
            thisX = bf.internal.interaction(thisDm{:});
        elseif all(~thisIsCategorical)
            % Interaction between continuous variables - extract from table
            thisX = prod(table2array(tbl.(thisTerms{:})),2);
            thisX = thisX - mean(thisX); % Remove mean
            isContinuous(i) = true;
        else % Mixture of categorical and continuous
            isContinuous(i) = true;
            % Prepare deisng matrix for the first term
            if thisIsCategorical(1)
                thisX = classreg.regr.modelutils.designmatrix(tbl,'PredictorVars',thisTerms{1},'intercept',false,'model','linear','responseVar',responseVar,'DummyVarCoding','full');
                others = 2:thisNrTerms;
                if any(ismember(thisTerms{others},allTerms))
                    % B is already included as a main effect, have to
                    % remove one level from A to avoid colinearity,
                    thisX = bf.internal.zeroSumConstraint(thisX);
                end
            else
                thisX = tbl.(thisTerms{1});
            end
            % Then multiply with successive terms
            for j=2:thisNrTerms
                if thisIsCategorical(j)
                    thisB = classreg.regr.modelutils.designmatrix(tbl,'PredictorVars',thisTerms{j},'intercept',false,'model','linear','responseVar',responseVar,'DummyVarCoding','full');
                    others = setdiff(1:thisNrTerms,j);
                    if any(ismember(thisTerms(others),allTerms))
                        % A is already included as a main effect, have to
                        % remove one level from B to avoid colinearity,
                        thisB = bf.internal.zeroSumConstraint(thisB);
                    end
                else
                    thisB = tbl.(thisTerms{j});
                end
                thisX  =  thisX.*thisB;
            end
            %Remove mean
            thisX  = thisX-mean(thisX);           
        end      
    else
        % A main term
        if  bf.internal.isCategorical(tbl,allTerms{i}) || p.Results.forceCategorical
            thisX = classreg.regr.modelutils.designmatrix(tbl,'PredictorVars',allTerms{i},cateoricalOpts{:});
            % Sum-to-zero contrasts that equates marginal priors across levels.
            if ~ismember(allTerms{i},treatAsRandom) && p.Results.zeroSumConstraint
                thisX = bf.internal.zeroSumConstraint(thisX);
            end
        else %Continuous
            thisX =  tbl.(allTerms{i});
            thisX  = thisX-mean(thisX);
            isContinuous(i)= true;
        end        
    end
    if isempty(thisX)
        error('The %s term does not vary in this table. Check your data table.', allTerms{i});
    end
    X{i} =thisX; % Store for later use
end

if nargout>1
    y= tbl.(responseVar);
end
end
