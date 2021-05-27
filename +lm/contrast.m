function [v,TA,TB] = contrast(m,A,B,defineDifference,scale)
% Specify two conditions A and B using cell arrays of parm/value pairs to
% retreive the contrast weights vector A-B.
%
% INPUT
% m - The linear model
% A - Condition A - A cell array of variable/value pairs.
% B - Condition B - B cell array of variable/value pairs.Can be empty, in
%       which case the function returns the vector defining A.
% defineDifference = When set tot true, the B cell array specifies only those
%           variables that are different in B. [true]
% scale -  Set this to true to scale the weights such that sum(abs(v)) ==2.
% This puts the contrast score on the same scale as the original means
% (Kirk, 1995). [false]
% 
% OUTPUT
% v = The contrast correspoonding to A-B
%
% BK - Mar 2021
if nargin<5
    scale = false;
    if nargin<4
        defineDifference =true;
    end
end

import lm.*
if isa(A,'table')
    % Use the conditions specified in the table.
    TA =A;
    defaultA = {};
else
    %% Create tables that define the conditions specified by A and (if requested ) B.
    varTypes  = m.VariableInfo.Class(m.VariableInfo.InModel);
    varNames = m.VariableNames(m.VariableInfo.InModel);
    defaultValues = m.Variables(1,m.VariableInfo.InModel);% First row in the data is the default
    TA= defaultValues;
    defaultA = setdiff(varNames,A(1:2:end))';
    TA = fillTable(varTypes,varNames,TA,A); % Replace default with values specified in A
end

% Same for condition B
if nargin >2 && ~isempty(B)
    if isa(B,'table')
        TB = B;
        defaultB = {};
    else
        if defineDifference
            TB = TA;  % Start with TA, replace what is in B
            defaultB = setdiff(varNames,cat(2,A(1:2:end),B(1:2:end)))';
            for i=1:2:numel(B)
                thisType = varTypes{strcmpi(B{i},varNames)};
                TB.(B{i}) = convert(B{i+1},thisType);
            end
        else
            TB=defaultValues; % Start with default, then replace what is in B.
            defaultB = setdiff(varNames,B(1:2:end))';
            TB = fillTable(varTypes,varNames,TB,B);
        end
    end
end
%% With these tables we can use the builtin functions to create a "designmatrix" for A and B
[~,varLocs] = ismember(TA.Properties.VariableNames,m.VariableNames);

terms = m.Formula.FELinearFormula.Terms(:,varLocs);
aTerms = terms;
dvCoding = lm.dummyVarCoding(m);
[vA,terms,cols2vars,cols2terms,colNames,termNames]  = classreg.regr.modelutils.designmatrix(TA,'Model',aTerms, ...
    'DummyVarCoding', dvCoding, ...
    'CategoricalVars',m.VariableInfo.IsCategorical(varLocs), ...
    'CategoricalLevels',m.VariableInfo.Range(varLocs)); %#ok<ASGLU>
if any(ismissing(vA))
    error('The A condition contains missing values, suggesting you specified levels that do not exist? (Case sensitive categoricals?)');
end
% Remove the terms that depend on the (arbitrary) default valuse
for i=1:numel(defaultA)
    defaultTerms = find(~cellfun(@isempty,regexp(termNames,defaultA{i})));
    if ~isempty(defaultTerms)
        out = ismember(cols2terms,defaultTerms);
        vA(out) = 0;
    end
end

%% Same for B if requested.
if nargin <3 || isempty(B)
    vB = zeros(size(vA));
    vB(1) = 1; % Intercept only
else
    bTerms = terms;
    [vB,~,cols2vars,cols2terms,colNames,termNames]  = classreg.regr.modelutils.designmatrix(TB,'Model',bTerms, ...
        'DummyVarCoding', dvCoding, ...
        'CategoricalVars',m.VariableInfo.IsCategorical(varLocs), ...
        'CategoricalLevels',m.VariableInfo.Range(varLocs)); %#ok<ASGLU>
    if any(ismissing(vB))
        error('The B condition contains missing values, suggesting you specified levels that do not exist? (Case sensitive categoricals?)');
    end
    for i=1:numel(defaultB)
        defaultTerms = find(~cellfun(@isempty,regexp(termNames,defaultB{i})));
        if ~isempty(defaultTerms)
            out = ismember(cols2terms,defaultTerms);
            vB(out) = 0;
        end
    end
end

% The contrast is the difference between the two vectors
% This removes the influence of the default values.
v =vA-vB;
if scale
    v = v ./(0.5*sum(abs(v)));
end
end

function TX= fillTable(varTypes,varNames,TX,X)
for i = 1:numel(varNames)
    ix= find(ismember(X(1:2:end),varNames{i}));
    if numel(ix)==1
        thisVal = convert(X{2*(ix-1)+2},varTypes{i});
        TX.(varNames{i}) = thisVal;
    elseif numel(ix)==0
        %nothing to do (probably a random effect)
    else
        error('?');
    end
end
end

function val = convert(val,type)
if ~isa(val,type)
    if isa(val,'char') && strcmpi(type,'categorical')
        val = categorical({val});
    elseif isa(val,'char') ||isa(val,'string') && strcmpi(type,'cell')
        val = {char(val)};
    else
        val =feval(type,val);
    end
end
end

