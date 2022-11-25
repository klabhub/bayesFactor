function [v,TA,TB] = contrast(m,A,B,defineDifference,scale)
% Specify two conditions A and B using cell arrays of parm/value pairs to
% retrieve the contrast weights vector A-B.
%
% INPUT
% m - The linear model
% A - Condition A - A cell array of variable/value pairs.
% B - Condition B - B cell array of variable/value pairs.Can be empty, in
%       which case the function returns the vector defining A relative to
%       the intercept only model
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
assert(islogical(scale)&& islogical(defineDifference),"scale and defineDifference parameters should be logical values");

import lm.*
if isa(A,'double') && isa(B,'double')
    % Both lready specfied as numeric contrasts
    v= A-B;
    TA= table;
    TB =table;
    return;
elseif isa(A,'table')
    % Use the conditions specified in the table.
    TA =A;    
else
    %% Create tables that define the conditions specified by A
    TA = fillTable(m,A);    
end

% Same for condition B
if nargin >2 && ~isempty(B) 
    if isa(B,'table')
        TB = B;        
    else
        if defineDifference
            TB = TA;  % Start with TA, replace what is in B
            varTypes  = m.VariableInfo.Class(m.VariableInfo.InModel);
            varNames = m.VariableNames(m.VariableInfo.InModel);  
            defaultB = setdiff(varNames,cat(2,A(1:2:end),B(1:2:end)))';
            for i=1:2:numel(B)
                thisType = varTypes{strcmpi(B{i},varNames)};
                TB.(B{i}) = convert(B{i+1},thisType);
            end
        else
             TB = fillTable(m,B);                 
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
end

% The contrast is the difference between the two vectors
% This removes the influence of the default values.
v =vA-vB;
if scale
    v = v ./(0.5*sum(abs(v)));
end
end

function TX= fillTable(m,propValSpecs)
    % Provide a linear model  (m) and a cell array specifying
    % property/value pairs that define a condition (i.e., one side of the
    % contrast)
    varTypes  = m.VariableInfo.Class(m.VariableInfo.InModel);
    varNames = m.VariableNames(m.VariableInfo.InModel);
    varRange =m.VariableInfo.Range(m.VariableInfo.InModel);
    nrVars= numel(varTypes);
    TX = table('Size',[1 nrVars],'VariableType',varTypes,'VariableNames',varNames);
    for i=1:nrVars
        [tf,ix] =ismember(varNames{i},propValSpecs(1:2:end));
        if tf
            value = propValSpecs{2*ix};
        else
            switch (varTypes{i})
                case {'categorical','string'}
                    value = varRange{i}(1); % First in range is the default
                case 'double'
                    value = 0; % A continuous variable that was not specified in the contrast (not in X) - set to zero 
                otherwise
                    error('Non categorical variable (%s)... not sure what to do here...',varNames{i});
            end
        end
        TX.(varNames{i}) =value;
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

