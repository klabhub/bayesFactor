function [T,X] = allConditions(m)
% Return a table in which each row is a specific condition in the lienar model
% Grouping variables are ignored -i.e. fixed effects only.
% INPUT
%  m - The linear (mixed) model
% OUTPUT
% T = table with columns representing the variables in the model
%     and each row correspondong to a specific condition (=combination of
%     factors)
% X = The same conditions, now coded as vectors.
%
% BK - Mar -2021


allVariables = m.Variables(:,m.VariableInfo.InModel);
fixedEffectVariables = allVariables(:,m.Formula.FELinearFormula.PredictorNames);
T= unique(fixedEffectVariables,'rows');

if nargout >1
    X = lm.contrast(m,T);
end    

end