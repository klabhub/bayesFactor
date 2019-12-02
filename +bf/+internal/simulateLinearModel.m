function out = simulateLinearModel(lm,effectSize,N,subject)
% Simulate data based on a linear model .
% 
% lm = A linearMixedModel that specifies the base model. This function will
% generate simulated data based on the predictions of this model.
% effectSize = 1: use linear model parameters as estimated in the model
%              0: use effect size 0; simulated the null model.
%              other values scale the std dev of the error term. So 2->
%              noise in the simultaed data has half the stdev.
% 
% N = How many observations to generate. Note that N=1 means
% one observation for each of the combinations of all fixed effects.
% subject = The name of the variable in the model that indicates the
% 'subject'. This is used to select a complete subset of data for a single 
% measurement. In other words,for N=1, the output table will have all the
% data collected for 1 subject. So, for a 3x2 design with 10 subjects you'd
% get 6 rows with N=1. Leaving subject empty would generate a 60 row output table. 
%
% OUTPUT
% out {1} =table with simulated data.
% out{2} = formula of the linear model.

% Create a complete data table for a single subject who behaves exactly as
% the model (no noise).
tbl =lm.Variables;
errorStd = std(lm.residuals);
if ~isempty(subject)
    uSubjects = unique(tbl.(subject));
    if iscell(uSubjects)
        % Char specified 
        stay = strcmpi(tbl.(subject),uSubjects{1});
    else
        % Num specified
        stay = tbl.(subject)==uSubjects(1);
    end
end
tbl = tbl(stay,:);

if effectSize==0
    y = zeros(height(tbl),1);
else
    y = predict(lm,tbl);
    errorStd = errorStd/effectSize; % Defining effectsize relative to stdev of the residuals.   
end
tbl.(lm.ResponseName) = y;
tbl = repmat(tbl,[N 1]); % Replicate N times.
% Now add noise as esimated by the model
tbl.(lm.ResponseName) = tbl.(lm.ResponseName) + errorStd*randn([height(tbl) 1]);
out = {tbl,char(lm.Formula)};
end
