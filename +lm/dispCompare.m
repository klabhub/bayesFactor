function dispCompare(c)
% Pretty print the results of the output of the compare function for a
% generalized linear mixed model.
% The line will include the change in AIC in the second compared to the
% first model (hence negative numbers means that the bigger model is a
% better one).
% The stats are the likelihood ratio (Chi-squared distribution).
% 
% The better model will be highlighted with '*'.
% BK - Nov 2021
    % c is the outcome of compare(lm1,lm2). Show the result
if isnan(c.pValue(2))
    % No loglikelihood p-value -  probably not nested. So use AIC only.
    if diff(c.AIC)<0
        % Alternnative model (#2) is best
        first='';second='*';
    else        
        first='*';second='';
    end
    fprintf('%s%s%s vs %s%s%s dAIC %d \n',first,c.Model(1),first,second,c.Model(2),second,round(diff(c.AIC)));
    return;
elseif c.pValue(2)<0.05
    % Alternnative model (#2) is best
    first='';second='*';
else
    second='';first='*';
end
fprintf('%s%s%s vs %s%s%s dAIC %d (X2 (%d) = %3.3f, p= %3.3g)\n',first,c.Model(1),first,second,c.Model(2),second,round(diff(c.AIC)),c.deltaDF(2),c.LRStat(2),c.pValue(2));
end
