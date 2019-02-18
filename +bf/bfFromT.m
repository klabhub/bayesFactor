function bf10 = bfFromT(T,df,n)
% Estimate the Bayes Factor from a T statistic associated with a 
% two independent sample T-test, using the BIC
% approximation of Falukenberry 2018. 
% Note that this is provided for completeness only; the bf.ttest function
% provides a better calculation of the Bayes factor (also based on T, df, and n)
% that does not make the BIC approximation.
%
% INPUT
% F - The T-value of the test
% df1 - Treatment degrees of freedom
% n - Sample size.
% 
% OUTPUT
% bf10 - Bayesfactor for the H1 (effect of treatment) over H0 (no effect)


bf01 = sqrt(n*(1+T.^2/df).^(-n));
bf10 = 1./bf01;