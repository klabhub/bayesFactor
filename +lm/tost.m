function [p,T,df,delta,ci,str] = tost(m,A,B,bounds,alpha)
% Post-hoc equivalence test for two conditions in a linear model.
% This is based on the two-one sided tests (TOST) approach. 
% See Schuirmann 1987,Lakens 2017 for a tutorial treatment.
%
% In TOST the null hypothesis is that the effect is larger than some upper
% bound or smaller than some lower bound. The alternative hypothesis is an
% effect that falls within the bounds.
%
% The bounds are chosen by the user to specify the smallest effect size of
% interest (SESOI) and by doing so the user essentially states that effects
% smaller than those bounds are not of interest.
%
% If the p-value returned by this function is smaller than some
% pre-specified significance level (usally 0.05) then you reject the null
% hypothesis that the effect is outside the specified bounds and
% concluded that there are no effects of interest (as defined by the SESOI
% bounds).
%
% The TOST procedure tests two null hypotheses using one-sided T-tests:
% H0 lower: A-B < min(bounds)
% - To *reject* this null we test whether A-B-min(bounds) >0 with a right-tailed T-test
% % H0 upper: A-B > max(bounds)
% - To *reject*this we test whether A-B-max(bounds) < 0 with a left tailed T-test
% If we can reject both of these then we conclude that A-B lies between the bounds
%
% INPUT
% m =  linear model
% A  =  Cell array specifying condition A (see lm.posthoc)
% B =  Cell array specifying condition B (see lm.posthoc)
% bounds= Bounds for the smallest effect size of interest (SESOI). For
% superiority/inferiority equivalence tests, one of the bounds can be inf
% or -inf.
% alpha -  Significance level.
%
% OUTPUT
% p,T,df associated with the test with highest p-value
% delta = the actual difference between A and B
% ci =  1-alpha confidence internval
% str = Char showing the results in a paper ready format.
% 
% EXAMPLE
% Fit a model with two fixed-effect predictors and a random
%             effect. Test whether 8 and 6 cylinder cars have essentially
%           equivalent MPG. Where 'essentially equivalent' is defined as 1 MPG.
%
% load carsmall
% T = table(MPG,Weight,Model_Year,Cylinders);
% T.Cylinders = nominal(T.Cylinders);
% glme = fitglme(T,'MPG ~ Weight + Cylinders + (1|Model_Year)','Distribution','Normal','DummyVarCoding','Effects')
% [p,stat,df,delta] = lm.tost(glme,{'Cylinders',8},{'Cylinders',6},[-1 1]);
% p = 0.23 -> Showing that we can ***not*** conclude they are equivalent.
% (i.e. not reject the null hypothesis that the true differences lie outside the
% [-1 1] interval).
%
% BK - Feb 2021


%%
% BK - Feb 2021
if nargin <5
    alpha = 0.05;
end
if numel(bounds) ~=2
    error('TOST requires both an upper and a lower equivalence bound.');
end
% Is A-B surprisingly larger than the lower bound
[pLB,tLB,df,delta,ciLB,strLB] = lm.posthoc(m,A,B,min(bounds),'right',alpha);
% Or surprisingly smaller than the upper bound
[pUB,tUB,~,~,ciUB,strUB] = lm.posthoc(m,A,B,max(bounds),'left',alpha);
% Keep the least significant of the two one-sided tests.
if pUB < pLB
    p = pLB;
    T = tLB;
    ci  =ciLB;
    str = strLB;
else
    p =pUB;
    T =tUB;
    ci = ciUB;
    str= strLB;
end



if false
    %% Test the code with the example in Lakens 2017, page 357
    %  Results are similar (this uses randn so it cannot be a 
    % perfect match). Lakens has t(182) = 2.69, p = .004
    nControl = 95;
    nOrganic = 89;
    control = 5.25 + 0.95*randn(nControl,1);
    organic = 5.22 + 0.83*randn(nOrganic,1);
    T= table([control;organic],[repmat("control",[nControl 1]);repmat("organic",[nOrganic 1])],'variablenames',{'judgment','food'});
    % Fit an LM
    m = fitlme(T,'judgment~food','DummyVarCoding','effects');
    bound = 0.384;
    dfManual = nControl +nOrganic -2;
    pooledSigma = sqrt(((nControl-1)*var(control) + (nOrganic-1)*var(organic))/dfManual);
    tL = (mean(control)-mean(organic) - -bound)/(pooledSigma*sqrt(1/nControl +1/nOrganic));
    tU = (mean(control)-mean(organic) - bound)/(pooledSigma*sqrt(1/nControl +1/nOrganic));
    pL = 1-tcdf(tL,dfManual);
    pU = tcdf(tU,dfManual);
    
    if pL<pU
        pManual = pU;
        statManual = tU;
    else
        pManual = pL;
        statManual = tL;
    end
    % Now compare the manual results with the klm code.
    [pKlm,statKlm,dfKlm,delta,ci,str] = lm.tost(m,{'food','organic'},{'food','control'},bound*[-1 1]);
    
    fprintf('\n Manual : t(%d)= %3.3f, p= %3.3f,  \n LM Package: %s\n',dfManual,statManual,pManual,str)
end
