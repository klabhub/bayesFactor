function [p,stat,df,delta,c,debugStr] = posthoc(m,A,B,predictedDelta,tail)
% Perform a posthoc comparison of condition A and B in a Linear Mixed Model, 
% using coefTest, but using cell arrays of parm/value pairs to specify conditions A and B
% instead of the contrast
%
%
% INPUT
% lm                = The linear model
% A                 =  Cell array specifying condition A
% B                 = Cell array specifying condition B
% predictedDelta    = Predicted difference between A and B. Defaults to 0.
% tail              = Specify the tail of the distribution. 'left','right' (use one sided T-tests)
%                       or 'both' (use F-test).
%
% OUTPUT
% p                 = The p-value associated with the test. (see coefTest)
% stat              =  F for two-tailed tests,  T for one sided
% df                = Error degrees of freedom (see coefTest)
% delta             = The difference.
% contrast          = The contrast used for this test.
% debugStr          = Contrast shown together with coefficient names to
%                       help understand why the contrast is the way it is...
% EXAMPLE
% Fit a model with two fixed-effect predictors and a random
%             effect. Test for the significance of the Cylinders term. The
%             p-value is the same as shown in the anova table.
% load carsmall
% T = table(MPG,Weight,Model_Year,Cylinders);
% T.Cylinders = nominal(T.Cylinders);
% glme = fitglme(T,'MPG ~ Weight + Cylinders + (1|Model_Year)','Distribution','Normal','DummyVarCoding','Effects')
% H0 8 and 6 cylinders are the same.
% p = lm.posthoc(glme,{'Cylinders',8},{'Cylinders',6})
%
% EXAMPLE
% Use cell arrays to select specific cells in a multifactorial design. E.g. in a design
% with stim and valid factors, each with multiple levels
% {'stim','VIS','valid',1} represents the cell in which the stim factor has level VIS and the
% valid factor has level 1.
%
%
% BK -  Jan 2021
% Mar 2021- rewrote to use lm.contrast

if nargin <5
    tail = 'both';
end
if nargin <4
    predictedDelta = 0;
end

c  = lm.contrast(m,A,B); % the linear contrast
delta = c*m.fixedEffects;
debugStr = strcat(m.CoefficientNames', ' : ' , cellstr(num2str(c')));

predictedDelta = predictedDelta(:).*ones(size(delta));
if ~all(size(delta)==size(predictedDelta))
    error('Predicted delta [nrRows nrCols] has to match the expected delta');
end

switch upper(tail)
    case 'BOTH'
        % Two-sided tests - use the built-in F test.
        [p,stat,~,df] = coefTest(m,c,predictedDelta);
    case {'LEFT','RIGHT'}
        cov  = c*m.CoefficientCovariance*c';
        stat = (delta-predictedDelta)/sqrt(cov);
        df  =  m.DFE;
        if strcmpi(tail,'LEFT')
            p = tcdf(stat,m.DFE);
        else
            p = 1-tcdf(stat,m.DFE);
        end
    otherwise
        error('%s is not a valid tail specification',tail);
end


end

