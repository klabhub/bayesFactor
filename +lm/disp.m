function disp(m,factors,showEffects)
% Convenience disp function to show Anova results of a Linear model in 
% standard notation for easy copy and paste.
%
% INPUT
%  m =  a linear model
% factors = a cell array of factors whose stats are to be shown. Defaults
% to all but the Intercept.
% showEffects = Show effect sizes too.[false]
% OUTPUT
%
% BK - Feb 2020
if nargin<2 || isempty(factors)
        factors =m.anova.Term(2:end);
end
if nargin <3
    showEffects = false;
end

if ischar(factors)
    if strcmpi(factors,'*')
        factors =m.anova.Term; % All including the intercept
    else
        factors = {factors};
    end
end

fprintf('%s\n',m.Formula.char);
[partialEta,partialEtaLB,partialEtaUB] = lm.partialEtaSquared('lm',m);
eta = ['partial ' char(hex2dec('03B7')) char(178)];
precision = '%3.2g';

for f=1:numel(factors)
    factor = factors{f};
    stay = strcmpi(m.anova.Term,factor);
    
    fmt = ['\t %s: F(%d,%d)=' precision ',p=' precision ' ' eta ' = ' precision ' CI: [' precision ',' precision ']'];    
    vars = {factor,m.anova.DF1(stay,1),m.anova.DF2(stay,1),m.anova.FStat(stay,1),m.anova.pValue(stay,1),partialEta(stay),partialEtaLB(stay),partialEtaUB(stay)};
    if m.anova.pValue(stay,1) <0.05
        style = 2; % Error output stream ; red
    else
        style =1; % stdout; 
    end
    if showEffects
        fmt = cat(2, fmt ,[' (Effect: ' precision ' CI [' precision ' ' precision '] )']);
        fe = m.fixedEffects;
        vars = cat(2,vars,{fe(stay),m.Coefficients.Lower(stay),m.Coefficients.Upper(stay)});
    end
    fprintf(style,[fmt '\n'],vars{:});    
end
