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
if nargin <3
    showEffects = false;
    if nargin<2 || isempty(factors)
        factors =m.anova.Term(2:end);
    end
end

if ischar(factors)
    if strcmpi(factors,'*')
        factors =m.anova.Term; % All including the intercept
    else
        factors = {factors};
    end
end

fprintf('%s\n',m.Formula.char);
for f=1:numel(factors)
    factor = factors{f};
    stay = strcmpi(m.anova.Term,factor);
    fmt = '\t %s: F(%d,%d)=%3.2f,p=%3.3g';
    vars = {factor,m.anova.DF1(stay,1),m.anova.DF2(stay,1),m.anova.FStat(stay,1),m.anova.pValue(stay,1)};
    if m.anova.pValue(stay,1) <0.05
        style = 2; % Error output stream ; red
    else
        style =1; % stdout; 
    end
    if showEffects
        fmt = cat(2, fmt ,' (%3.3f )');
        fe = m.fixedEffects;
        vars = cat(2,vars,{fe(stay)});
    end
    fprintf(style,[fmt '\n'],vars{:});    
end
