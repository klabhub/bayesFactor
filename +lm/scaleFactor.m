function [scale,units] = scaleFactor(m,mode)
% Return the value and units of a scale factor for effects in linear models
% m - The linear model
% mode:
%       INTERCEPT - scale effects to the intercept grand mean.
%      RANDOMSTD - Scale to teh standard deviation of the random effects.
%   RAW - Do not scale (i.e. scale =1, units ='');
% 
%OUTPUT
% scale= The scale factor by which the caller should divide the effects.
% units = The units of the resulting effect. 
%
% Example: see lm.disp
%BK - Sept 2021

if nargin<2
    mode ='RAW';
end


    if contains(mode,'INTERCEPT','IgnoreCase',true)
        % Scale to the fixed effect intercept
        % How big is the effect relative to the grand mean.
        scale  = 0.01*m.Coefficients.Estimate(strcmpi(m.CoefficientNames,'(Intercept)'));
        units = '%%';
    elseif contains(mode,'RANDOMSTD','IgnoreCase',true)
        % Scale to the standard deviation of the random effects
        % "How big is the effect relative to the variation
        % across subjects?'
        scale = 0.01*std(m.randomEffects);
        units = '%%';
    elseif contains(mode,'RAW','IgnoreCase',true)
        % Effect in raw model units'
        scale =1;
        units = '';
    else
        error('Unknown scaling mode %s' ,mode);
    end


end
