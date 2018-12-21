function rouderFigures(figs,nrSets)
% This function is a test of the BayesFactor package implementation.
% It runs two analyses that are reported in the Rouder et al. 2012 paper.
% and produces the figures from that paper.
% Soruce :
%    Rouder, J.N., Morey, R.D., Speckman, P.L., and Province, J.M. (2012).
%       Default Bayes factors for ANOVA designs. J. Math. Psychol. 56, 356–374.
%
%  INPUT
% figs = Which figure to generate. Defaults  to both [4 5].
% nrSets = how many bootstrap sets for figure 4. [10] . Set to 1000 for a
% robust rosult.
%
% BK - 2018

if nargin<2
    nrSets =10; % WQuick and dirty.
    if nargin <1
        figs = [4 5];
    end
end

% Quick simulation illustrating the issue with having too many random
% effects and interactions.
if any(figs==5)
    rouderFigure5;
end
% This will take a while as it bootsraps the BF for different effect sizes.
% set nrSets =10 for a quick and dirty check, 1000 to really run what
% Rouder et al did.
if any(figs==4)
    rouderFigure4(nrSets);
end

if any(figs==2)
    rouderFigure2;
end
end


function rouderFigure5
% Runs the simulations of Rouder et al. 2009 in Figure 5, to illustrate the
% differences between fixed and random effects.
load rouder2012Data
bf= bayesFactor;
% Both ori and freq fixed
bfFullFixed= bf.linearMixedModel(data,'rt~ori*freq');
bfBothFixed= bf.linearMixedModel(data,'rt~ori+freq');
bfOriFixed= bf.linearMixedModel(data,'rt~ori +ori:freq');
bfFreqFixed= bf.linearMixedModel(data,'rt~freq+ ori:freq');

bf10  = nan(3,4); % 3 factors (ori,freq,int) and 4 kinds of effects (fixed, mixed, mixed, random)
bf10(1,1) = bfFullFixed/bfFreqFixed; % Main or int effect of ori
bf10(2,1) = bfFullFixed/bfOriFixed; % Main or int effect of freq
bf10(3,1) = bfFullFixed/bfBothFixed; % Interaction


% Ori fixed, freq random
bfFullMixed= bf.linearMixedModel(data,'rt~ori*freq','treatAsRandom',{'freq'});
bfBothMixed= bf.linearMixedModel(data,'rt~ori+freq','treatAsRandom',{'freq'});
bfOriMixed= bf.linearMixedModel(data,'rt~ori +ori:freq','treatAsRandom',{'freq'});
bfFreqMixed= bf.linearMixedModel(data,'rt~freq+ ori:freq','treatAsRandom',{'freq'});

bf10(1,2) = bfFullMixed/bfFreqMixed;
bf10(2,2) = bfFullMixed/bfOriMixed;
bf10(3,2) = bfFullMixed/bfBothMixed;


% Ori random , freq fixed
bfFullMixed= bf.linearMixedModel(data,'rt~ori*freq','treatAsRandom',{'ori'});
bfBothMixed= bf.linearMixedModel(data,'rt~ori+freq','treatAsRandom',{'ori'});
bfOriMixed= bf.linearMixedModel(data,'rt~ori +ori:freq','treatAsRandom',{'ori'});
bfFreqMixed= bf.linearMixedModel(data,'rt~freq+ ori:freq','treatAsRandom',{'ori'});

bf10(1,3) = bfFullMixed/bfFreqMixed;
bf10(2,3) = bfFullMixed/bfOriMixed;
bf10(3,3) = bfFullMixed/bfBothMixed;

% Both random
bfFullRandom =bf.linearMixedModel(data,'rt~ori*freq','treatAsRandom',{'freq','ori'});
bfBothRandom= bf.linearMixedModel(data,'rt~ori+freq','treatAsRandom',{'freq','ori'});
bfOriRandom= bf.linearMixedModel(data,'rt~ori+ori:freq','treatAsRandom',{'freq','ori'});
bfFreqRandom= bf.linearMixedModel(data,'rt~freq+ori:freq','treatAsRandom',{'freq','ori'});

bf10(1,4) = bfFullRandom/bfFreqRandom;
bf10(2,4) = bfFullRandom/bfOriRandom;
bf10(3,4) = bfFullRandom/bfBothRandom;


figure(5);
clf;

h = bar(1:3,bf10);
set(gca,'XTick',1:3,'XTIckLabel',{'Orientation','Frequency','Interaction'},'YScale','Log','YTick',[0.1 1 10 100 ])
set(gca,'Ylim',[0.1 200])
for i=1:4
    h(i).FaceColor=  (i-1)*0.3333*ones(1,3);
    h(i).BaseValue = 1;
end
legend(h,{'Both Fixed','Orientation Fixed','Frequency Fixed','Both Random'},'Location','NorthEast')
ylabel 'Bayes Factor vs. Null Model'

end
%% Rouder figure 4
function rouderFigure4(nrSets)
% Runs the simulations of Rouder et al. 2009 in Figure 4.
% This basically shows the ability to extract Main and Interaction effects
% in a 2-way ANOVA.
load rouder2012Data
bf=bayesFactor;
effects = [0   0   0
    0.2 0   0
    0.5 0   0
    1   0   0
    0.2 0.4 0
    0.5 0.4 0
    1   0.4 0
    0.2 0.2   0
    0.5 0.5   0
    1   1     0
    0.4 0.4 0.2
    0.4 0.4 0.5];

nrEffects = size(effects,1);
% Create a design matrix from the data, only to simulate fake rt's with
% different effects
X= classreg.regr.modelutils.designmatrix(data,'intercept',false,'responsevar','rt','DummyVarCoding','effects','PredictorVars',{'ori','freq'},'model','interactions');
for j=1:nrEffects
    rt = X *effects(j,:)';
    for i=1:nrSets
        tmp  =data;
        tmp.rt =rt + randn([size(X,1) 1]);
        bfFull = bf.linearMixedModel(tmp,'rt~ori*freq');
        bfMain = bf.linearMixedModel(tmp,'rt~ori+freq');
        bfOri(i,j)  = bf.linearMixedModel(tmp,'rt~ori');
        bfFreq(i,j)  = bf.linearMixedModel(tmp,'rt~freq');
        bfInteraction(i,j) = bfFull/bfMain;
    end
end

%%
figure(4);clf
effectIx = {1:4,5:7,8:10,11:12};
for i=1:4
    thisPlotEffects  = effectIx{i};
    subplot(1,4,i)
    if i==4
        x = effects(thisPlotEffects,3);
    else
        x = effects(thisPlotEffects,1);
    end
    [meOri]= median(bfOri(:,thisPlotEffects));
    hOri=  plot(x,meOri,'o-','MarkerSize',10,'Color','k','MarkerFaceColor','k');
    hold on
    meFreq = median(bfFreq(:,thisPlotEffects));
    hFreq=  plot(x,meFreq,'o-','MarkerSize',10,'Color','k','MarkerFaceColor',[0.7 0.5 0]);
    meInt = median(bfInteraction(:,thisPlotEffects));
    
    hInt=  plot(x,meInt,'o-','MarkerSize',10,'Color','k','MarkerFaceColor','w');
    if i==1
        legend([hOri, hFreq,hInt],{'Orientation','Frequency','Interaction'},'Location','NorthWest');
    end
    set(gca,'YScale','Log','YLim',[0.1 1e6],'YTick',[0.1 1 10 100 1000 1e6])
    ylabel 'Median Bayes Factor for Effect'
end
end

function rouderFigure2
% Shows the relationship between T test results, the Bayes Factor and
% sample size.
T = 1:0.05:10;
N = round(logspace(log10(5),log10(5000),100));
criticalT = nan(numel(N),3);
cntr =0;
for n=N
    cntr= cntr+1;
    for t= T
        % We call the bayesFactor.ttest function directly with
        % the results of a t-test (e.g. as returned by ttest.m)
        stats.tstat = t;
        stats.df    = n-1;
        stats.N     = n;  
        stats.tail  = 'both';
        stats.p     = NaN; % no tused
        % Call the class member
        bf = bayesFactor.ttest([],'stats',stats);
        % Check whether we;ve reach any of the thresholds.
        % If so, store.
        if bf>3 && isnan(criticalT(cntr,1))
            criticalT(cntr,1) = t;
        end
         if bf>10 && isnan(criticalT(cntr,2))
            criticalT(cntr,2) = t;
         end
         if bf>30 && isnan(criticalT(cntr,3))
            criticalT(cntr,3) = t;
            break; % Goto next n
         end        
    end
end

%%
figure(2);
plot(N,criticalT)
xlabel 'Sample Size'
ylabel 'Critical t-value'
legend('BF=3','BF=10','BF=30')
set(gca,'XScale','Log','YLim',[1 6],'YTick',2:6,'XTick',[5 20 50 200 1000 5000],'XLim',[4 6000])
end
