function testSuite %% bfTest

y = [0.7 0.71 0.68 0.57 0.6 0.68 0.57 0.67 0.62 0.75
    0.53 0.52 0.82 0.78 0.6 0.67 0.76 0.65 0.49 0.63
    0.87 0.72 0.85 0.82 0.81 0.88 0.84 0.94 0.94 0.68
    0.97 0.9 0.67 0.66 0.6 0.61 1.07 0.89 0.71 0.6]';
[ nrSubjects,nrConditions]= size(y);
data = table(y(:),[true(2*nrSubjects,1);false(2*nrSubjects,1)],repmat([true(nrSubjects,1);false(nrSubjects,1)],[2 1]),repmat((1:nrSubjects)',[nrConditions 1]),'VariableNames',{'rt','ori','freq','subject'});

rouderFigure5(data);

end



function rouderFigure5(data)
% Runs the simulations of Rouder et al. 2009 in Figure 5, to illustrate the
% differences between fixed and random effects.

% Both ori and freq fixed
bfFullFixed= bf.lme(data,'rt',{'ori','freq'},{},'interactions','all','nDimsForMC',3);
bfBothFixed= bf.lme(data,'rt',{'ori','freq'},{},'interactions','none','nDimsForMC',3);
bfOriFixed= bf.lme(data,'rt',{'ori'},{},'interactions',{'ori:freq'},'nDimsForMC',3);
bfFreqFixed= bf.lme(data,'rt',{'freq'},{},'interactions',{'ori:freq'},'nDimsForMC',3);

bf10  = nan(3,4); % 3 factors (ori,freq,int) and 4 kinds of effects (fixed, mixed, mixed, random)
bf10(1,1) = bfFullFixed/bfFreqFixed; % Main or int effect of ori
bf10(2,1) = bfFullFixed/bfOriFixed; % Main or int effect of freq
bf10(3,1) = bfFullFixed/bfBothFixed; % Interaction


% Ori fixed, freq random
bfFullMixed= bf.lme(data,'rt',{'ori'},{'freq'},'interactions','all','nDimsForMC',3);
bfBothMixed= bf.lme(data,'rt',{'ori'},{'freq'},'interactions','none','nDimsForMC',3);
bfOriMixed= bf.lme(data,'rt',{'ori'},{},'interactions',{'ori:freq'},'nDimsForMC',3);
bfFreqMixed= bf.lme(data,'rt',{},{'freq'},'interactions',{'ori:freq'},'nDimsForMC',3);

bf10(1,2) = bfFullMixed/bfFreqMixed; 
bf10(2,2) = bfFullMixed/bfOriMixed; 
bf10(3,2) = bfFullMixed/bfBothMixed;


% Ori random , freq fixed
bfFullMixed= bf.lme(data,'rt',{'freq'},{'ori'},'interactions','all','nDimsForMC',3);
bfBothMixed= bf.lme(data,'rt',{'freq'},{'ori'},'interactions','none','nDimsForMC',3);
bfOriMixed= bf.lme(data,'rt',{},{'ori'},'interactions',{'ori:freq'},'nDimsForMC',3);
bfFreqMixed= bf.lme(data,'rt',{'freq'},{},'interactions',{'ori:freq'},'nDimsForMC',3);

bf10(1,3) = bfFullMixed/bfFreqMixed; 
bf10(2,3) = bfFullMixed/bfOriMixed; 
bf10(3,3) = bfFullMixed/bfBothMixed; 

% Both random
bfFullRandom = bf.lme(data,'rt',{},{'ori','freq'},'interactions','all','nDimsForMC',3);
bfBothRandom= bf.lme(data,'rt',{},{'ori','freq'},'interactions','none','nDimsForMC',3);
bfOriRandom = bf.lme(data,'rt',{},{'ori'},'interactions',{'ori:freq'},'nDimsForMC',3);
bfFreqRandom= bf.lme(data,'rt',{},{'freq'},'interactions',{'ori:freq'},'nDimsForMC',3);

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
function rouderFigure4(data)
% Runs the simulations of Rouder et al. 2009 in Figure 4.
% This basically shows the ability to extract Main and Interaction effects
% in a 2-way ANOVA.

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
nrSets = 1000;
tic
X= classreg.regr.modelutils.designmatrix(data,'intercept',false,'responsevar','rt','DummyVarCoding','effects','PredictorVars',{'ori','freq'},'model','interactions');
parfor j=1:nrEffects
    rt = X *effects(j,:)';
    for i=1:nrSets
        tmp  =data;
        tmp.rt =rt + randn([size(X,1) 1]);
        bfFull = bf.lme(tmp,'rt',{'freq','ori'},'interactions','all','sharedPriors','within','nDimsForMC',3);
        bfMain = bf.lme(tmp,'rt',{'freq','ori'},'interactions','none','sharedPriors','within','nDimsForMC',3);
        bfOri(i,j)  = bf.lme(tmp,'rt',{'ori'});
        bfFreq(i,j)  = bf.lme(tmp,'rt',{'freq'});
        bfInteraction(i,j) = bfFull/bfMain;
    end
end

%%
figure(1);clf
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
