function [bf,nrSubjects] = bfPower
% Power Analysis for Sequential Bayes Factor designs

minimumBf = 6;
maximumBf = 1/6;

%% Example 1
nrRepeats = 10000;
ns =[20 100];
bf = nan(nrRepeats,2,2); %[H0 H1]

nrNs = numel(ns);
sigmaEffectSize = 0;
meanEffects = [0 0.5];  % H0 H1
sigma = 1;
for h= 1:2
    meanEffectSize = meanEffects(h);
    for r=1:nrRepeats
        for n =1:nrNs
            effectSize  = meanEffectSize+ sigmaEffectSize*randn;
            samplesA  =  effectSize + sigma* randn(ns(n),1);
            samplesB  =  sigma*randn(ns(n),1);
            
            bf(r,h,n) =bfFrom2SamplesT([samplesA samplesB]);
            
        end
    end
end
acceptH0= 100*squeeze(nanmean(bf<maximumBf));
acceptH1= 100*squeeze(nanmean(bf>minimumBf));
%%
figure(1);
clf
for n=1:2
    for h= 1:2
        subplot(2,2,4-2*h+n);
        logBfEdge = -2:0.25:8;
        
        [N,X]= hist(log10(bf(:,h,n)),logBfEdge);%,'Normalization','Probability');
        plot(10.^X,N./sum(N))
        set(gca,'XScale','log');%'XTick',log(bfEdge),'XTickLabels',num2str(bfEdge',2));
        hold on
        plot(([maximumBf maximumBf]),ylim,'k:','LineWidth',2)
        plot(([minimumBf minimumBf]),ylim,'k:','LineWidth',2)
        xlabel 'Bayes Factor '
        ylabel 'Density'
        if h==1
            % H0
            str =['False Pos: ' num2str(acceptH1(h,n)) '% True Neg: ' num2str(acceptH0(h,n)) '% Inconclusive: '  num2str(100-acceptH1(h,n)-acceptH0(h,n)) '%'];
        else
            % H1
            str =['False Neg: ' num2str(acceptH0(h,n)) '% True Pos: ' num2str(acceptH1(h,n)) '% Inconclusive: '  num2str(100-acceptH1(h,n)-acceptH0(h,n)) '%'];
        end
        title(char(['n = ' num2str(ns(n)) ', \delta = ' num2str(meanEffects(h)) ' under H' num2str(h-1)],...
            str) );
    end
end

return;

%% Example 2


nrRepeats = 1000;
nMin = 30;
nStep = 1;
nStop = 500;
bf = nan(nrRepeats,2); %[H0 H1]
nrSubjects = nan(nrRepeats,2); % [H0 H1]

meanEffectSize = 0.5;
sigmaEffectSize = 0.0;




%H1
for r=1:nrRepeats
    for n =nMin:nStep:nStop
        effectSize  = meanEffectSize+ sigmaEffectSize*randn;
        samplesA  =  effectSize + randn(n,1);
        samplesB  =  randn(n,1);
        
        thisBf =bfFrom2SamplesT([samplesA samplesB]);
        if thisBf > minimumBF || thisBf < maximumBF
            bf(r,1) = thisBf;
            nrSubjects(r,1) = n;
            break;
        end
    end
end


%H0
for r=1:nrRepeats
    for n =nMin:nStep:nStop
        %effectSize  = sigmaEffectSize*randn;
        
        samplesA  = 0 + randn(n,1);
        samplesB  = 0 + randn(n,1);
        
        thisBf =bfFrom2SamplesT([samplesA samplesB]);
        if thisBf < maximumBF || thisBf >minimumBF
            bf(r,2) = thisBf;
            nrSubjects(r,2) = n;
            break;
        end
        
    end
end

%%
figure(1);
clf
subplot(2,2,1)

histogram(bf(:,1),'binlimits',[0 100],'Normalization','Probability');
title([ num2str(mean(bf(:,1)>100)) '% > 100 removed']);
subplot(2,2,2)
histogram(bf(:,2),'binlimits',[0 10],'Normalization','Probability');

subplot(2,2,3);
histogram(nrSubjects(:,1),'Normalization','Probability')
hold on;
histogram(nrSubjects(:,2),'Normalization','Probability')
title(['#Samples under H1 and H0: Median: [' num2str(nanmedian(nrSubjects)) '] 80%:  [' num2str(prctile(nrSubjects,80)) ']' ]);
positiveRate = 100*nanmean(bf>minimumBF)
negativeRate = 100*nanmean(bf<maximumBF)

end
function bf = bfFrom1SamplesT(samples,cauchyScale)
if nargin<2
    cauchyScale = sqrt(2)/2;
end
% One-sample T
n = size(samples,1);
T = mean(samples,1)*sqrt(n)/std(samples,0,1);
bf = t1smpbf(T,n,cauchyScale);
end



function bf = bfFrom2SamplesT(samples,cauchyScale)
% BF for 2 sample T-test with unequal numbers and variance.
% samples  - [N 2] each column has the samples for each group. Missing
% values are nan.
%
if nargin<2
    cauchyScale = sqrt(2)/2;
end

dims =size(samples,2);
assert(dims==2,'Two columns with samples needed');
var = nanvar(samples);
n = sum(~isnan(samples));
meanDiff  = diff(nanmean(samples));
T = meanDiff./sqrt(sum(var./n));

bf = t2smpbf(T,n(1),n(2),cauchyScale);
end