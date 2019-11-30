%%
% Simulate a t statistic that passes 0.05 two-tailed for N=20;
N =20;
p = 0.01;
stats.tstat = tinv(1-p/2,N-1);
stats.df = N-1;
stats.p = p;
stats.tail  ='both';
stats.N = N;
scale = 0.1:0.1:5;
nrScales= numel(scale);
bfT  = nan(1,nrScales);
for i=1:nrScales
    bfT(i) = bf.ttest('stats',stats,'scale', scale(i));
end
    

figure(1) 
R=2;C=2;
subplot(R,C,1);
x = 0:0.1:10;
y = bf.internal.cauchyPdf(x,scale(1));
thisP = tpdf(scale(1),
plot(x,y./sum(y))
xlabel 'Effect Size'
ylabel 'Prior probabiity'
title (['Prior. Scale = ' num2str(scale(1))]);

subplot(R,C,2);
y = bf.internal.cauchyPdf(x,scale(end));
plot(x,y./sum(y))
xlabel 'Effect Size'
ylabel 'Prior probabiity'
title (['Prior. Scale = ' num2str(scale(end))]);
subplot(R,C,[3 4]);
plot(scale,bfT)
hold on
xlabel 'Prior Scale'
ylabel 'Bayes Factor'
title (sprintf('BF for N = %d , T = %2.2f (p= %2.2f)',N,stats.tstat,p));
plot(sqrt(2)/2*ones(1,2),ylim,'r')
plot(sqrt(2)*ones(1,2),ylim,'g')