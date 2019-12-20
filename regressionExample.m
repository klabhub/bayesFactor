%% Example with continuous covariates (i.e. regression)

% Compare the Liang et al formula (bf.bfFromR2) with the explciit
% integration implemented in the anova function (so that it can handle
% continuous covariates as well as mixtures of continuous and categorical
% covariates).
x=(1:20)';
nrSamples = numel(x);
nrParameters=1;
noiseLevel = 5;   
for i=1:100
    y = 1*x+noiseLevel*randn([nrSamples 1]);
    [b,~,~,~,stats] = regress(y,[ones(size(x)) x]);    
    R2 =stats(1);
    % Use the Liang et al formula directly
    bfFormula(i) = bf.bfFromR2(R2,nrSamples,nrParameters);
    % Create a table to do integration in the anova function
    T = table(y,x);
    bf10(i)=bf.anova(T,'y~x');
end
%%
figure(1);
clf;
subplot(1,2,1);
plot(bf10,bfFormula,'.')
hold on
set(gca,'XScale','Log','YScale','Log')
plot(xlim,xlim,'k')
axis square;
xlabel 'BF Integrated'
ylabel 'BF Liang Formula'
subplot(1,2,2);
hist(log10(bf10)-log10(bfFormula))
xlabel 'Difference (log10)'
ylabel '# Models'