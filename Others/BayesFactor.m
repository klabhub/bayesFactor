function varargout = BayesFactor(varargin)
% BAYESFACTOR MATLAB code for BayesFactor.fig
%      BAYESFACTOR, by itself, creates a new BAYESFACTOR or raises the existing
%      singleton*.
%
%      H = BAYESFACTOR returns the handle to a new BAYESFACTOR or the handle to
%      the existing singleton*.
%
%      BAYESFACTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BayesFactor.M with the given input arguments.
%
%      BAYESFACTOR('Property','Value',...) creates a new BAYESFACTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BayesFactor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BayesFactor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BayesFactor

% Last Modified by GUIDE v2.5 04-Jun-2013 14:23:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @BayesFactor_OpeningFcn, ...
    'gui_OutputFcn',  @BayesFactor_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end

% --- Executes just before BayesFactor is made visible.
function BayesFactor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BayesFactor (see VARARGIN)

% Choose default command line output for BayesFactor
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
UpdateData(handles)
UpdateResults(handles)
% UIWAIT makes BayesFactor wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = BayesFactor_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end


% UpdateData is called if the user changes the experimental data. It then
% calls UpdateResults.
function ntrialsY_Callback(hObject, eventdata, handles)
UpdateData(handles)
end
function avmY_Callback(hObject, eventdata, handles)
UpdateData(handles)
end
function ntrialsX_Callback(hObject, eventdata, handles)
UpdateData(handles)
end
function avmX_Callback(hObject, eventdata, handles)
UpdateData(handles)
end


% UpdateResults is called if the user changes the parameters of the prior.
function DeltaPiMin_Callback(hObject, eventdata, handles)
UpdateResults(handles)
end
function DeltaPiMax_Callback(hObject, eventdata, handles)
UpdateResults(handles)
end
function DeltaPiPDF_Callback(hObject, eventdata, handles)
UpdateResults(handles)
end



function UpdateData(handles)
% Called every time the data changes, ie MY, NY, MO, NO.
% Updates the GUI and resets the prior:


% Read off values from the GUI:
[MY,NY,MO,NO,observedD,mu,Dmax,Dmin,Deltamin,Deltamax] = GetData(handles);


% Write these values into the GUI so the user can see how things change:
set(handles.obtained,'string',num2str(observedD, '%7.5g'));
set(handles.PiMean,'string',num2str(mu, '%7.5g'));
set(handles.propYoung,'string',num2str(MY/NY, '%7.5g'));
set(handles.propOld,'string',num2str(MO/NO, '%7.5g'));

% Set default prior to be a half-Gaussian with SD = half Dmax. Again, user
% can override this if they want.
sigma = 0.5*abs(Deltamin);
strlong = sprintf('exp(-0.5*x.^2/%f^2) .* (x>=0)',sigma);
set(handles.DeltaPiPDF,'string',strlong);

% Now see how the results change:
UpdateResults(handles);

end

function [MY,NY,MO,NO,observedD,mu,Dmax,Dmin,Deltamin,Deltamax] = GetData(handles)
% Read off values from the GUI:
% Total number of trials done by all younger participants:
NY=str2num(get(handles.ntrialsX,'string'));
% Total number of SUCCESSFUL trials done by all younger participants:
MY=str2num(get(handles.avmX,'string'));
% Total number of trials done by all older participants:
NO=str2num(get(handles.ntrialsY,'string'));
% Total number of SUCCESSFUL trials done by all older participants:
MO=str2num(get(handles.avmY,'string'));

% Check these values make sense; if not, reset them to something sensible
if length(NY)>1 || isempty(NY) || NY<=0 ||  isnan(NY) || isinf(NY) || ~isreal(NY)
    set(handles.ntrialsX,'string',MY);
end
if MY>NY
    set(handles.avmX,'string',NY);
end
if length(NO)>1 || isempty(NO) || NO<=0 ||  isnan(NO) || isinf(NO) || ~isreal(NO)
    set(handles.ntrialsY,'string',MO);
end
if MO>NO
    set(handles.avmY,'string',NO);
end
if length(MY)>1 || isempty(MY) || MY<0 ||  isnan(MY) || isinf(MY) || ~isreal(MY)
    set(handles.avmX,'string',NY);
end
if length(MO)>1 || isempty(MO) || MO<0 ||  isnan(MO) || isinf(MO) || ~isreal(MO)
    set(handles.avmY,'string',NO);
end
% Read off the values from the GUI again in case we changed them:
NY=str2num(get(handles.ntrialsX,'string'));
MY=str2num(get(handles.avmX,'string'));
NO=str2num(get(handles.ntrialsY,'string'));
MO=str2num(get(handles.avmY,'string'));


% Observed difference in proportion correct, D:
observedD = MY/NY - MO/NO;

% Observed mean proportion correct, averaged over all participants and
% ASSUMED to be equal to the mean success probability averaged across
% participants:
mu = ( MY + MO ) / (NY + NO);

% Total number of successful trials summed over both age-groups:
M = MY+MO;

% Given our assumption that the true mean success probability (averaged across age-groups) is the
% observed mean proportion of successes (averaged across age-groups), the
% data also constrain the possible range of differences, as follows:
Deltamin = max([ -2*mu 2*(mu-1) ]) ;
Deltamax = min([ 2*(1-mu) 2*mu ]) ;
set(handles.DeltaPiMin,'string',num2str(Deltamin));
set(handles.DeltaPiMax,'string',num2str(Deltamax));

% As well as the above limits on the difference in *probabilities*, we can
% also ask about the maximum possible difference in *proportions*, taking
% into account the sample sizes and the total number correct. I call these
% Dmax and Dmin.
% NB Dmax and Dmin were not actually used in the end code, but I record them for interest.
% Maximum possible difference in proportion correct between age-groups (young minus old),
% given the total numbers of trials performed in each age-group, NY and NO, and the total number of correct trials, M:
Dmax = (M/NY).*(M<=NY) + (1-(M-NY)/NO).*(M>NY);
% Minimum possible difference in proportion correct between age-groups (young minus old),
% given the total numbers of trials performed in each age-group, NY and NO, and the total number of correct trials, M:
Dmin = -M/NO.*(M<=NO) + ((M-NO)/NY-1).*(M>NO);
% (NB if sample sizes are equal, NO=NY, Dmin=-Dmax and Dmin=Deltamin, Dmax=Deltamax).
% Worked example to explain what I mean.
% Suppose MY=0, NY=20; MO=3, NO=10. Mean proportion correct is mu = 0.1. So the limits on the probabilities, piY and piO,
% are 0 and 0.2, ie Deltamin=-0.2 and Deltamax=0.2. In this example
% Dmax=0.15 and Dmin=-0.3. So (i) Dmax is lower than Deltamax. That's
% because of the unequal sample sizes. To get the biggest difference in success *proportions* while keeping
% NY,NO and M the same, we would have to say that the old group had 0
% successes out of 10 trials and the young group had 3 successes out of 20
% trials. That would result in an observed difference in *proportions* of
% 3/20-0/10 = 0.15.
% Also, (ii) Dmin is smaller than Deltamin. In this case, Dmin is actually
% the observed data, where the young group has 0 successes out of 20 trials
% and the old group has 3/10. That's an observed difference in proportions
% of 0/20-3/10=-0.3. But if we assume the young group has a true success
% *probability* of 0, then the old group must have a true success
% probability of 0.2 in order for the mean still to be 0.1




end

function UpdateResults(handles)
% This is called EITHER when new experimental data has been entered (a
% change in MY, NY, MO or NO) OR when the prior has been changed (a change
% in the prior function or its range).

% Read off experimental data from the GUI:
[MY,NY,MO,NO,observedD,mu,Dmax,Dmin,Deltamin,Deltamax]  = GetData(handles);

% Read off current settings for prior from the GUI:
strlong = get(handles.DeltaPiPDF,'string');
ProbDelta = @(x) eval(strlong);


% Calculate normalisation factor for Pr{Delta}
normfactor = quadgk(ProbDelta, Deltamin, Deltamax,'abstol',1e-16','reltol',1e-6,'waypoints',0);

% Draw a bar chart of the data
axes(handles.axesdata);
pY = MY/NY;
pO = MO/NO;
bar(1,pY,'facecol',[0.992 0.918 0.796]);
hold on
bar(2,pO,'facecol',[0.757 0.867 0.776]);
ylim([0 1])
ylabel('proportion correct');
set(handles.axesdata,'xtick',[1 2],'xticklabel',{'Y' 'O'});
title('Data')
% MArk on 95% confidence intervals using binomial statistics
[LoY UpY]=BinoConf_Score(MY,NY);
[LoO UpO]=BinoConf_Score(MO,NO);
errorbar([1 2],[pY pO],[pY pO]-[LoY LoO],[UpY UpO]-[pY pO],'k','linest','none')
hold off


% Draw the plot of prior and posterior probabilities
axes(handles.axes1);
DP = [Deltamin : (Deltamax-Deltamin)/1e3 : Deltamax];
% First plot the prior assumed for Delta-pi:
Prior = ProbDelta(DP)/normfactor;
plot(DP,Prior,'r-')
hold on
xlabel('Difference in success probability, \Delta')
ylabel('Probability')
% Now plot the Pr{data|theory} for different values of Delta-pi
% So that it is visible, normalise it to the same peak as the prior
mx = max(Prior);
pdata = GetProbDiffGivenTheory(DP);
plot(DP,pdata/max(pdata)*mx,'b-');
% Mark on limits:
plot(Deltamin*[1 1],ylim,'k:')
plot(Deltamax*[1 1],ylim,'k:')
% Mark on observed D:
plot(observedD*[1 1],ylim,'b');
hold off



set(handles.Lnull,'string','calculating....')
set(handles.Ltheory,'string','calculating....')
set(handles.BayesFactor,'string','calculating....')
drawnow

% Under the null hypothesis, Delta=0:
Likelihoodnull =  GetProbDiffGivenTheory(0)
set(handles.Lnull,'string',num2str(Likelihoodnull))
drawnow

if Likelihoodnull==0
    disp('such a large difference is vanishingly unlikely to be observedD under the null hypothesis')
    Bayesfactor = Inf;
    set(handles.Ltheory,'string','Not calculated')
else
    Likelihoodtheory = quadgk(@integrand, Deltamin, Deltamax,'abstol',1e-16','reltol',1e-6,'waypoints',0)/normfactor
    set(handles.Ltheory,'string',num2str(Likelihoodtheory))
    Bayesfactor = Likelihoodtheory/Likelihoodnull
end
set(handles.BayesFactor,'string',num2str(Bayesfactor))


    function f = integrand(Delta)
        % Integral you need to do to get Likelihoodtheory
        % integral d(Delta) Pr{ data | Delta etc} Pr{Delta}
        f = GetProbDiffGivenTheory(Delta) .* ProbDelta(Delta);
    end

    function totalprob = GetProbDiffGivenTheory(Delta)
        % NB - Delta can be a vector, in which case the returned prob is
        % also a vector
        piY = mu + Delta/2;
        piO = mu - Delta/2;
        % People do a task with sum-nx trials, probability piY of being correct on any given trial.
        % The average proportion correct is sum( mx_i ) / sum(nx_i).
        
        % We have a fixed value of "observedD", given by
        % observedD = MY/NY - MO/NO;
        % The possible values of MY are:
        poss_MY = [0:NY];
        totalprob = zeros(size(Delta));
        for jp=1:length(poss_MY)
            this_MY = poss_MY(jp);
            % if MY were this, what would MO have to be to get the
            % same/similar  observedD difference?
            % Rearranging observedD = MY/NY - MO/NO;
            this_MO = this_MY*NO/NY - observedD*NO ;
            % Strictly if this_MO is not an integer, I guess you can't
            % get the observed difference in sample probabilities with the
            % MY under consideration
            % See if it is an integer:
            if abs(this_MO - round(this_MO)) < 1e-6 ... % count this as an integer -there might be rounding error
                    && this_MO>=0 && this_MO<=NO
                this_MO = round(this_MO);
                % The probability of getting the observed mean score in the X
                % condition is:
                PrY = binopdf( this_MY, NY, piY);
                % and in Y is
                PrO = binopdf( this_MO, NO, piO);
                totalprob = totalprob + PrY.*PrO;
            end
        end
        if any(isnan(totalprob))
            stop
        end
        
    end



end


function [Lo Up]=BinoConf_Score(m,n,varargin)
%[Lo Up]=BinoConf_Score(m,n,alpha)
%m=hits, n=total of trials, alpha=0.05 for CI 95%, alpha = 0.32 for 68% CI corresponding to +/-1SD of normal
%Binomial confidence interval
%Score confidence interval Edwin B. Wilson (1927)
%According to Agresti and Coull (1998) this is the best confidence interval,
%it yields coverage probabilities close to nominal confidence levels, even
%for very small sample sizes.
%
%Agresti, A. and Coull, B. A. (1998). Approximate is better than "exact"
%for interval estimation of binomial proportions", The American Statistician, 52(2), 119-12
%Wilson, E. B. (1927). Probable inference, the law of succession,
%and statistical inference. J. Amer. Statist. Assoc. 22 209–212.
%

%28-07-2008
%Dr. Ignacio Serrano-Pedraza

if nargin==2
    alpha = 0.05;
else
    alpha = varargin{1};
end

p=m./n;
z=norminv(alpha/2);
A=p+((z^2)./(2*n));
B=z*sqrt((p.*(1-p)+z^2./(4*n))./n);
Lo=(A+B)./(1+z^2./n);
Up=(A-B)./(1+z^2./n);
end