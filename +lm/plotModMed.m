function G = plotModMed(results,varargin)
% Plot the results of a moderated mediation analysis in a graph.
%
% INPUT
% results = struct output from lm.modmed()
% sigMediatorOnly - Show only mediators with a significant indirect path
%                       (false).
% effectSize - Show effect size of the indirect path. Can be  'none', or
%               'kappa2'  (see lm.modmed())
%
% OUTPUT
%  graph only
%  Significant paths are shown with solid lines, non-significant paths with
%  dashed lines.
%
%  Each path is shown with the estimate of the weight and its confidence
%  limits as estimated in the results struct.
%
%  Each column corresponds to the results obtained for one value of the
%  moderator (shown at the top).
%
% TODO:
%   Add better visualization for mediated moderation.
%
% BK  - Sept 2021

p =inputParser;
p.addParameter('sigMediatorOnly',false,@(x) (islogical(x) || isnumeric(x) ));
p.addParameter('anyModIsSig',true,@islogical);
p.addParameter('effectSize','none',@ischar);
p.parse(varargin{:});



if ischar(results.parms.mediator);mediators={results.parms.mediator};else; mediators = results.parms.mediator;end
nrMediators = numel(mediators);
% The nodes are ordered as TRREATMENT OUTCOME MEDIATORS, with the latter
% ordered as the names in mediator
nodeNames= cat(2,{[results.parms.treatment] [ results.parms.outcome]},mediators(:)');
TREATMENT =1;
OUTCOME =2;
nrModeratorValues= size(results.a,2);
styles = {':','-'};  % Nonsig , sig
% Styles for graphs
props.LineWidth =2;
props.ArrowSize= 20;
props.Interpreter ='none';
props.EdgeFontSize =11;
props.NodeFontSize =11;
props.NodeFontWeight ='bold';

G = cell(2,nrModeratorValues);
%Loop over moderator values
for mo =1:nrModeratorValues
    % Show the overall effect (independent of mediators) at the top :
    % c-path
    subplot(10,nrModeratorValues,mo); % 1:10 ratio for overall and mediated effects
    if results.parms.bootstrap >1
        isSignificant = prod(sign(results.clim.c(1,mo,:)))>0;
    else
        isSignificant  = false;
    end
    thisStyle = styles(isSignificant+1);
    G{1,mo}= digraph(TREATMENT,OUTCOME,results.c(mo),nodeNames(1:2));
    
    edgeLabels = sprintf('c: %3.2g ',results.c(1,mo));
    if results.parms.bootstrap
        edgeLabels = sprintf('%s CI [%3.2g %3.2g]',edgeLabels,results.clim.c(1,mo,:));
    end
    h =  plot(G{1,mo},'EdgeLabel',{edgeLabels},'XData',[0 1],'YData',[0 0],'LineStyle',thisStyle);
    set(h,props);
    if ~isempty(results.parms.moderator)
        if isnumeric(results.moderatorValues)
            title (sprintf('Moderator %s : %+3.2f',results.parms.moderator,results.moderatorValues(mo)));
        else
            title (sprintf('Moderator %s : %s',results.parms.moderator,char(results.moderatorValues(mo))));
        end
    end
    axis off
    
    % Show all mediated effects in one graph
    subplot(10,nrModeratorValues,(nrModeratorValues+mo):nrModeratorValues:10*nrModeratorValues);
    s = nan(1+2*nrMediators,1);
    t = nan(1+2*nrMediators,1);
    weight = nan(1+2*nrMediators,1);
    edgeLabels = cell(1+2*nrMediators,1);
    %c' path
    s(1) = TREATMENT;
    weight(1) = results.cPrime(mo);
    t(1) = OUTCOME;
    
    
    edgeLabels{1} = sprintf('c'': %3.2g ',results.cPrime(mo));
    if results.parms.bootstrap
        isSignificant(1) = prod(sign(results.clim.cPrime(1,mo,:)))>0;
        edgeLabels{1} = sprintf('%s CI [%3.2g %3.2g]',edgeLabels{1},results.clim.cPrime(1,mo,:));
    end
    % a-paths
    for i=1:nrMediators
        edgeLabels{1+i} = sprintf('a: %3.2g ',results.a(i,mo));
        if results.parms.bootstrap
            isSignificant(1+i) = prod(sign(results.clim.a(i,mo,:)))>0;
            edgeLabels{1+i} = sprintf('%s CI [%3.2g %3.2g]',edgeLabels{1+i},results.clim.a(i,mo,:));
        end
        s(1+i) = TREATMENT;
        t(1+i) = 2+i;
        weight(1+i) = results.a(i,mo);
        
    end
    % b-paths
    for i=1:nrMediators
        edgeLabels{1+nrMediators+i} = sprintf('b: %3.2g ',results.b(i,mo));
        if results.parms.bootstrap
            isSignificant(1+nrMediators+i) = prod(sign(results.clim.b(i,mo,:)))>0;
            edgeLabels{1+nrMediators+i} = sprintf('%s CI [%3.2g %3.2g]',edgeLabels{1+nrMediators+i},results.clim.b(i,mo,:));
        end
        
        s(1+nrMediators+i) = 2+i;
        t(1+nrMediators+i) = OUTCOME;
        weight(1+nrMediators+i) = results.b(i,mo);
    end
    % Construct the directed graph
    G{2,mo}= digraph(s,t,weight,nodeNames);
    nodeLabels =cell(G{2,mo}.numnodes,1);
    % Add labels with ab estimates and (optionally) effect size
    for n=1:G{2,mo}.numnodes
        if n>2
            nodeLabels{n} = sprintf('%s ab: %3.2g ', nodeNames{n},results.ab(n-2,mo));
            if results.parms.bootstrap>0
                nodeLabels{n} = sprintf('%s CI [%3.2g %3.2g]', nodeLabels{n},results.clim.ab(n-2,mo,:));
            end
            switch p.Results.effectSize
                case 'none'
                case 'kappa2'                    
                    nodeLabels{n} = sprintf('%s k2 %3.2g ', nodeLabels{n},results.kappa2(n-2,mo));
                    if results.parms.bootstrap>0
                        nodeLabels{n} = sprintf('%s CI [%3.2g %3.2g]', nodeLabels{n},results.clim.kappa2(n-2,mo,:));
                    end
            end
        else
            nodeLabels{n} = nodeNames{n};
        end
    end
    % Store in Nodes table
    
        
    if results.parms.bootstrap ==0 && p.Results.sigMediatorOnly~=0
        error('No significance estimates. Canno use sigMediatorOnly');
    end
        
        
    if islogical(p.Results.sigMediatorOnly) &&  p.Results.sigMediatorOnly
        %Use clim  as computed
        abSignificantPerMod = prod(sign(results.clim.ab),3)>0;                       
    elseif p.Results.sigMediatorOnly ~=0
        
        lims = 100*[p.Results.sigMediatorOnly/2 1-p.Results.sigMediatorOnly/2];
        lower = find(results.bs.bins<=lims(1),1,'last');
        upper = find(results.bs.bins>=lims(2),1,'first');
        abSignificantPerMod =  prod(sign(results.bs.ab(:,:,[lower upper])),3)>0;        
    else
        abSignificantPerMod = true(nrMediators,nrModeratorValues,1);
    end
        if p.Results.anyModIsSig
        abSignificant = [true; true; any(abSignificantPerMod,2)];
        else
            abSignificant = [true; true; abSignificantPerMod(:,mo,:)];
        end
    
    G{2,mo}.Nodes= addvars(G{2,mo}.Nodes,...
        abSignificant,...
        [0;0;results.ab(:,mo)],...
        nodeLabels, 'NewVariableNames',{'abSignificant','ab','labels'});
    % Remove mediator nodes whose ab path is not signficant
    if p.Results.sigMediatorOnly
        nodeId = findnode(G{2,mo},nodeNames(~G{2,mo}.Nodes.abSignificant));
        in=[];out=[];
        for n=nodeId'
            in = [in; inedges(G{2,mo},n)]; %#ok<AGROW>
            out = [out; outedges(G{2,mo},n)]; %#ok<AGROW>
        end
        isSignificant([in;out])=[];
        edgeLabels([in ;out]) = [];
        G{2,mo} = rmnode(G{2,mo},nodeId);
    end
    
    % Pick style dependent on significance of the path
    thisStyle = styles(isSignificant+1);
    % Postition the mediators in the graph such that the strongest ab path
    % is on top.
    X = [0 1 0.5*ones(1,G{2,mo}.numnodes-2)];
    [~,ordered] = sort(abs(G{2,mo}.Nodes.ab(3:end)),'ascend');
    [~,y] =sort(ordered);
    Y  = [0 0 y'];
    
    h =  plot(G{2,mo},'EdgeLabel',edgeLabels,'NodeLabel',G{2,mo}.Nodes.labels,'XData',X,'YData',Y,'LineStyle',thisStyle');
    set(h,props);
    axis off
    ylim([-0.5 G{2,mo}.numnodes-2+0.5])
end

sgtitle(sprintf('#Obs: %d #Bootstrap:%d alpha:%3.3f DummyCoding:%s MatchMissing:%d nrMediators:%d RandomEffects:%s',results.lm6.NumObservations, results.parms.bootstrap,results.parms.alpha,results.parms.dummyVarCoding,results.parms.matchMissing,numel(results.parms.mediator),results.parms.randomEffects))
end