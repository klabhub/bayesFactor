function plotModMed(results,varargin)
% Plot the results of a moderated mediation analysis in a graph
p =inputParser;
p.addParameter('sigMediatorOnly',false,@islogical);
p.addParameter('effectSize','none',@ischar);
p.parse(varargin{:});


% Adjacency matrix: Treatment Outcome Mediators Moderator
if ischar(results.parms.mediator);mediators={results.parms.mediator};else; mediators = results.parms.mediator;end
nrMediators = numel(mediators);
nodeNames= cat(2,{[results.parms.treatment] [ results.parms.outcome]},mediators);
TREATMENT =1;
OUTCOME =2;
nrModeratorValues= size(results.a,2);
styles = {':','-'};
props.LineWidth =2;
props.ArrowSize= 20;
props.Interpreter ='none';
props.EdgeFontSize =11;
props.NodeFontSize =11;
props.NodeFontWeight ='bold';
   
for mo =1:nrModeratorValues
    subplot(10,nrModeratorValues,mo);
    isSignificant = prod(sign(results.clim.c(1,mo,:)))>0;
    

    thisStyle = styles(isSignificant+1);
    G= digraph(TREATMENT,OUTCOME,results.c(mo),nodeNames(1:2));    
   edgeLabels = {sprintf('c: %3.2f CI [%3.2f %3.2f]',G.Edges.Weight,results.clim.c(1,mo,:))};
   h =  plot(G,'EdgeLabel',edgeLabels,'XData',[0 1],'YData',[0 0],'LineStyle',thisStyle);
   set(h,props);
   title (sprintf('Moderator %s : %+3.2f',results.parms.moderator,results.moderatorValues(mo)));
    axis off
    subplot(10,nrModeratorValues,(nrModeratorValues+mo):nrModeratorValues:10*nrModeratorValues);
    s = nan(1+2*nrMediators,1);
    t = nan(1+2*nrMediators,1);
    
    weight = nan(1+2*nrMediators,1);
    edgeLabels = cell(1+2*nrMediators,1);
    s(1) = TREATMENT;
    weight(1) = results.cPrime(mo);
    t(1) = OUTCOME;
    isSignificant(1) = prod(sign(results.clim.cPrime(1,mo,:)))>0;
    edgeLabels{1}  =sprintf('c'': %3.2f CI [%3.2f %3.2f]',results.cPrime(mo),results.clim.cPrime(1,mo,:));
    
    for i=1:nrMediators
        isSignificant(1+i) = prod(sign(results.clim.a(i,mo,:)))>0;        
        s(1+i) = TREATMENT;
        t(1+i) = 2+i;
        weight(1+i) = results.a(i,mo);  
        edgeLabels{1+i}  =sprintf('a: %3.2f CI [%3.2f %3.2f]', results.a(i,mo),results.clim.a(i,mo,:));
        
    end
    for i=1:nrMediators
        isSignificant(1+nrMediators+i)= prod(sign(results.clim.b(i,mo,:)))>0;        
        s(1+nrMediators+i) = 2+i;
        t(1+nrMediators+i) = OUTCOME;
        weight(1+nrMediators+i) = results.b(i,mo);        
        edgeLabels{1+nrMediators+i}  =sprintf('b: %3.2f CI [%3.2f %3.2f]', results.b(i,mo),results.clim.b(i,mo,:));
    end
    
    G= digraph(s,t,weight,nodeNames); 
    nodeLabels =cell(G.numnodes,1);
    for n=1:G.numnodes
        if n>2
            nodeLabels{n} = sprintf('%s ab: %3.2f CI [%3.2f %3.2f]', nodeNames{n},results.ab(n-2,mo),results.clim.ab(n-2,mo,:));
            switch p.Results.effectSize
                case 'none'
                case 'kappa2'
                    nodeLabels{n} = sprintf('%s k2 %3.2g CI [%3.2g %3.2g]', nodeLabels{n},results.kappa2(n-2,mo),results.clim.kappa2(n-2,mo,:));
            end
        else
            nodeLabels{n} = nodeNames{n};
        end
    end
            
    G.Nodes= addvars(G.Nodes,[true; true; prod(sign(results.clim.ab(:,mo,:)),3)>0],...
                                            [0;0;abs(results.ab(:,mo))],...
                                            nodeLabels, 'NewVariableNames',{'abSignificant','ab','labels'});
    
    if p.Results.sigMediatorOnly                
       nodeId = findnode(G,nodeNames(~G.Nodes.abSignificant));  
       in=[];out=[];
       for n=nodeId'
        in = [in; inedges(G,n)]; %#ok<AGROW>
        out = [out; outedges(G,n)]; %#ok<AGROW>
       end
        isSignificant([in;out])=[];    
        edgeLabels([in ;out]) = [];
       G = rmnode(G,nodeId);                          
    end
    
    thisStyle = styles(isSignificant+1);
    
    X = [0 1 0.5*ones(1,G.numnodes-2)];
    Y = [0 0 (1:G.numnodes-2)];
    [~,ordered] = sort(G.Nodes.ab(3:end),'ascend');
    Y  = Y([1 2 ordered'+2]);
    
   h =  plot(G,'EdgeLabel',edgeLabels,'NodeLabel',G.Nodes.labels,'XData',X,'YData',Y,'LineStyle',thisStyle');
    set(h,props);
    axis off
    ylim([-0.5 G.numnodes-2+0.5])
end

end