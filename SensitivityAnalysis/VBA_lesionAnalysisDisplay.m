function VBA_lesionAnalysisDisplay(results,nodes)

nx = numel(results.lesion);
nObs = numel(results.normal.kernel) ;
nu = numel(results.normal.kernel(1).landmarks);


if nargin > 1

cpt=1;
for iObs = 1:nObs
	for iInput = 1:nu
        
        subplot(nu,nObs,cpt);
        
        base = results.normal.kernel(iObs).landmarks(iInput).aMax;     
        for iNode = 1:nx
         effect(iNode) = results.lesion(iNode).kernel(iObs).landmarks(iInput).aMax;
        end
        beta  = 1*(effect-base)/(base);
        
        theta = zeros(1,results.out.dim.n_theta);
        connect  = grapher_connectivityPattern(results.out,theta);
        
        grapher_staticDcmDisplay(nodes,beta,connect,[]);        
        
        cpt = cpt+1;
    end  
end
else
    
cpt = 1;
for iRegion = 1:nx
    for iObs = 1:nObs
        subplot(nx,nObs,cpt);
        hold on
        effects = [[results.normal.kernel(iObs).landmarks.aMax]' [results.lesion(iRegion).kernel(iObs).landmarks.aMax]'];
        effects_sd = sqrt([[results.normal.kernel(iObs).sigma.landmarks]' [results.lesion(iRegion).kernel(iObs).sigma.landmarks]']);
        h=bar(effects);
        errorbar(get_xBars(h),effects(:),effects_sd(:),'LineStyle','none','Color','k');
        set(h(1),'FaceColor',[32 56 84]/255);
        set(h(2),'FaceColor',[74 126 187]/255);
        for i=1:nu, labels{i} = sprintf('u_%d',i); end
        set(gca,'Xtick',1:nu,'XTickLabel',labels);
        xlabel('kernel')
        ylabel('max amplitude')
        title(sprintf('lesion on node %d',iRegion));
        cpt=cpt+1;
%         ylim([-2.5 2.5])
    end
end
end
end

function x=get_xBars(h)
    f=findall(h,'Type','Patch');
    xall = get(f,'XData');
    xall = horzcat(xall{:});
    xall = xall(:);
    xall=reshape(xall,4,numel(xall)/4);
    x = mean(xall);
    
end