function VBA_sensitivityAnalysisDisplay(results,nodes)


nResps = size(results.timeseriesBeta,2);
n_u = results.out.dim.u;
idxTheta = results.theta.idx;

cpt = 1;
figure('Position',[50,50, 350*n_u, 350*nResps])

for iObs = 1:nResps
    for iu = 1:n_u  
        
        betas = reshape(results.kernelLandmarksBeta(iObs,iu,:),1,1+numel(idxTheta)) ;
        betasSe = reshape(results.kernelLandmarksBetaSe(iObs,iu,:),1,1+numel(idxTheta)) ;
        
        f=subplot(nResps,n_u,cpt);
        
        if nargin == 1 
        %% simple bar graph
%         bar(betas,'FaceColor',[53 91 135.5]/255); 
        bar(betas(2:end),'FaceColor',[53 91 135.5]/255); 
        hold on;
        errorbar(betas(2:end),betasSe(2:end),'LineStyle','none','Color','k');

        ylabel(sprintf('response %d',iObs));
        xlabel('connection')
%         set(gca,'XTickLabel',{'0',results.theta.lbl{:}});
        set(gca,'XTickLabel',{results.theta.lbl{:}});
        title(sprintf('kernel u_%d',iu));
        
        else
        %% DCM display
        theta = zeros(1,results.out.dim.n_theta);
        theta(idxTheta) = betas(2:end)/betas(1);
        connect  = grapher_connectivityPattern(results.out,theta);
        connect = connect/max(abs(connect(:)));
        
        grapher_staticDcmDisplay(nodes,[],connect,f)
%         colormap(flipud(cbrewer('div','RdYlGn',100)));


        end

        cpt = cpt +1;
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
