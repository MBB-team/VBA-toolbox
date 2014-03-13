function VBA_sensitivityAnalysisDisplay(results,nodes)


nResps = size(results.kernel,2);
n_u = results.out.dim.u;
idxTheta = results.theta.idx;

cpt = 1;
figure('Position',[50,50, 350*n_u, 350*nResps])

for iObs = 1:nResps
    for iu = 1:n_u  
        
        betas = reshape(results.kernelLandmarksBeta(iObs,iu,:),1,1+numel(idxTheta)) ;
%         betasSe = reshape(results.kernelLandmarksBetaSe(iObs,iu,:),1,1+numel(idxTheta)) ;
        
        f=subplot(nResps,n_u,cpt);
        
        if nargin == 1 
        %% simple bar graph
%         bar(betas,'FaceColor',[53 91 135.5]/255); 
        bar(betas(2:end),'FaceColor',[53 91 135.5]/255); 
        hold on;
%         errorbar(betas(2:end),betasSe(2:end),'LineStyle','none','Color','k');

        ylabel(sprintf('response %d',iObs));
        xlabel('connection')
%         set(gca,'XTickLabel',{'0',results.theta.lbl{:}});
        set(gca,'XTickLabel',{results.theta.lbl{:}});
        title(sprintf('kernel u_%d',iu));
        
        else
        %% DCM display
        theta = zeros(1,results.out.dim.n_theta);
        theta(idxTheta) = betas(2:end)/sign(betas(1));
%         theta(idxTheta) =  theta(idxTheta)./sqrt(betasSe(2:end));%/abs(betas(1));
        connect  = grapher_connectivityPattern(results.out,(theta));
        connect = connect/max(abs(connect(:)));
        grapher_staticDcmDisplay(nodes,[],connect,f)
                set(f,'CLim',1*[-1 1])

%         set(f,'Position',[(iu-1)/n_u (nResps-iObs)/nResps .9/n_u .9/nResps])
%         set(f,'visible','off');
%         colormap(flipud(cbrewer('div','RdYlGn',100)));


        end

        cpt = cpt +1;
    end
end


  
% figure('Color','w'); 
% nTheta = numel(results.theta.idx);
% cpt=1;
% for iInput=1:n_u
%     for iTheta = 1:nTheta
%         subplot(n_u,nTheta,cpt);
%         cpt = cpt +1;
%         
%         colors=get_colors(results.thetaArray(:,iTheta));
%         
%         for i=1:min(300,numel(results.kernel))
%     
%             ts = results.kernel(i).timeline;
%             xs = results.kernel(i).timeseries;
%             ts = [ts fliplr(ts)];
%             xs = [xs-.01 fliplr(xs)+.01] ;
%             cs = results.thetaArray(i,iTheta) * ones(size(ts));
%             
%             l=line(ts,xs(iInput,:),'Color',colors(i,:));
% %             patch(ts,xs(iInput,:),cs,'EdgeColor','flat','EdgeAlpha',.2,'FaceColor','none');
%             hold on; 
%         end
%         xlim([min(ts)-.05,max(ts)+.05]);
%         xlabel('time (s)');
%         ylabel(sprintf('kernel %d',iInput));
%         title(sprintf('connection %s',results.theta.lbl{iTheta}));
%         
%      end
%         
% end

end

    
function  colors=get_colors(betas)
    
    cmap = colormap('jet');
    range = linspace(min(betas),max(betas),size(cmap,1)) ;
    colors = interp1(range,cmap,betas);
    
 
 end
    

% 
% function x=get_xBars(h)
%     f=findall(h,'Type','Patch');
%     xall = get(f,'XData');
%     xall = horzcat(xall{:});
%     xall = xall(:);
%     xall=reshape(xall,4,numel(xall)/4);
%     x = mean(xall);
%     
% end


