function h=VBA_susceptibilityAnalysisDisplay2(results,norm,nodes,temperature)


nResps = numel(results.contributions_w);
n_u = numel(results.contributions_u{1});%results.out.dim.u;

idxTheta = results.theta.idx;

cpt = 1;
h=figure('Position',[50,50, 200*n_u, 200*nResps]);


for iObs = 1:nResps
    if nargin < 2 || isempty(norm) || strcmp(norm,'none')
        betas_obs = results.contributions_w{iObs};
        betas_obs=betas_obs/(2*mean(betas_obs(:)));
    elseif strcmp(norm,'output')
        betas_obs = results.contributions_normoutput{iObs};
%         betas_obs = results.contributions_w{iObs};      
%         betas_obs = max(results.contributions_w{iObs},0);      
%         betas_obs = betas_obs./repmat(sum(betas_obs),numel(idxTheta),1);
%         betas_obs = betas_obs/(max(betas_obs(:)));
    elseif strcmp(norm,'param')
%         betas_obs = results.contributions_normparam{iObs}; 
        betas_obs = results.contributions_w{iObs}; 
        betas_obs = max(results.contributions_w{iObs},0);      
        betas_obs = betas_obs./repmat(sum(betas_obs,2),1,n_u);
        betas_obs = betas_obs/(max(betas_obs(:)));
    else
        error('norm should be {''none'',''output'',''param''}');
    end
%     for iu = [3:4 7]
    cmax = max(betas_obs(:));
    if strcmp(norm,'param')
        cmax=1.5*cmax;
    end
    for iu = 1:n_u 
        
        betas = betas_obs(:,iu);
                
%         f=subplot(nResps,n_u,cpt);
        f=subplot('Position',[(iu-1)/n_u (nResps-iObs)/nResps 1/n_u 1/nResps]);
        
        if nargin == 1 
        %% simple bar graph
        bar(betas,'FaceColor',[53 91 135.5]/255); 
        hold on;

        ylabel(sprintf('response %d',iObs));
        xlabel('connection')
        set(gca,'XTickLabel',{results.theta.lbl{:}});
        title(sprintf('kernel u_%d',iu));
        
        else
        %% DCM display
        
        try
    inF = out.options.inF{1};
catch
    inF = out.options.inF;
end
n_nodes = numel(inF.n5);

%% ________________________________________________________________________
%  load DCM matrices
A = inF.A;
B = inF.B;
C = inF.C;
D = inF.D;


first_order = A~= 0 ;
for i=1:numel(B)
    first_order = max(first_order,B{i}~= 0);
end
first_order = first_order - diag(diag(first_order));
for i=1:numel(D)
    second_order{i} = D{i}~= 0;
    second_order{i} = second_order{i} ;% - diag(diag(second_order{i}));
end

mov.structure.first_order = first_order ;
mov.structure.second_order = second_order  ;

%%
        theta = zeros(1,results.out.dim.n_theta);
        theta(idxTheta) = betas;
        theta(theta<0) = 0;
        theta = theta/max(theta);
        boltz = @(x, beta) exp(x/beta) / nansum(exp(x(:)/beta));
%         connect  = 30*grapher_connectivityPattern(results.out,theta,'mean').^2;
        connect  = grapher_connectivityPattern2(results.out,theta);
        connect = connect/max(connect(:));

        connect = boltz(connect,temperature);
        connect = connect/max(connect(:));

        
        
        f=VBA_DCMgrapher_init(nodes,first_order,second_order) ;

        grapher_staticDcmDisplay(nodes,[],connect,f);
        set(f,'CLim',[0 1])
        
        
        
        if strcmp(norm,'param')
            colormap([.99*[1 1 1]; cbrewer('seq','YlOrRd',100).^.7]);
        elseif strcmp(norm,'output')
            cm = [.99*[1 1 1]; (cbrewer('seq','YlGnBu',100)).^1];
            colormap(cm);
        else
            colormap(flipud(colormap('bone')));
        end
%         set(f,'Position',[(iu-1)/n_u (nResps-iObs)/nResps .9/n_u .9/nResps])
%         set(f,'visible','off');
%         colormap(flipud(cbrewer('div','RdYlGn',100)));

        
        end

        cpt = cpt +1;
    end
end

%% pretty colorbar
hc = colorbar();
set(hc,'Location','East');
set(hc,'Ytick',[])
pPlot = get(gca,'Position');
set(hc,'Position',[pPlot(1)+.99*pPlot(3) .3*pPlot(4) .02*pPlot(3) .4*pPlot(4)])

end

    
function  colors=get_colors(betas)
    
    cmap = colormap('jet');
    range = linspace(min(betas),max(betas),size(cmap,1)) ;
    colors = interp1(range,cmap,betas);
    
 
end
 

