function sensitivityDisplay(results)

% close all;

%== parameters of interest extract


%

[n_x, n_theta]=size(results.perturb);
if nargin<2
    thetas=1:n_theta;
else
    n_theta=length(thetas);
end

n_u = size(results.effect{1});
lims_u = [-.05*n_u(:) 1.05*n_u(:)];

figure('Position',[-1000 (900-75*(n_theta+1)) n_x*150 (n_theta+1)*150]);
cpt=1;

% effects
for i=1:n_x;
    
        subplot(n_theta+1,n_x,cpt);
        
        if size(results.contrast,1) ==2
        
        h=mesh(results.predict{i});
        hold on

        [h1 h2]=plotEffect(results.effect{i},results.SigmaEffect{i}) ;
        xlabel('u1');
        ylabel('u2');
        pretty(h,[.8 0 0])
        pretty2(h1, h2,[.5 0 0])
        ze{i} = get(gca,'ZLim');
        
        else
            
            plot(results.effect{i},'.','Color',[.8 0 0]);
            hold on
            h=plot(results.predict{i});
            pretty(h,[.8 0 0]);

        end
        
        cpt=cpt+1;
        
        title(sprintf('\\fontsize{10}\\bf y=%d',i));
     

end

% results.perturb=cellfun(@(x) x/trace(x),results.perturb,'UniformOutput',false);

% contributions


for j=thetas
    
for i=1:n_x;
    
    z(1)= min(min(cellfun(@(x) min(x(:)),results.perturb(i,:)))) ; % normilize scale accross connections for each observation
    z(2)= max(max(cellfun(@(x) max(x(:)),results.perturb(i,:)))) ;
    cs=results.w{i};
    
        subplot(n_theta+1,n_x,cpt);
        
        if size(results.contrast,1) ==2
        h0=mesh(zeros(size(results.perturb{i,j})));
        hold on
        set(h0,'LineStyle',':','FaceColor','none','EdgeColor',[.1 .1 .1]);
        h=mesh(results.perturb{i,j});
        xlabel('u1');
        ylabel('u2');
         pretty(h,abs(cs(j)),z)     
        
        if i==1
            zlabel(sprintf('\\fontsize{12} \\bf %s',results.labels{j}));
        end
        else
           h0=plot(zeros(size(results.perturb{i,j})),':','Color',[.1 .1 .1]);
           hold on
           
           h=plot(results.perturb{i,j}) ;
           xlabel('c1');
           pretty(h,abs(cs(j)),z);
            if i==1
            ylabel(sprintf('\\fontsize{12} \\bf %s',results.labels{j}));
            end
        end
        
        cpt=cpt+1;
        
    end
end




%=
hc=colorbar();
set(hc,'Units','normalized', 'position', [0.95 0.11 0.01 1/5]);
set(hc,'YTickMode','manual','YTick',[]);

end

function [h h2]=plotEffect(eff,sig)
        h=mesh(eff);
        h2={};        
%         xs=get(h,'XData');
%         ys=get(h,'YData');

%         for xi=1:length(xs)
%             for yi=1:length(ys)
%                 z1=eff(xi,yi)-sqrt(sig(xi,yi));
%                 z2=eff(xi,yi)+sqrt(sig(xi,yi));
%                 h2{end+1}=plot3([ys(yi) ys(yi)],[xs(xi) xs(xi)],[z1 z2]);
%             end
%         end
end

function pretty(h,c,z)

set(gca,'XGrid','off')
set(gca,'XTick',[])
set(gca,'YGrid','off')
set(gca,'YTick',[])

if isempty(get(h,'ZData'))
    try
        z_ =   max(abs(z));
        ylim([-z_ z_]);
        set(gca,'YTick',[0])
        set(gca,'YGrid','off')
    catch
        set(gca,'YTick',get(gca,'YLim'));
    end
    if length(c)==1
        cm=flipud(colormap('Bone'));
        ci=round(10+(c*(size(cm,1)-10)));
        colormap(cm);
        set(h,'Color',cm(ci,:));
        caxis([0 1])
    else
        set(h,'Color',c);
    end
    set(h,'LineWidth',2);
    
else
    set(gca,'Box','on')
    set(h,'EdgeColor',[0 0 0]);
    set(h,'FaceAlpha',1)
    set(gca,'ZGrid','off')

    if length(c)==1
        set(h,'FaceColor','flat');
        CD = get(h,'CData');
        CD(:) = c;
        set(h,'CData',CD);
        colormap(flipud(colormap('Bone')));
        caxis([0 1])
    else
        set(h,'FaceColor',c);
    end
    try
        z_ =   max(abs(z));
        zlim([-z_ z_]);
        set(gca,'ZTick',[0])
        set(gca,'ZGrid','off')
    catch
        set(gca,'ZTick',get(gca,'ZLim'));    
    end
end






% margin
xl = get(gca,'XLim');
dx=.05*(xl(2)-xl(1));
set(gca,'XLim',[xl(1)-dx xl(2)+dx]);
yl = get(gca,'YLim');
dy=.05*(yl(2)-yl(1));
set(gca,'YLim',[yl(1)-dy yl(2)+dy]);

end

function pretty2(h,h2,c)

    set(h,'EdgeColor',c);
    set(h,'FaceColor','none');
    set(h,'Marker','.');
    set(h,'MarkerSize',15);
    set(h,'LineWidth',.5);
    for i=1:length(h2)
        set(h2{i},'Color',c);
        set(h2{i},'LineWidth',2);
    end
    
end