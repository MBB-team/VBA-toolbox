% This demo simulates learning agents playing against each other
% This demo considers the following types of learning agents:
% - k-ToM: Bayesian mentalizing agents (cooperative or competitive)
% - Inf: influence learners (cooperative or competitive)
% - k-BSL: Bayesian Sequence Learners
% - RB: Random Biased policy
% - WSLS: win-stay/lose-switch
% - RL: reinforcement learners
% - HGF: Hierarchical Gaussian Filter
% - meta-ToM: mixture of k-ToM and k-BSL experts
% All these only differ wrt their evolution and observation functions.
% Their optional inputs are set automatically using prepare_agent.m.
% The game payoff table can be arbitrary, but two default possibilites are
% considered here: competitive ("hide-and-seek") versus cooperative
% (coordination game). Hereafter, this is referred to as the "interaction
% type". NB: This only makes a difference for mentalizing agents, i.e.
% k-ToM and Inf learners.


clear all
close all
clc

% set the learning styles and interaction types of the first player
styles1 = {'RB','WSLS','RL','0-ToM','Inf','1-ToM','2-ToM','3-ToM','BSL','HGF','metaToM'};
modes1 = {'comp','comp','comp','comp','comp','comp','comp','comp','comp','comp','comp'};

% set the learning styles and interaction types of the second player
styles2 = {'RB','0-ToM','1-ToM','2-ToM'};
modes2 = {'comp','comp','comp','comp'};

% game payof tables (competitive and cooperative interaction types)
compGame = cat(3,[1,-1;-1,1],[-1,1;1,-1]); % competitive game
coopGame = cat(3,[1,-1;-1,1],[1,-1;-1,1]); % cooperative game

nt = 50; % number of trials
noisy = 1; % flag for noisy computations
stdx = exp(-1); % std-dev of computational noise (player #1 only!)

tic

%Nmc = 5000; % # Monte-Carlo simulations
Nmc = 5; 
fprintf( 'WARNING: this is a demo! The number of Monte-Carlo simulation should largely be increased (Nmc~5000) for real simulations\n') 

Perf = zeros(length(styles1),length(styles2),Nmc);
for imc = 1:Nmc
    imc
    for i=1:length(styles1)
        % prepare player 1
        if isequal(modes1{i},'comp')
            info1.payoffTable = compGame;
        else
            info1.payoffTable = coopGame;
        end
        [info1] = prepare_agent(styles1{i},info1.payoffTable,1);
        for j=1:length(styles2)
            % prepare player 2
            if isequal(modes2{j},'comp')
                info2.payoffTable = compGame;
            else
                info2.payoffTable = coopGame;
            end
            [info2] = prepare_agent(styles2{j},info2.payoffTable,2);
%             % sample around evol/obs parameters
%             info1.phi = info1.phi + randn(size(info1.phi));
%             info2.phi = info2.phi + randn(size(info2.phi));
%             info1.theta = info1.theta + randn(size(info1.theta));
%             info2.theta = info2.theta + randn(size(info2.theta));
            % simulate game
            [rew,y1,y2] = runGame_2players(info1,info2,nt,0);
            % summarize game
            Perf(i,j,imc) = mean(rew(1,:)); % player #1's perf
            % volterra decomposition
            [out] = get_VolterraInGames([y1;y2],8,1);
            Aself(i,j,imc) = out.A(1);
            Aother(i,j,imc) = out.A(2);
            % if noisy computations
            if noisy
                try
                    [rew,y1,y2] = runGame_2players(info1,info2,nt,stdx);
                    Perf2(i,j,imc) = mean(rew(1,:)); % (noisy) player #1's perf
                catch
                    Perf2(i,j,imc) = NaN;
                end
            end
        end
    end
    
end

toc

%save ToMgames_ASD2.mat

% plot performance patterns
mP = mean(Perf,3);
sP = std(Perf,[],3);
hf = figure('color',[1 1 1],'name','mean Perf');
ha = axes('parent',hf,'nextplot','add');
imagesc(mP,'parent',ha)
set(ha,'xtick',1:length(styles2),'ytick',1:length(styles1),'xticklabel',styles2,'yticklabel',styles1)
axis(ha,'tight')
axis(ha,'square')
drawnow
rotateXLabels(ha,90);
hf = figure('color',[1 1 1],'name','std Perf');
ha = axes('parent',hf,'nextplot','add');
imagesc(sP,'parent',ha)
set(ha,'xtick',1:length(styles2),'ytick',1:length(styles1),'xticklabel',styles2,'yticklabel',styles1)
axis(ha,'tight')
axis(ha,'square')

hf = figure('color',[1 1 1],'name','Performance patterns');
if noisy
    na = 2;
    mP2 = VBA_nanmean(Perf2,3);
    sP2 = sqrt(VBA_nanvar(Perf2,3));
else
    na = 1;
end
col = 'bgry';
if noisy
    miP = min([mP(:);mP2(:)]);
    maP = max([mP(:);mP2(:)]);
else
    miP = min(mP(:));
    maP = max(mP(:));
end
for i=1:length(styles1)
    ha = subplot(na,length(styles1),i,'parent',hf,'nextplot','add');
    title(ha,[styles1{i},' ',modes1{i}])
    if noisy
        ha2 = subplot(na,length(styles1),length(styles1)+i,'parent',hf,'nextplot','add');
        title(ha2,[styles1{i},' ',modes1{i},' [noisy]'])
    end
    if i==1
        pos = get(ha,'position');
    end
    for j=1:4
        md = mP(i,j);%-mP(i,1);
        stdd = sqrt(sP(i,j).^2);%+sP(i,1).^2);
        bar(j,md,'facecolor',col(j),'parent',ha)
        errorbar(j,md,stdd/sqrt(Nmc),'k.','parent',ha)
        if noisy
            md2 = mP2(i,j);%-mP(i,1);
            stdd2 = sqrt(sP2(i,j).^2);%+sP(i,1).^2);
            bar(j,md2,'facecolor',col(j),'parent',ha2)
            errorbar(j,md2,stdd2/sqrt(Nmc),'k.','parent',ha2)
        end
    end
    set(ha,'ylim',[miP,maP],'xtick',[])
    if noisy
        set(ha2,'ylim',[miP,maP],'xtick',[])
    end
end


return

% plot Volterra analyses
% mAs = mean(Aself,3);
% hf = figure('color',[1 1 1],'name','mean Aself');
% ha = axes('parent',hf,'nextplot','add');
% imagesc(mAs,'parent',ha)
% set(ha,'xtick',1:length(styles2),'ytick',1:length(styles1),'xticklabel',styles2,'yticklabel',styles1)
% axis(ha,'tight')
% axis(ha,'square')
% 
% mAo = mean(Aother,3);
% hf = figure('color',[1 1 1],'name','mean Aother');
% ha = axes('parent',hf,'nextplot','add');
% imagesc(mAo,'parent',ha)
% set(ha,'xtick',1:length(styles2),'ytick',1:length(styles1),'xticklabel',styles2,'yticklabel',styles1)
% axis(ha,'tight')
% axis(ha,'square')

hf = figure('color',[1 1 1],'name','Volterra decomp.');
miAo = min(mAo(:));
maAo = max(mAo(:));
miAs = min(mAs(:));
maAs = max(mAs(:));
sAo = std(Aother,[],3)/sqrt(Nmc);
sAs = std(Aself,[],3)/sqrt(Nmc);
col = 'bgr';
for i=1:length(styles1)
    ha = subplot(2,length(styles1),i,'parent',hf,'nextplot','add');
    title(ha,[styles1{i},' ',modes1{i}])
    ha2 = subplot(2,length(styles1),length(styles1)+i,'parent',hf,'nextplot','add');
    title(ha2,[styles1{i},' ',modes1{i}])
    if i==1
        pos = get(ha,'position');
        pos2 = get(ha2,'position');
    end
    for j=1:3
        bar(j,mAo(i,j),'facecolor',col(j),'parent',ha)
        errorbar(j,mAo(i,j),sAo(i,j),'k.','parent',ha)
        bar(j,mAs(i,j),'facecolor',col(j),'parent',ha2)
        errorbar(j,mAs(i,j),sAs(i,j),'k.','parent',ha2)
    end
    set(ha,'ylim',[miAo,maAo],'xtick',[])
    set(ha2,'ylim',[miAs,maAs],'xtick',[])
end


