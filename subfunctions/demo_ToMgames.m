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
styles1 = {'RB','0-ToM','1-ToM','2-ToM'};
modes1 = {'comp','comp','comp','comp'};

% set the learning styles and interaction types of the second player
styles2 = {'RB','0-ToM','1-ToM'};
modes2 = {'comp','comp','comp'};

% game payof tables (competitive and cooperative interaction types)
compGame = cat(3,[1,-1;-1,1],[-1,1;1,-1]); % competitive game
coopGame = cat(3,[1,-1;-1,1],[1,-1;-1,1]); % cooperative game

nt = 150; % number of trials

tic

% --only for meta-ToM learner--
% hf = figure('name','metaToM belief');
% for i=1:length(styles2)
%     ha(i) = subplot(length(styles2),1,i,'parent',hf,'nextplot','add');
% end

Nmc = 100; % # Monte-Carlo simulations
Perf = zeros(length(styles1),length(styles2),Nmc);
for imc = 1:Nmc
    imc
    for i=1:length(styles1)
        % prepare player 1
        if isequal(modes1{i},'comp')
            payoffTable1 = compGame;
        else
            payoffTable1 = coopGame;
        end
        [f_p1,g_p1,theta1,phi1,inF1,inG1,x01] = prepare_agent(styles1{i},payoffTable1,1);
        for j=1:length(styles2)
            % prepare player 2
            if isequal(modes2{j},'comp')
                payoffTable2 = compGame;
            else
                payoffTable2 = coopGame;
            end
            [f_p2,g_p2,theta2,phi2,inF2,inG2,x02] = prepare_agent(styles2{j},payoffTable2,2);
            % sample around evol/obs parameters
            phi1 = phi1;% + randn(size(phi1));
            phi2 = phi2;% + randn(size(phi2));
            theta1 = theta1;% + randn(size(theta1));
            theta2 = theta2;% + randn(size(theta2));
            % first move
            x1 = [];
            x2 = [];
            x1(:,1) = x01;
            x2(:,1) = x02;
            g1(1) = feval(g_p1,x1(:,1),phi1,NaN(3,1),inG1);
            tmp = VBA_sample('multinomial',struct('p',[g1(1);1-g1(1)],'n',1),1,0);
            y1(1) = tmp(1);
            g2(1) = feval(g_p2,x2(:,1),phi2,NaN(3,1),inG2);
            tmp = VBA_sample('multinomial',struct('p',[g2(1);1-g2(1)],'n',1),1,0);
            y2(1) = tmp(1);
            rew(1,1) = payoffTable2(2-y1(1),2-y2(1),1);
            rew(2,1) = payoffTable2(2-y1(1),2-y2(1),2);
            % next moves
            for t=2:nt
                % build input
                if t==2
                    u1 = [y2(t-1);y1(t-1);NaN];
                    u2 = [y1(t-1);y2(t-1);NaN];
                else
                    u1 = [y2(t-1);y1(t-1);y2(t-2)];
                    u2 = [y1(t-1);y2(t-1);y1(t-2)];
                end
                % learn
                x1(:,t) = feval(f_p1,x1(:,t-1),theta1,u1,inF1);
                x2(:,t) = feval(f_p2,x2(:,t-1),theta2,u2,inF2);
                % act
                g1(t) = feval(g_p1,x1(:,t),phi1,u1,inG1);
                tmp = VBA_sample('multinomial',struct('p',[g1(t);1-g1(t)],'n',1),1,0);
                y1(t) = tmp(1);
                g2(t) = feval(g_p2,x2(:,t),phi2,u2,inG2);
                tmp = VBA_sample('multinomial',struct('p',[g2(t);1-g2(t)],'n',1),1,0);
                y2(t) = tmp(1);
                % perf
                rew(1,t) = payoffTable2(2-y1(t),2-y2(t),1);
                rew(2,t) = payoffTable2(2-y1(t),2-y2(t),2);
            end
            Perf(i,j,imc) = mean(rew(1,:)); % perf of player #1
            % volterra decomposition
            %            [out] = get_VolterraInGames([y1;y2],8,1);
            %            Aself(i,j,imc) = out.A(1);
            %            Aother(i,j,imc) = out.A(2);
            % --only for meta-ToM learner--
            %            Pi(j,:,imc) = sigmoid(x1(1,:));
            %            plot(ha(j),Pi(j,:,imc))
            %            drawnow
        end
    end
    
end

toc

% plot performance patterns
mP = mean(Perf,3);
hf = figure('color',[1 1 1],'name','mean Perf');
ha = axes('parent',hf,'nextplot','add');
imagesc(mP,'parent',ha)
set(ha,'xtick',1:length(styles2),'ytick',1:length(styles1),'xticklabel',styles2,'yticklabel',styles1)
axis(ha,'tight')
axis(ha,'square')
hh = rotateXLabels(ha,90);

sP = std(Perf,[],3);
hf = figure('color',[1 1 1],'name','std Perf');
ha = axes('parent',hf,'nextplot','add');
imagesc(sP,'parent',ha)
set(ha,'xtick',1:length(styles2),'ytick',1:length(styles1),'xticklabel',styles2,'yticklabel',styles1)
axis(ha,'tight')
axis(ha,'square')


return

% plot Volterra analyses
mAs = mean(Aself,3);
hf = figure('color',[1 1 1],'name','mean Aself');
ha = axes('parent',hf,'nextplot','add');
imagesc(mAs,'parent',ha)
set(ha,'xtick',1:length(styles2),'ytick',1:length(styles1),'xticklabel',styles2,'yticklabel',styles1)
axis(ha,'tight')
axis(ha,'square')

mAo = mean(Aother,3);
hf = figure('color',[1 1 1],'name','mean Aother');
ha = axes('parent',hf,'nextplot','add');
imagesc(mAo,'parent',ha)
set(ha,'xtick',1:length(styles2),'ytick',1:length(styles1),'xticklabel',styles2,'yticklabel',styles1)
axis(ha,'tight')
axis(ha,'square')

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


