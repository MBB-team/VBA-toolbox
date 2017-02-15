% This demo simulates agents playing against each other

clear all
close all
clc

models = {'RB','0-ToM','1-ToM','2-ToM','3-ToM','WSLS','RL','HGF','Inf','BSL'};


% hide-and-seek game
payoffTable = cat(3,[1,-1;-1,1],[-1,1;1,-1]); % game payoff matrix
nt = 50; % number of trials


Nmc = 5000; % # Monte-Carlo simulations
Perf = zeros(length(models),length(models),Nmc);
for imc = 1:Nmc
    imc
   for i=1:length(models)
       % prepare player 1
       [f_p1,g_p1,theta1,phi1,inF1,inG1,x01] = prepare_agent(models{i},payoffTable,1);
       for j=i+1:length(models)
           [f_p2,g_p2,theta2,phi2,inF2,inG2,x02] = prepare_agent(models{j},payoffTable,2);
           % first move
           x1 = [];
           x2 = [];
           x1(:,1) = x01;
           x2(:,1) = x02;
           % sample around evol/obs parameters
           phi1 = phi1 + randn(size(phi1));
           phi2 = phi2 + randn(size(phi2));
           theta1 = theta1 + randn(size(theta1));
           theta2 = theta2 + randn(size(theta2));
           for t=1:nt
               % act
               if t<2
                   u1 = NaN;
                   u2 = NaN;
               else
                   u1 = y2(t-1);
                   u2 = y1(t-1);
               end
               g1(t) = feval(g_p1,x1(:,t),phi1,[0,u1],inG1);
               tmp = VBA_sample('multinomial',struct('p',[g1(t);1-g1(t)],'n',1),1,0);
               y1(t) = tmp(1);
               g2(t) = feval(g_p2,x2(:,t),phi2,[0,u2],inG2);
               tmp = VBA_sample('multinomial',struct('p',[g2(t);1-g2(t)],'n',1),1,0);
               y2(t) = tmp(1);
               % learn
               x1(:,t+1) = feval(f_p1,x1(:,t),theta1,[y2(t),y1(t),u1],inF1);
               x2(:,t+1) = feval(f_p2,x2(:,t),theta2,[y1(t),y2(t),u2],inF2);
               % perf
               rew(1,t) = payoffTable(2-y1(t),2-y2(t),1);
               rew(2,t) = payoffTable(2-y1(t),2-y2(t),2);
           end
           Perf(i,j,imc) = mean(rew(1,:));
           Perf(j,i,imc) = mean(rew(2,:));
       end
   end
    
end

save ToMgames2.mat

mP = mean(Perf,3);
hf = figure('color',[1 1 1],'name','mean Perf');
ha = axes('parent',hf,'nextplot','add');
imagesc(mP,'parent',ha)
plot(ha,[0.5,4.5],[0.5,0.5],'k')
plot(ha,[0.5,4.5],[10.5,10.5],'k')
plot(ha,[0.5,0.5],[0.5,10.5],'k')
plot(ha,[4.5,4.5],[0.5,10.5],'k')
set(ha,'xtick',1:length(models),'ytick',1:length(models),'xticklabel',models,'yticklabel',models)
axis(ha,'tight')
axis(ha,'square')


sP = std(Perf,[],3);
hf = figure('color',[1 1 1],'name','std Perf');
ha = axes('parent',hf,'nextplot','add');
imagesc(sP,'parent',ha)
plot(ha,[0.5,4.5],[0.5,0.5],'k')
plot(ha,[0.5,4.5],[10.5,10.5],'k')
plot(ha,[0.5,0.5],[0.5,10.5],'k')
plot(ha,[4.5,4.5],[0.5,10.5],'k')
set(ha,'xtick',1:length(models),'ytick',1:length(models),'xticklabel',models,'yticklabel',models)
axis(ha,'tight')
axis(ha,'square')
