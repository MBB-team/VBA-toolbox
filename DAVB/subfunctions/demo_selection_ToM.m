% demo for evolutionary selection of ToM sophistication

% ToM sophistication is a phenotype or trait that is defined in terms of a
% "ToM level", which scores the depth of recursive inference on other's
% beliefs and preferences given their observed behaviour. This recursion is
% induced by reciprocal interations, which are well captured in the context
% of iterated games. Level-0 players adapt their behaviour without caring
% for other's intentions or beliefs, by tracking the evolving likelihood of
% their opponent's choices. Level-1 players assume thay face level-0
% players, and attempt to recognize the priors that shape the way they
% learn, and hence their behavioural tendencies. Level-k players (where
% k>1) have to idenfity the level of their opponent, in addition to her
% priors. Theoretically, one can show that there is no limit to the ToM
% recursion or sophistication level.
% Let x_k be the frequency of level-k players within the population. This
% demo uses standard replicator dynamics of evolutionary game theory to
% derive the adaptive fitness of the trait, ie:
%
% dx_k/dt = x_k ( f_k(x) + sum_j{x_j*f_j(x)} )
%
% where f_k(x) is the malthusian fitness of ToM level-k agents.
% This fitness is obtained according to the expected outcome of games
% (here: a mixture of "battle of the sexes" and "hide and seek"), for each
% type of agent or trait. In our case, the ToM level changes the learning
% rule of agents between successive trials of the iterated game.
% Let lambda be the probability of interacting under the rules of the game
% "hide and seek", and tau the game duration. Each couple of agents (ie
% each pair of ToM levels) can be characterized in terms of the expected
% outcome of both players after tau iterations and for a given lambda.This
% induces a outcome matrix Q(lambda,tau), the i^th line of which gives the
% expected reward delivered to a level-i player, when playing against each
% ToM level.
% Note: the malthusian fitness of a ToM level-k agent is its payoff,
% averaged across the population. The chance of encountering any member of
% type k' is proportional to its frequency wihin the population. This
% yields:
%
% f(x) = Q(lambda,tau)*x
%
% where the competition is really against learning strategies or traits.
% This means the rate of change of trait's frequencies within the
% population is a function of both the expected payoff (given expected
% actions of each type of agents) and the current distribution of trait
% frequencies.




% 0-first fit analytical functions to expected outcome tables
games = {BSScore,HSScore,SHScore};
games_names = {'battle of the sexes','hide and seek','stag-hunt'};
mo = [0.5625,0.5,1]; % mean outcome for each game

% regularize MCMC expected payoffs?
regularize = 1;

m = 7; % max ToM level
if regularize  % max number of iterations in repeated games
    T = 2^11;
else
    T = 512;
end
ng = 3; % number of games



try
    Q;
catch
    
    clear all
    close all
    clc
    
    % load expected outcome tables
    load('HSScorew4_RL_Nash')
    
    
    Q = zeros(m,m,T,ng);
    vQ = zeros(m,m,512,ng);
    extrema = zeros(ng,2);
    
    for i=1:ng
        hf(i) = figure('color',[1 1 1],'name',games_names{i});
        for k=1:m
            for l=1:m
                y0 = squeeze(games{i}.P1(k).P2(l).par1(1).par2(1).mean -mo(i).*[1:512]);
                vy0 = squeeze(games{i}.P1(k).P2(l).par1(1).par2(1).sd.^2);
                y0 = y0(:);
                vy0 = vy0(:);
                if regularize
                    if i==2 && k==l
                        gx = zeros(T,1);
                    else
                        u = [1:512]'-1;
                        disp('-----')
                        disp([games_names{i},': level ',num2str(k-1),' against level ',num2str(l-1)])
                        [gx] = fitFixedForm(y0,u,[1:T]'-1,0);
                    end
                else
                    gx = y0;
                end
                
                Q(k,l,:,i) = gx./[1:T]';
                vQ(k,l,:,i) = vy0./[1:512]';
                
                extrema(i,1) = min([extrema(i,1);gx]);
                extrema(i,2) = max([extrema(i,2);gx]);
                ha(k,l,i) = subplot(m,m,k+(l-1)*m,...
                    'parent',hf(i),'nextplot','add');
                title(ha(k,l,i),['level ',num2str(k-1),' against level ',num2str(l-1)])
                plot(ha(k,l,i),1:T,gx,'k')
                plotUncertainTimeSeries(y0',vy0',1:512,ha(k,l,i));
                set(ha(k,l,i),'xlim',[1 T])
                getSubplots
                drawnow
            end
        end
        set(vec(ha(:,:,i)),'ylim',extrema(i,:))
    end
    
end

% return

legstr = {'k=0','k=1','k=2','k=3','k=4','RL','Nash'};

m = 5;

% Parameters of the simulation
n_t = 2^10;%32e2;
dt = 2e-0;
f_fname = @log_replicator;
g_fname = @g_odds;
alpha   = Inf;
sigma   = Inf;
phi     = [];
theta   = [0.5,0]; % P(1): lambda (frequency of coop game), P(2): tau
x0      = zeros(m,1); % log-odds for each ToM level
u = [];

in.Q = Q(1:m,1:m,:,:);
in.tau = [1:512];
in.coop_game = 3; % 1: "battle of sexes" , 3: "stag-hunt"
in.dt = dt;
in.f_fitness = @fitness_ToM;

% options.checkGrads = 1;
options.inF         = in;
options.inG         = in;
options.verbose     = 0;
dim.n_theta         = 2;
dim.n_phi           = 0;
dim.n               = m;
dim.n_t             = n_t;
dim.p             = dim.n;
% fill in options structure (including default priors)
[options,u,dim] = VBA_check(sparse(dim.n,n_t),[],f_fname,g_fname,dim,options);


gcoop = [1,3];
ngc = length(gcoop);
if regularize
    gtau = 2.^[0:1:11];
else
    gtau = 2.^[0:1:9];
end
gtau = 16;%2^9;
ngt = length(gtau);
glambda = 0.5;%0:0.1:1;
ngl = length(glambda);

XX = zeros(m,ngt*ngl,ngc);
gridX = zeros(1,ngt*ngl,ngc);
gridY = zeros(1,ngt*ngl,ngc);
Xeq = cell(ngt,ngl,ngc);
feq = cell(ngt,ngl,ngc);

N = 128;

% Loop over game durations

disp(' ')
disp('----- Looping over game duration and mixing frequencies -----')
disp(' ')

dbstop if error

for k=1%:2
    
    k
    tmp = gcoop(k)
    options.inF.coop_game = gcoop(k);
    options.inG.coop_game = gcoop(k);
    
    for j=1:ngl
        
        X = zeros(m,ngt);
        
        for i=1:ngt
            
            disp(['lambda= ',num2str(glambda(j)),', tau= ',num2str(gtau(i))])
            
            % choose coop game frequency and duration
            theta = [glambda(j),gtau(i)];
            
            % find multiple equilibria
            [eq,tmp,out,ha,ha2] = findEquilibria(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,dim,N,1);
            legend(ha2,legstr)
            Xeq{i,j,k} = eq;
            feq{i,j,k} = out.lambda;
            
            % simulate replicator dynamics from maxEnt distribution
            [y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);
            
            % store ESS
            X(:,i) = y(:,end);
            
            XX(:,i+(j-1)*ngt,k) = X(:,i);
            gridX(i+(j-1)*ngt,k) = gtau(i);
            gridY(i+(j-1)*ngt,k) = glambda(j);
            
            
            %         % evaluate stability of steady state
            %         J = numericDiff(@f_replicator,1,y(:,end),theta,[],in);
            %         ev = eig((J-eye(dim.n))./dt);
            %         hf = figure('color',[1 1 1],'name','ESS: stability analysis');
            %         ha = subplot(2,1,1,'parent',hf,'nextplot','add');
            %         rectangle('Position',[-1,-1,2,2],'curvature',[1 1],'parent',ha)
            %         for i=1:length(ev)
            %             plot(ha,real(ev(i)),imag(ev(i)),'r+')
            %         end
            %         axis(ha,'equal')
            %         grid(ha,'on')
            %         ha = subplot(2,1,2,'parent',hf,'nextplot','add');
            %         plot(ha,y')
            %         pause
            
        end
        
        %         if ngt > 1
        %             hf = figure('color',[1 1 1],'name',['lambda=',num2str(glambda(j))]);
        %             h(1) = subplot(2,1,1,'parent',hf,'nextplot','add');
        %             plot3(h(1),X(2,:),X(3,:),X(4,:),'k')
        %             plot3(h(1),X(2,:),X(3,:),X(4,:),'k.')
        %             plot3(h(1),X(2,end),X(3,end),X(4,end),'ro')
        %             grid on
        %             xlabel(h(1),'k=1')
        %             ylabel(h(1),'k=2')
        %             zlabel(h(1),'k=3')
        %             h(2) = subplot(2,1,2,'parent',hf,'nextplot','add');
        %             plot(h(2),gtau,X')
        %             legend(h(2),{'k=0','k=1','k=2','k=3'})
        %         end
        
    end
    
    dx = mean(diff(log(gtau)./log(2)));
    dy = mean(diff(glambda));
    
    if ngt > 1 %&& ngl > 1
        hf = figure(...
            'color',[1 1 1],...
            'name',[games_names{gcoop(k)},' : susceptibility to lambda and tau']);
        ha = axes('parent',hf,'nextplot','add');
        multipie(XX(:,:,k)+eps,dy.*log(gridX(:,:,1))./log(2),dx.*gridY(:,:,1),ha);
        legend(ha,legstr(1:m))
        xlabel(ha,'tau (game duration)')
        ylabel(ha,'lambda (frequency of cooperative game)')
        set(ha,...
            'xtick',dy.*log(gtau)./log(2),...
            'xticklabel',round(gtau),...
            'ytick',glambda,'yticklabel',glambda)
    end
    
    
    
%     save results_selection_ToM4.mat
    
    
    
end




