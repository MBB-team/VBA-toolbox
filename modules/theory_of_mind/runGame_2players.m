function [rew,y1,y2] = runGame_2players(info1,info2,nt,stdx)
% simulate an interative game btw 2 players
% IN:
%   - info1/info2: structures containing the following fields:
%       .f_p: handle of evolution function
%       .g_p: handle of observation function
%       .theta: evolution params
%       .phi: observation params
%       .x0: initial conditions
%       .inF: optional input to the evolution function
%       .inG: optional input to the observation function
%   - nt: number of game trials
%   - stdx:  % std-dev of computational noise (player #1 only!)
% OUT:
%   - rew: 2xnt matrix of players' action outcomes
%   - y1/y2: 1xnt vectors of players' actions

% player 1
f_p1 = info1.f_p;
g_p1 = info1.g_p;
theta1 = info1.theta; 
phi1 = info1.phi;
x01 = info1.x0;
inG1 = info1.inG;
inF1 = info1.inF;
payoffTable1 = info1.payoffTable;

% player 2
f_p2 = info2.f_p;
g_p2 = info2.g_p;
theta2 = info2.theta; 
phi2 = info2.phi;
x02 = info2.x0;
inG2 = info2.inG;
inF2 = info2.inF;
payoffTable2 = info1.payoffTable;

% first move
x1 = [];
x2 = [];
x1(:,1) = x01 + stdx.*randn(size(x01));
x2(:,1) = x02;
g1(1) = g_p1(x1(:,1),phi1,NaN(3,1),inG1);
y1(1) = VBA_random ('Bernoulli', g1(1));
g2(1) = g_p2(x2(:,1),phi2,NaN(3,1),inG2);
y2(1) = VBA_random ('Bernoulli', g2(1));

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
    x1(:,t) = f_p1(x1(:,t-1),theta1,u1,inF1);
    x2(:,t) = f_p2(x2(:,t-1),theta2,u2,inF2);
    % perturb learning
    x1(:,t) = x1(:,t) + stdx.*randn(size(x1(:,t)));
    % act
    g1(t) = g_p1(x1(:,t),phi1,u1,inG1);
    y1(t) = VBA_random ('Bernoulli', g1(t));

    g2(t) = g_p2(x2(:,t),phi2,u2,inG2);
    y2(t) = VBA_random ('Bernoulli', g2(t));

    % perf
    rew(1,t) = payoffTable1(2-y1(t),2-y2(t),1);
    rew(2,t) = payoffTable2(2-y1(t),2-y2(t),2);
end


