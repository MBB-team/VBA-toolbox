function [ gx ] = ObsRecGen(x,P,u,in)
% Marie devaine wrote this
% ObsRecGen obsevration function general for recursive ToM
%   This one is the simplest one and does not include bias
% x are the hidden states
% P is the temperature
% u is useless here can be left empty
% in contains important information such as the level of the player, its
% role and the game nb total de para
player= in.player; % 1 or 2: role of the player
ntotPara= in.npara;
lev=in.lev;
game=in.game;
alpha=0.36;
indlev=in.indlev;

if lev==0
    muT=x(1);
    SigT=exp(x(2));
    Pi= sigmoid( muT/(sqrt(1+alpha*SigT)));
else
    vecp=sigmoid(x(1:(lev-1)));
    vecp(end+1)=max(0,1-sum(vecp));
    vecV=zeros(lev,1);
    vecPi=zeros(lev,1);
    for j=1:lev
        f=x(indlev(j).f);
        df=x(indlev(j).df)';
        
        for i=1:ntotPara
            Sig=exp(x(indlev(j).Par(2*i)));
            vecV(j)=vecV(j)+Sig*df(i)^2;
        end
        vecPi(j)=sigmoid(f/sqrt(1+alpha*vecV(j)));
 
    end
    Pi= vecp*vecPi;
end
if length(P)==1
gx=sigmoid(fplayer(Pi,exp(P(1)),player,game));
else %bias
  gx=sigmoid(fplayer(Pi,exp(P(1)),player,game)+P(2));
end
end
