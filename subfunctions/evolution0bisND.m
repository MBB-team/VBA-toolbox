function [staT]=evolution0bisND(x,theta, u,in)
%evolution function for 0-ToM player without doubled hidden states
%theta is log(sigma) is a parameter of the player, mu and Sig are the sufficient 
%statistics of last trial to be updated
% u is the new observation (last move of other player
% to avoid numerical problems Sig is first put in the log space and then
% passed back through an exp
% theta is also transform through 1.4* sigmoid() so that no explosion is
% created when simulated by higher ToM-levels
if isempty(u)||isnan(u(1)) % signals a missed trial
    staT=x;
else    
    mu=x(1);
    LSig=x(2); %sigma in log space
    Sig=exp(LSig);
    SigT=1/(1/(1.4*sigmoid(theta)+Sig)+sigmoid(mu)*(1-sigmoid(mu)));
%    SigT=1/(1/(exp(theta)+Sig)+sigmoid(mu)*(1-sigmoid(mu))+1e-2);
    muT=mu+SigT*(u(1)-sigmoid(mu));
    LSigT=log(SigT);
    
    staT=[muT;LSigT];%mu and LSig are the copies
end    