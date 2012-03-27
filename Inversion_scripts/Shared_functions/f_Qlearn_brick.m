 function  [dq] = f_Qlearn_brick( q,P,r,in )
%%% update of Qvalue a scalar value given past Qvalue and reward
% INPUT 
% - q : Qvalue at trial t-1
% - P : learning rate
% - r : reward 
alpha = sigm(P,struct('INV',0));
dq = alpha*(r-q); % just the update term

