function  [ fx] = f_Qlearn_NQ_sim( x_t,P,u_t,in )
%%% update of a scalar Qvalue given past Qvalue and reward
% INPUT 
% - x_t : [Qpost,Qpre,...,QpostNa,QpreNa]
% - P : learning rate
% - u_t : action / rewards 
 a = u_t(1); % choice made (1 to Na)
 r = u_t(1+a); % reward
 ic = 2*(a-1)+1; % index of Qvalue of chosen action for update
 alpha = sigm(P,struct('INV',0)); % learning rate through sigmoid mapping
 
 fx = x_t; % copy all Qvalues
 fx(2:2:end) = fx(1:2:end); % Copy Qvalues prior to update
 fx(ic) = x_t(ic) + alpha*(r-x_t(ic)); % update 