function  [ gx ] = g_softmax_NQ( x_t,P,u_t,in )
%%% Softmax decision for Na alternatives
% INPUT
% - x_t : (Na*2)*1 vector, [Qpost1;Qpre1;....QpostNa;QpreNa]
% - P : Scalar, inverse temperature of the softmax decision rule
% - u_t : (1+Na)*1 vector, [1]Index of selected action and [2,3]rewards for both alternatives
% OUTPUT
% - gx : Vector, probability of all actions
beta = exp(P(1)); %  inverse temperature, through exp mapping
v = exp(beta*x_t(2:2:end)); %Proba on copy
gx=v./sum(v);
%gx = v(u_t(1)); % Proba of performed action
