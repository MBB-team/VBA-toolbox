function  [ gx ] = g_softmax_2Q_modu( x_t,P,u_t,in )
%%% Softmax emission law for RL with 2 Qvalues
% INPUT
% - x_t : 2*1 vector, [Q1;Q2]
% - P : Scalar, inverse temperature of the softmax decision rule
% - u_t : 2*1 vector, Action and reward
% OUTPUT
% - gx : Scalar, P(a=a1|x_t), probability of action 1


% Transformation on parameters
P(1) =  P(1) + P(2)*u_t(3) + P(3)*u_t(4) ; % modulation of beta by input 1,2,3

try 
    if isequal(in.param_transform.type,'exponential')
        beta = exp(P(1));
    elseif isequal(in.param_transform.type,'modified sigmoid')
        beta = in.param_transform.a + (in.param_transform.b-in.param_transform.a)./(1+exp(-P(1)));                
    end
catch
    beta = P(1);
end



dQ = (x_t(1)-x_t(2));
gx =sig( -beta*dQ );


function y=sig(x)
y = 1/(1+exp(-x));
y(y<eps) = eps;
y(y>1-eps) = 1-eps; 
