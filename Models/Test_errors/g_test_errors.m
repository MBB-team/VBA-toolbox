
function  [ gx] = g_test_errors( x_t,P,u_t,in )
% observation of unidimensional position of moving item through sigmoid function
% parameterized by an inverse temperature
% INPUT
% - x_t : 2*1 vector, [speed,position]
% - P : Scalar, inverse temperature of the softmax
% - u_t : irrelevant
% OUTPUT
% - gx : Scalar, mapped position
gx = sigm(exp(P(1))*x_t(2));

%----------- TESTING ERRORS
try
    e = in.error;
catch
    e = 0;
end

if e == 1 % out of bounds : parameter
    P(10);
elseif e == 2 % out of bounds : input
    u_t(10);
elseif e == 3 % out of bounds : hidden state
    x_t(10);
elseif e == 4 % erroneous calculation
    [1,1]+[1;1];
elseif e == 5 % output nan
    gx = gx*nan;
elseif e == 6 % wrong dimensions
    gx = [gx;gx];

end



