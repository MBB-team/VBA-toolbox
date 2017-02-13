function  y = invsigmoid(x)
% inverse-sigmoid mapping [see sigmoid.m].
% function  y = invsigmoid(x)
% The inverse sigmoid mapping is defined as follows: y = log(x/(1-x)).
y=log(x./(1-x));
% deal with numerical issues
y(x<=1e-4)=-9.2102;
y(x>=(1-1e-4))=9.2102;