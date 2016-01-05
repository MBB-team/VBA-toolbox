function y=bernoulli(p,n)
% BERNOULLI.M
%bernoulli(p,n)
% This function generates n independent draws of a Bernoulli
% random variable with probability of success p. (p should be a column
% vector)
% first, draw n uniform random variables
x = rand(n,length(p));
% set y=1 if x is less than p. This gives probability p of success
y = (x <= p');
% end function definition