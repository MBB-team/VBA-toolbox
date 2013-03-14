function la = bp_la(u, p, rf)

% Unpack parameters
mu2_0 = p(1);
sa2_0 = p(2);
mu3_0 = p(3);
sa3_0 = p(4);
ka    = p(5);
om    = p(6);
th    = p(7);

% Regularization count
rc = 0;

% Regularization list
rl = {};

% Add dummy "zeroth" trial
u = [0, u];
n = length(u);

% Initialize updated quantities

% Representations
mu1 = NaN(1,n);
mu2 = NaN(1,n);
pi2 = NaN(1,n);
mu3 = NaN(1,n);
pi3 = NaN(1,n);

% Other quantities
pi1hat = NaN(1,n);
pi2hat = NaN(1,n);
pi3hat = NaN(1,n);
w2     = NaN(1,n);
da1    = NaN(1,n);
da2    = NaN(1,n);

% Representation priors
% Note: first entries of the other quantities remain
% NaN because they are undefined and are thrown away
% at the end; their presence simply leads to consistent
% trial indices.
mu1(1) = sgm(mu2_0, 1);
mu2(1) = mu2_0;
pi2(1) = 1/sa2_0;
mu3(1) = mu3_0;
pi3(1) = 1/sa3_0;

% Pass through representation update loop
for k = 2:1:n

  %%%%%%%%%%%%%%%%%%%%%%
  % Effect of input u(k)
  %%%%%%%%%%%%%%%%%%%%%%
  
  % 1st level
  % ~~~~~~~~~
  % Precision of prediction
  pi1hat(k) = 1/(sgm(mu2(k-1), 1)*(1 -sgm(mu2(k-1), 1)));
  
  % Update
  mu1(k) = u(k);
  
  % Prediction error
  da1(k) = mu1(k) -sgm(mu2(k-1), 1);
  
  % 2nd level
  % ~~~~~~~~~
  % Precision of prediction
  pi2hat(k) = 1/(1/pi2(k-1) +exp(ka *mu3(k-1) +om));
  
  % Updates
  pi2(k) = pi2hat(k) +1/pi1hat(k);
  
  mu2(k) = mu2(k-1) +1/pi2(k) *da1(k);
  
  % Volatility prediction error
  da2(k) = (1/pi2(k) +(mu2(k) -mu2(k-1))^2) *pi2hat(k) -1;
  
  
  % 3rd level
  % ~~~~~~~~~
  % Precision of prediction
  pi3hat(k) = 1/(1/pi3(k-1) +th);
  
  % Weighting factor
  w2(k) = exp(ka *mu3(k-1) +om) *pi2hat(k);
  
  % Updates
  pi3(k) = pi3hat(k) +1/2 *ka^2 *w2(k) *(w2(k) +(2 *w2(k) -1) *da2(k));
  
  if pi3(k) <= 0 || isnan(pi3(k))
      pi3(k) = rf*pi3(k-1);
      rc = rc+1;
      rl{end+1} = ['pixh(', num2str(k), ')'];
  end

  mu3(k) = mu3(k-1) +1/2 *1/pi3(k) *ka *w2(k) *da2(k);

end

% Remove representation priors
mu1(1)  = [];
mu2(1)  = [];
pi2(1)  = [];
mu3(1)  = [];
pi3(1)  = [];

% Remove other dummy initial values
pi1hat(1) = [];
pi2hat(1) = [];
pi3hat(1) = [];
w2(1)     = [];
da1(1)    = [];
da2(1)    = [];

% Create result data structure
la = struct;

la.mu1     = mu1;
la.mu2     = mu2;
la.sa2     = 1./pi2;
la.mu3     = mu3;
la.sa3     = 1./pi3;

la.sa1hat  = 1./pi1hat;
la.sa2hat  = 1./pi2hat;
la.sa3hat  = 1./pi3hat;
la.w2      = w2;
la.da1     = da1;
la.da2     = da2;

la.rc      = rc;
la.rl      = rl;

return;
