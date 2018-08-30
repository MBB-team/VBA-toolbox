function [theta,phi] = defaultHRFparams

% this function outputs efault HRF params, as derived from inverting the
% Balloon model given the canonical SPM HRF.

theta = [
    0.1305
    0.0665
   -0.4951
   -0.1185
         0
    0.4904
    ];

phi = [
    -0.0093
   -0.2772
   ];
