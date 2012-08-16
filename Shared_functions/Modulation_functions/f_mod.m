function [ fx ] = f_mod( x,Pe,u,in )
% Evolution function to manage parametric modulators
% INPUT
% - x : hidden states
% - Pe : parameters of the modulated model 
% - u : input
% - in : 
%    .f_fname : evolution function handle of the initial model
%    .dim : dimensions of the initial model
%    .theta : cell of info about the modulation of each parameter of the
%    inital model
%       .indpi : location in Pe of the unmodulated parameter
%       .indpm : indices of parameter and input modulator
%    .inF

P = zeros(in.dim.n_theta,1); % parameters for the initial model
for i_p = 1:in.dim.n_theta % for each parameter of the initial model
    theta = in.theta{i_p};
    P(theta.indp) = Pe(theta.indpi);
    for i_m = 1:length(theta.indpm)
        P(theta.indp) = P(theta.indp) + Pe(theta.indpm(i_m))*u(theta.indu(i_m));
     %   P(i_p) =  P(i_p) + Pe(theta.indpm(i_m))*u(theta.indu(i_m));
    end
end

fx = zeros(size(x));
    [fi] = feval(...
        in.f_fname,... % the function to be called for session i
        x,... % the indices of the hidden states concerned by the evolution function for session i
        P,... % the indices of the evolution parameters for session i
        u,... % the indices of the data for session i
        in.inF); 
fx = fi; 


