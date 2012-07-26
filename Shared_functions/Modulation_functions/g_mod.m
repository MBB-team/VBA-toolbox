function [ gx ] = g_mod( x,Pe,u,in )
% Evolution function to manage parametric modulators
% INPUT
% - x : hidden states
% - Pe : parameters of the model containing modulation parameters
% - u : input
% - in : 
%    .f_fname : evolution function handle of the initial model
%    .dim : dimensions of the initial model
%    .phi : cell of info about the modulation of each parameter of the
%    inital model
%       .indpe : location in Pe of the unmodulated parameter
%       .ind_mod : indices of parameter and input modulator
%    .inF

P = zeros(in.dim.n_phi,1); % parameters for the initial model
for i_p = 1:length(in.dim.n_phi) % for each parameter of the initial model
    phi = in.phi{i_p};
    P(phi.indp) = Pe(phi.indpe);
    for i_m = 1:length(phi.indpm)
        P(i_p) =  P(i_p) + Pe(phi.indpm(i_m))*u(phi.indu(i_m));
    end
end

    [gx] = feval(...
        in.g_fname,... % the function to be called for session i
        x,... % the indices of the hidden states concerned by the evolution function for session i
        P,... % the indices of the evolution parameters for session i
        u,... % the indices of the data for session i
        in.inG); 
