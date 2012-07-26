function  [ fx] = f_Qlearn_Pessiglione2006_2Q_M2( x_t,P,u_t,in )
%%% update of a scalar Qvalue given past Qvalue and reward
% INPUT 
% - x_t : Qvalues at trial t-1
% - P : learning rate
% - u_t : action / reward 
a = u_t(1)+1;
fx = x_t; % copy all
fx(2:2:end) = fx(1:2:end); % copy 
f_fname = @f_Qlearn_brick;

        I = [2*(a-1)+1]; % exact indice of Qvalue in hidden states vectors
        [fi] = feval(...
        f_fname,... % the function to be called
        x_t(I),... % the indices of the hidden states concerned by the evolution function 
        P(1),... % the indices of the evolution parameters
        u_t(2),... % the indices of the data
        []);       
        fx(I) = fx(I) + fi; % if multiple updates, their contributions are added
        
        %%%%%%%%%%%%%
% 
% function  [ fx] = f_Qlearn_2Q( x_t,P,u_t,in )
% %%% update of a scalar Qvalue given past Qvalue and reward
% % INPUT 
% % - x_t : Qvalues vector at trial t-1
% % - P : learning rate (through mapping)
% % - u_t : action / reward 
% a = u_t(1)+1; % choice made (1 or 2)
% r = u_t(2); % reward
% ic = 2*(a-1)+1; % index of Qvalue of chosen action for update
% alpha = sigm(P,struct('INV',0)); % learning rate through sigmoid mapping
% 
% fx = x_t; % copy all Qvalues
% fx(2:2:end) = fx(1:2:end); % Copy Qvalues prior to update
% fx(ic) = x_t(ic) + alpha*(r-x_t(ic)); % update 
% 
%        