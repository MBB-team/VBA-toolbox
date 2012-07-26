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
        P(u_t(3)),... % the indices of the evolution parameters
        u_t(2),... % the indices of the data
        []);       
        fx(I) = fx(I) + fi; % if multiple updates, their contributions are added
        
