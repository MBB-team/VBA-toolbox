function  [ fx] = f_Qlearn_2Q( x_t,P,u_t,in )
%,dfdx,dfdP
%%% update of Qvalue a scalar value given past Qvalue and reward
% INPUT 
% - x_t : Qvalue at trial t-1
% - P : learning rate
% - u_t : reward 
r = u_t(2);
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