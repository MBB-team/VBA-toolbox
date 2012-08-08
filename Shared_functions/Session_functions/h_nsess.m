function [fb] = h_nsess(y,t,in)
%fb = zeros(in.dim_fb,1);
fb = zeros(size(in.indfb));
for i=1:in.nsess
    [fbi] = feval(...
        in.sess(i).h_fname,... % the function to be called for session i
        y(in.sess(i).indy),... % the indices of output concerned
        t,...
        in.sess(i).inH); % the indices of the data
 %   fb( ((i-1)*in.dim_singlefb+1) : i*in.dim_singlefb )  = fbi; 
    fb( in.sess(i).indfb )  = fbi; 
end
