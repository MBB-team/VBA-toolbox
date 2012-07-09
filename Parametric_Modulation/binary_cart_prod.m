
function o = binary_cart_prod(N_p)
 o = zeros(2^N_p,N_p);
 for i = 0:2^N_p-1
    u = str2num(dec2bin(i)')';
    o(i+1,end-length(u)+1:end) =  u;
 end
end