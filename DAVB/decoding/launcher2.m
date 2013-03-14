function r=launcher2(N)


links = [1 0; 0 1; 1 1];
decode = [1 0; 0 1];


% r=[];
load('negfeedback_newdat_full_2');
for noise = [1] %  5 25 125 625
    for reps = [5 10 15]
        for li=1:3;
            for di=1:2
            As = 1 * [1 -1] .* links(li,:);
            Bs = 1.5 * [1 -1] .* ([1 1]-links(li,:));
            hAs = decode(di,:);
            r=[r launcher(N,As,Bs,hAs,noise,reps)];
            save('negfeedback_newdat_full_2','r');
            end
        end
    end
end
   
   


end

