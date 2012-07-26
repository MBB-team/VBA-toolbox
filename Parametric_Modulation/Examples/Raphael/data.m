clear all
load sources_betas


i_subject = 20;


Y = [];
U = [];


for i_session = 1:6

    y = Betas{i_subject}{i_session}.M1;
    u = [Betas{i_subject}{i_session}.leftTPJ;...
        Betas{i_subject}{i_session}.mpfc;...
        Betas{i_subject}{i_session}.caudate];

    U = [U;u];
    Y = [Y;y];
    
end