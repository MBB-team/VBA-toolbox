% Polya's urn

close all

% Fig 1

gridn =2.^[1:10];
nn = length(gridn);

OBF = zeros(nn,1);
EP = zeros(nn,1);

a0 = 10;

options.binomial = 1;
dim.n_phi = 1;
dim.n = 0;
dim.n_theta = 0;
n=1024;
y = [ones(1,n),zeros(1,n)];
g_fname = @g_Id_phi;
[posterior,out] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);
% pause
% 
% for i=1:nn
%     n = gridn(i);
%     nk = [n+1;n];
%     OBF(i) = PolyaOBF(a0,nk);
%     EP(i) = VB_PPM(nk(1)+a0,nk(2)+a0,0.5,0,'beta');
% %     pause
% end
% figure,plot(log(OBF),'.')
% figure,plot(EP,'.')  