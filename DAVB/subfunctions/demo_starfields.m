% this demo simulates the response of the visual cortex to sparse visual
% inputs. Here, the cortex is hierarchically organized into three layers,
% with increasing receptive fields sizes. Critically, the visual inputs are
% such that the first level receptive field are too small to integrate the
% information spatially.

clear variables
% close all

n = 8;

in.A1 = kron(eye(n),[[ones(4,1);zeros(4,1)],[zeros(4,1);ones(4,1)]]);
in.A2 = kron(eye(n),ones(2,1));
in.mu3_0 = zeros(n,1);
in.df = 1e-2;

[n1,n2] = size(in.A1);
n3 = size(in.A2,2);

smooth = 1
if smooth
    in.L1 = lap(n1)'*lap(n1);
    in.L2 = lap(n2)'*lap(n2);
    in.L3 = lap(n3)'*lap(n3);
else 
    in.L1 = eye(n1);
    in.L2 = eye(n2);
    in.L3 = eye(n3);
end

s3 = [.1,1,10];

N = 2^8;

dispUpdates = 0;

m1 = zeros(length(s3),N);
m2 = zeros(length(s3),N);
m3 = zeros(length(s3),N);
pe1 = zeros(length(s3),N);
pe2 = zeros(length(s3),N);
pe3 = zeros(length(s3),N);
pe4 = zeros(length(s3),N);

for j=1:length(s3)
    
    P = .1*[s3(j);1;1;1];
    
    
    for i=1:N
        
        
        % simulate data under visual cortex model
        u = visualU(P,in);
        
        % let the visual cortex do its VB update
        x = visualVB([],P,u,in);
        
        
        
        i0 = 0;
        mu1 = x(i0+1:i0+n1,:);
        i0 = i0+n1;
        mu2 = x(i0+1:i0+n2,:);
        i0 = i0+n2;
        mu3 = x(i0+1:i0+n3,:);
        i0 = i0+n3;
        
        PE1 = log(x(i0+1:i0+n1,:).^2+eps);
        i0 = i0+n1;
        PE2 = log(x(i0+1:i0+n1,:).^2+eps);
        i0 = i0+n1;
        PE3 = log(x(i0+1:i0+n2,:).^2+eps);
        i0 = i0+n2;
        PE4 = log(x(end,:).^2+eps);
        
        
        % display update?
        if dispUpdates
            hf = figure('color',[1 1 1]);
            ha = subplot(4,2,1,'parent',hf);
            plot(ha,mu1')
            title(ha,'mu1')
            ha = subplot(4,2,2,'parent',hf);
            plot(ha,PE1')
            title(ha,'PE1 = u-mu1')
            ha = subplot(4,2,3,'parent',hf);
            plot(ha,mu2')
            title(ha,'mu2')
            ha = subplot(4,2,4,'parent',hf);
            plot(ha,PE2')
            title(ha,'PE2 = mu1-A1*mu2')
            ha = subplot(4,2,5,'parent',hf);
            plot(ha,mu3')
            title(ha,'mu3')
            ha = subplot(4,2,6,'parent',hf);
            plot(ha,PE3')
            title(ha,'PE3 = mu2-A2*mu3')
            ha = subplot(4,2,8,'parent',hf);
            plot(ha,PE4')
            title(ha,'PE4 = mu3-mu3_0;')
            pause
        end

        
        m1(j,i) = sum(mu1(:,end));
        m2(j,i) = sum(mu2(:,end));
        m3(j,i) = sum(mu3(:,end));
        
        pe1(j,i) = sum(PE1(:,end));
        pe2(j,i) = sum(PE2(:,end));
        pe3(j,i) = sum(PE3(:,end));
        pe4(j,i) = sum(PE4(:,end));
        
    end
    
end

hf = figure('color',[1 1 1]);
ha = subplot(4,2,1,'parent',hf);
plotUncertainTimeSeries(mean(m1,2),var(m1,[],2)/N,[],ha);
title(ha,'mu1')
ha = subplot(4,2,2,'parent',hf);
plotUncertainTimeSeries(mean(pe1,2),var(pe1,[],2)/N,[],ha);
title(ha,'PE1 = u-mu1')
ha = subplot(4,2,3,'parent',hf);
plotUncertainTimeSeries(mean(m2,2),var(m2,[],2)/N,[],ha);
title(ha,'mu2')
ha = subplot(4,2,4,'parent',hf);
plotUncertainTimeSeries(mean(pe2,2),var(pe2,[],2)/N,[],ha);
title(ha,'PE2 = mu1-A1*mu2')
ha = subplot(4,2,5,'parent',hf);
plotUncertainTimeSeries(mean(m3,2),var(m3,[],2)/N,[],ha);
title(ha,'mu3')
ha = subplot(4,2,6,'parent',hf);
plotUncertainTimeSeries(mean(pe3,2),var(pe3,[],2)/N,[],ha);
title(ha,'PE3 = mu2-A2*mu3')
ha = subplot(4,2,8,'parent',hf);
plotUncertainTimeSeries(mean(pe4,2),var(pe4,[],2)/N,[],ha);
title(ha,'PE4 = mu3-mu3_0')


hf = figure('color',[1 1 1]);
ha = subplot(4,2,1,'parent',hf);
plot(ha,mu1')
title(ha,'mu1')
ha = subplot(4,2,2,'parent',hf);
plot(ha,PE1')
title(ha,'PE1')
ha = subplot(4,2,3,'parent',hf);
plot(ha,mu2')
title(ha,'mu2')
ha = subplot(4,2,4,'parent',hf);
plot(ha,PE2')
title(ha,'PE2')
ha = subplot(4,2,5,'parent',hf);
plot(ha,mu3')
title(ha,'mu3')
ha = subplot(4,2,6,'parent',hf);
plot(ha,PE3')
title(ha,'PE3')

ha = subplot(4,2,8,'parent',hf);
plot(ha,PE4')
title(ha,'PE4')



% [mu1;mu2;mu3;PE1;PE2;PE3;vec(Si1);vec(Si2);vec(Si3);F];
