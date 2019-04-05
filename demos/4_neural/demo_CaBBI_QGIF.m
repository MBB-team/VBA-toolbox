% This script inverts the QGIF model for the given fluorescence trace using
% CaBBI method: [ Rahmati. et al. 2016, PLoS Comput Biol 12(2) ]
%
% To user:
% I) in Step 1, determine the directory to the fluorescence trace (a real row vector);
%    as default the code will download and use our published data for CaBBI 
% II) in Step 2, determine the sampling rate (Hz) used for recording the fluorescence trace
% III) run the code
% IV) the main outputs of the CaBBI are listed in Step 11. 
%
% Further options:
% 1) In Step 10, if necessary, you can change the spike-detection threshold which is used to extract 
%    spike-times from inferred voltage traces. Such change may increase the spike-detection accuracy. 
% 2) In Step 6, if necessary, you can change some internal parameters of the optimization process.
% 3) For further description of the toolbox options, please see:
%    http://mbb-team.github.io/VBA-toolbox/wiki/Controlling-the-inversion-using-VBA-options/
%    http://mbb-team.github.io/VBA-toolbox/wiki/VBA-output-structure/


clear all; close all

%% Step 1. Importing the fluorescence trace
% import the fluorescence trace
dataFolder = [pwd filesep 'Data' filesep];
addpath('dataFolder'); 
Fluor_trace_name = 'fluorescence_data8';                                   % file name of the fluorescence trace
Fluor_trace_dir = [dataFolder,Fluor_trace_name,'.mat'];
 
if exist(Fluor_trace_dir,'file') == 0                                      % check whether the fluorescence data have been already downloaded
    url = 'https://github.com/VahidRahmati/Sample_Data/blob/master/Sample_Data.zip?raw=true';
    urlwrite(url,[pwd filesep 'Data.zip']);                                       % download the fluorescence data
    unzip('Data.zip',[pwd filesep 'Data']); 
end

ydata = importdata(Fluor_trace_dir);                                       % the fluorescence trace before pre-processing (a real row vector)


%% Step 2. Determining the sampling rate used to record the fluorescence trace
sampling_rate = 22.6;                                                        % in [Hz], for our data enter 22.6
display(['the sampling rate used to record the fluorescence trace is << ',num2str(sampling_rate),' Hz >>']);


%% Step 3. Pre-processing
% prparing the data
warning('on')
if ~all(isfinite(ydata))                                                   % check whether there is any NaN and/or Inf in the given fluorescence trace
    ydata(~isfinite(ydata)) = [];                                          % reovme the NaN and/or Inf elements from the trace
    warning('the given fluorescence traces had some NaN and/or Inf elements, which were removed')
end

if ~isreal(ydata)                                                          % check whether the given fluorescence trace contains any complex number
    error('..... the given fluorescence traces contains complex numbers .....')
end

% normalizing
if size(ydata,1) ~= 1; ydata = ydata'; end                                 % fluorescence trace should be a row vector
ydata = (ydata - mean(ydata))./(max(mean(ydata),1));                       % relativer fluorescncee trace
ydata = ydata - median(ydata);                                             % removing the median

% removing the slow drifts 
xdata = linspace(0,1,length(ydata));                                                   % vector of frames
degree = '4';      
if exist('fittype')
    % Polynomial degree (4 or 3)
    ft = fittype( ['poly',degree] );                                           % define the model to fit the slow drifts 
    opts = fitoptions( ft );
    opts.Lower = [-Inf -Inf -Inf -Inf -Inf];
    opts.Robust = 'Bisquare';
    opts.Upper = [Inf Inf Inf Inf Inf];
    warning('off')
    [fitresult, gof] = fit( xdata', ydata', ft, opts );                          % fitting the polynomial model to fluorecence trace
    warning('on')
    Coeffs = coeffvalues(fitresult);
    if degree == '4'
        fittedcurve = Coeffs(1)*xdata.^4 + Coeffs(2)*xdata.^3 + Coeffs(3)*xdata.^2 + Coeffs(4)*xdata + Coeffs(5);
    elseif degree == '3'
        fittedcurve = Coeffs(1)*xdata.^3 + Coeffs(2)*xdata.^2 + Coeffs(3)*xdata.^1 + Coeffs(4);
    end
else
    % less robust fit as fallback
    fittedcurve = polyval(polyfit(xdata,ydata,4),xdata) ;
end

Fluor_trace = ydata - fittedcurve;                                         % pre-processed fluorescence trace 


%% Step 4. Scaling the data
scale       = 1;                                                           % preferably set it a number between e.g. 0.4 and 1
Fluor_trace = scale* Fluor_trace/max(Fluor_trace);
xdata = 1:length(ydata);                                                   % vector of frames
plot(xdata, scale*ydata/max(ydata) ,'b',xdata, Fluor_trace(1:numel(Fluor_trace)),'r');
legend('original','pre-processed')
xlabel('time [frame]')
ylabel('a.u.')
title('fluorescence trace')


%% Step 5. The generative model 
% assigning the generative model
f_fname = @f_CaBBI_QGIF;                                                   % evolution function
g_fname = @g_CaBBI;                                                        % observation function

% setting the parameters of the QGIF model used for inversion
nFrames = numel(Fluor_trace);                                              % number of recording frames
dt = 0.2;                                                                  % [ms], time-step size of the QGIF model
inF.dt = dt; 
inF.k_Ca0 = 0.002;                                                         % scale parameter of [Ca2+] kinetics
TauCa_realTime = 7500;                                                     % [ms], the prior for decay time-constant of calcium transients in real time
dt_real = 1000/sampling_rate;                                              % [ms], frame duration (temporal precision) of fluorescence data
inF.Tau_Ca0 = TauCa_realTime*dt/dt_real;                                   % efective prior mean for tau_Ca in inversion; see Eqn. 17
inG.k_F0 = 5;                                                              % scale parameter of the fluorescence observations
inG.ind = 2;                                                               % index of [Ca2+] variable


%% Step 6. Setting the internal parameters of the optimization (VB-Laplace) method
options.backwardLag = 2;                                                   % lag parameter (selcet either 2 or 3); increasing it e.g. to 3 (or higher) "may" provide more accurate inversions 
options.updateHP    = 1;                                                   % this allows the precison (hyper) parameters to be learnt from data
options.MaxIterInit = 0;  
options.MaxIter     = 2;  % or e.g. 3                                      % maximum number of optimization iterations; increasing it will lead to a  
                                                                           % complete optimization process, although it may take long time to be converged)
                                                                           
options.DisplayWin  = 0;                                                   % set to zero to speed up the inversion
options.verbose     = 1;                                                   % set to zero to speed up the inversion
options.decim       = 1;                                                   % this determine the micro-time resolution; increasing it e.g. to 2 (or higher) "may" provide more accurate inversions


%% Step 7. Assigning the prior densities: (mu --> mean), (Sigma --> variance)
% prior densities for initial conditions 
priors.muX0         = [ -63.62; 5.65];                                     % [V(0), [Ca2+](0)]
priors.SigmaX0(1,1) = 100;
priors.SigmaX0(2,2) = 25;

% prior densities for evolution parameters
priors.muTheta         = [0; 0; 0]; 
priors.SigmaTheta(1,1) = 5;                                                % for X_Tau_Ca
priors.SigmaTheta(2,2) = 100;                                              % for Ca_base 
priors.SigmaTheta(3,3) = 5;                                                % for X_k_Ca

% prior densities for observation parameters
priors.muPhi         = [min(Fluor_trace)/2; 0]; 
priors.SigmaPhi(1,1) = 0.25;                                               % for d_F
priors.SigmaPhi(2,2) = 1;                                                  % for X_k_F 

% prior densities for precision parameters: Gamma distributions
priors.a_alpha = 1e0;
priors.b_alpha = 1e0;
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;

% building the precision matrix for evolution variables
 priors.iQx = cell(nFrames,1);
for t = 1:nFrames 
    dq(1)           = 5e-3;                                                % on V 
    dq(2)           = 10;                                                  % on [Ca2+]
    priors.iQx{t}   = diag(dq);
end

% building the precision vector for fluorescence observations    
priors.iQy = cell(nFrames,1);
for t=1:nFrames
    priors.iQy{t} = (3);                                                   % on F
end

% setting the option structure
options.inF     = inF;
options.inG     = inG;
options.priors  = priors;

% assigning the dimensions for parameter space
dim.n       = 2;                                                           % number of neural variables
dim.n_phi   = 2;                                                           % number of observation parameters
dim.n_theta = 3;                                                           % number of evolution parameters


%% Step 8. Inversion: inverting the QGIF model for the given fluorescence trace
[posterior, out] = VBA_NLStateSpaceModel(Fluor_trace,zeros(1,numel(Fluor_trace)),f_fname,g_fname,dim,options);


%% Post-processing
%% Step 9. Extracting the estimated Voltage and [Ca2+] traces (i.e. the posterior means)

V_inferred = posterior.muX(1,1:end);                                       % inferred Voltage (posterior mean)
Ca_inferred = posterior.muX(end,1:end);                                    % inferred [Ca2+] kinetics (posterior mean)

% Plotting the estimated voltage and [Ca2+] traces
Calcium_basal = 50;                                                        % in [nM], the basal [Ca2+] concentration  

% Plot the pre-processed fluorescence trace
figure
subplot(411); plot(1:nFrames, Fluor_trace,'g','LineWidth',2)
legend('fluorescence')
set(gca,'xLim',[0 nFrames])
axis tight

% Plot the posterior mean of [Ca2+] trace
subplot(412); plot(1:nFrames, Ca_inferred+Calcium_basal,'r','LineWidth',2)
legend('inferred [Ca2+]')
set(gca,'xLim',[0 nFrames])
axis tight

% plot the posterior mean of voltage trace
subplot(413); plot(1:nFrames, V_inferred,'b','LineWidth',2)
set(gca,'xLim',[0 nFrames])
axis tight


%% Step 10. Detecting spike times (or event onsets) from inferred voltage trace
Detection_threshold  = 0;                                                  % the spike (or event) detection threshold
                                                                           % (if necessary, change the detection threshold to e.g. -10)                                                                                                                                                     
idxSpikes = zeros(1,numel(V_inferred));
count  = 1;
for j=1:numel(V_inferred)-1
    
    if( V_inferred(j+1)>Detection_threshold && V_inferred(j)<Detection_threshold )
        
        if (count > 1) && ( (j-indices(count-1)) > ceil(6/dt) )             
            idxSpikes(j) = j;
            indices(count) = j;
            count = count+1;
        end
        
        if count == 1
            idxSpikes(j) = j;
            indices(count) = j;
            count = count+1;
        end
        
    end
    
end
idxSpikeToShow = idxSpikes;
idxSpikeToShow(idxSpikeToShow~=0) = 1;

subplot(414); plot(1:nFrames,idxSpikeToShow,'b','LineWidth',2)
legend('inferred spiketimes');
xlabel('time [frame]')
set(gca,'xLim',[0 nFrames])
axis tight

subplot(413); hold on; plot(1:nFrames, Detection_threshold*ones(1,nFrames),'-.k');
legend('inferred voltage','threshold'); hold off

% extracted spike (or event) times
SpikeTimes_Frame = idxSpikeToShow(idxSpikeToShow~=0)';                     % reconstructed spike (or event) times, in [frame]
SpikeTimes_millisecond = SpikeTimes_Frame*1000/sampling_rate;              % reconstructed spike (or event) times, in [ms]


%% Step 11. Main outputs of CaBBI method
CaBBI.nFrames           = nFrames;                                         % the number of recordings frames of the given fluorescence trace                          
CaBBI.sampling_rate     = sampling_rate;                                   % the sampling rate used for recording the fluorescence trace, in [Hz]
CaBBI.Fluor_trace       = Fluor_trace;                                     % the pre-processed fluorescebce trace    
CaBBI.V_inferred        = V_inferred;                                      % inferred Voltage trace (posterior mean) 
CaBBI.Ca_inferred       = Ca_inferred;                                     % inferred [Ca2+] trace (posterior mean)
CaBBI.SpikeTimes_Frame  = SpikeTimes_Frame;                                % reconstructed spike (or event) times, in [frame]
CaBBI.SpikeTimes_millisecond  = SpikeTimes_millisecond;                    % reconstructed spike (or event) times, in [ms]

%% cleanup

delete([pwd filesep 'Data.zip']);  % delete the fluorescence data
rmdir([pwd filesep 'Data'],'s');    