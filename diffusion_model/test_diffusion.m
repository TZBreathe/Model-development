%% Load experiment and OCV lookup 

clear;
load brOCV;
load experiment;
dataEnd=11000;
dataBegin=1;

currData=currentSeries(6950:dataEnd,2);
socData=socRefSeries(6950:dataEnd,2);
voltageData=voltageSeries(6950:dataEnd,2);
timeData=1:length(currData);
ocvData=BrOcv;


%% Fitting model

% Initial values and bounds

k0=[0.01195; 100; 0.2; 1000; 0.75; 7000; 1.2]; %R0, I0, alpha, tauD, beta,gama, kd
lowBound = [k0(1)/2;  1;       0.01;  800;     0.4;    1000;  k0(7)/2];
upBound = [k0(1)*2;   300;     0.4;   2000;    0.8;   20000;  k0(7)*2];

iniPopSpread=2;
initGaPopSize = 50;
randInitPopFactor = max(min(randn(length(k0),initGaPopSize)/2,1),-1);
initGaPop = k0./iniPopSpread.*randInitPopFactor + repmat(k0,1,initGaPopSize);
genAlgOpts.InitialPopulationMatrix = initGaPop';		
genAlgOpts.MaxTime = 60*5;

maxOptTime = 60*5; % Set optimisation time
maxHybridFminIters = 500; % Max number of IPA iterations

globObj = @(k)Diffusion_Param_Optim_Function(k,currData,timeData,socData,voltageData,ocvData);
% IPA optimisation options
hybridopts = optimoptions('fmincon','Display','none','MaxIterations',maxHybridFminIters);
				
				% GA optimisation options
genAlgOpts = optimoptions('ga','Display','iter','InitialPopulationMatrix',...
				initGaPop','UseParallel', true,'MaxTime',maxOptTime,'MaxGenerations',1000.*length(k0),...
				'HybridFcn',{@fmincon,hybridopts});
				
% Perform the optimisation with the bounds
[xOptTmp,fOpt(1,1),~,~] = ga(globObj,length(k0),[],[],[],[],lowBound,upBound,[],genAlgOpts);


%% Run and plot
params=xOptTmp;
[Vsim,soc,Vdiff]=diffusion_model_run(params,currData,timeData,socData,ocvData);

figure();
hold on;
plot(socData,Vsim,'bl');
plot(socData,voltageData);
hold off;

figure();
plot(timeData,Vsim-voltageData);

x=linspace(0,1,100);

figure();
hold on
for timestep=1:length(currData)
Vbv(timestep)=0.0256*2*asinh(currData(timestep)/(params(2)*((socData(timestep)+params(3))^1)*(1-socData(timestep)+params(3)))^1);
end
plot(socData,Vbv);
plot(socData,Vdiff,'red');

hold off


figure();
hold on;
plot(socData,socData); plot(socData,soc,'red');
plot(socData,currData/100);
hold off;

%% Test run

load experiment_validation;
load Opt_params
dataEnd=11000;
dataBegin=1;

currData=currentSeries(6950:dataEnd,2);
socData=socRefSeries(6950:dataEnd,2);
voltageData=voltageSeries(6950:dataEnd,2);
timeData=1:length(currData);

params=xOptTmp;
[Vsim,soc,Vdiff]=diffusion_model_run(params,currData,timeData,socData,ocvData);

figure();
hold on;
plot(socData,Vsim,'bl');
plot(socData,voltageData);
hold off;






