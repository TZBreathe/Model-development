%% Load experiment and OCV lookup 

load brOCV;
load experiment;

currData=currentSeries(6950:7950,2);
timeData=currentSeries(6950:7950,1)-currentSeries(6950,1);
socData=socRefSeries(6950:7950,2);
voltageData=voltageSeries(6950:7950,2);
timeData=1:length(currData);
ocvData=BrOcv;


%% Fitting model

% Initial values and bounds

k0=[0.01195 0.0095 4.38 1200 1.1]'; %R0, R1, tauRC, tauD, kd
lowBound = [[k0(1:3)/3];k0(4)/1.5;k0(5)/3];
upBound = [[k0(1:3)*3];k0(4)*1.5;k0(5)*1.5];

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


%%
params=xOptTmp;
[Vsim,soc]=diffusion_model_run(params,currData,timeData,socData,ocvData);

hold on
plot(timeData,Vsim,'bl');

plot(timeData,voltageData);


plot(socData,socData); plot(socData,soc,'red');
plot(timeData,currData/3)





