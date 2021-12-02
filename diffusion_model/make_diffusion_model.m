%% Load experiment and OCV lookup 

clear;
load brOCV;
load lincc_25;

% lets use Elbrus 15A 0 mV as training
currData=lincc_25{5}(2:1400,1);
socData=lincc_25{5}(2:1400,3);
voltageData=lincc_25{5}(2:1400,2);

timeData=1:length(currData)';
ocvData=BrOcv;


%% Fitting model

% Initial values and bounds

k0=[0.01195; 100; 0.01; 1000; 0.75; 7000; 1.2]; %R0, I0, alpha, tauD, beta,gama, kd
lowBound = [k0(1)/2;  1;       0.01;   800;     0.5;    800;  k0(7)/2];
upBound = [k0(1)*2;   300;     0.2;   2000;    0.8;  15000;  k0(7)*2];

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

save Opt_params.mat xOptTmp

%% Run and plot
params=xOptTmp;

[Vsim,soc]=diffusion_model_run(params,currData,timeData,socData,ocvData);

figure();
hold on;
plot(socData,Vsim,'bl');
plot(socData,voltageData);
xlabel('SoC');
ylabel('Voltage');

yyaxis right
plot(socData,currData,'red');
ylabel('Current','color','red');
legend('Model','Exp','Current','location','southeast');
hold off;


% figure();
% hold on;
% plot(timeData,Vsim-voltageData);
% xlabel('Time');
% ylabel('Voltage error');
% hold off;
% error_fit=rms(Vsim-voltageData);
% x=linspace(0,1,100);

% figure();
% hold on
% for timestep=1:length(currData)
% Vbv(timestep)=0.0256*2*asinh(currData(timestep)/(params(2)*((socData(timestep)+params(3))^1)*(1-socData(timestep)+params(3)))^1);
% end
% plot(socData,Vbv);
% plot(socData,Vdiff,'red');
% 
% hold off


% figure();
% hold on;
% plot(socData,socData); plot(socData,soc,'red');
% plot(socData,currData/100);
% hold off;

%% Test run
% let's test with Elbrus 20A -20mv

currData_t=lincc_25{7}(2:1400,1);
socData_t=lincc_25{7}(2:1400,3);
voltageData_t=lincc_25{7}(2:1400,2);
timeData=1:length(currData_t)';

params=xOptTmp;
[Vsim,socr]=diffusion_model_run(params,currData_t,timeData,socData_t,ocvData);
error_val=mean(abs(Vsim-voltageData_t));



figure();
hold on;
plot(socData_t,Vsim,'bl');
plot(socData_t,voltageData_t);
xlabel('SoC');
ylabel('Voltage');

yyaxis right
plot(socData_t,currData_t,'red');
ylabel('Current','color','red');
legend('Model','Exp','Current','location','southeast');
hold off;

% r=linspace(0,1,20);
% hold on;
% for i=1:200:1000
%     plot(r,socr(i,:));
% end
% xlabel('r');
% ylabel('SOC');
% legend('1s', '200s','400s','600s','8000s');
% hold off;
