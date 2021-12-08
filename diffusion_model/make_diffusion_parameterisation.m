%% Load experiment and OCV lookup 

clear;
load brOCV;
load CC_25;


currData=smooth(currentSeries{1}(22000:end),100);
socData=socRefSeries{1}(22000:end);
for i=1:length(socData)
    if socData(i)>1
        socData(i)=1;
    end
end
voltageData=voltageSeries{1}(22000:end);

% currData=currentSeries{7};
% socData=socRefSeries{7};
% voltageData=voltageSeries{7};

% currData=currentSeries{6}(1:505);
% socData=socRefSeries{6}(1:505);
% voltageData=voltageSeries{6}(1:505);

% currData=currentSeries{5}(1:780);
% socData=socRefSeries{5}(1:780);
% voltageData=voltageSeries{5}(1:780);

% currData=currentSeries{4}(1:1290);
% socData=socRefSeries{4}(1:1290);
% voltageData=voltageSeries{4}(1:1290);

% currData=currentSeries{3}(1:2840);
% socData=socRefSeries{3}(1:2840);
% voltageData=voltageSeries{3}(1:2840);

% currData=currentSeries{2}(1:6500);
% socData=socRefSeries{2}(1:6500);
% voltageData=voltageSeries{2}(1:6500);

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
genAlgOpts.MaxTime = 60*10;

maxOptTime = 60*10; % Set optimisation time
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

save param_01C.mat xOptTmp

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

error_fit=mean(abs(Vsim-voltageData))
 

%% save lut

load param_01C.mat
LUT25(:,1)=xOptTmp;
load param_05C.mat
LUT25(:,2)=xOptTmp;
load param_1C.mat
LUT25(:,3)=xOptTmp;
load param_2C.mat
LUT25(:,4)=xOptTmp;
load param_3C.mat
LUT25(:,5)=xOptTmp;
load param_4C.mat
LUT25(:,6)=xOptTmp;
load param_5C.mat
LUT25(:,7)=xOptTmp;




 %% Test run
% % 
load LUT25.mat
load lincc_25.mat;
params=LUT25;
currData_t=lincc_25{7}(10:1400,1);
socData_t=lincc_25{7}(10:1400,3);
voltageData_t=lincc_25{7}(10:1400,2);
timeData=1:length(currData_t)';
% 

 Vsim=diffusion_model_run_lut(params,currData_t,timeData,socData_t,ocvData);
% error_val=mean(abs(Vsim-voltageData_t));
% 
% 
% 
figure();
hold on;
plot(socData_t,Vsim,'bl');
plot(socData_t,voltageData_t);
xlabel('SoC');
ylabel('Voltage');

% yyaxis right
% plot(socData_t,currData_t,'red');
% ylabel('Current','color','red');
% legend('Model','Exp','Current','location','southeast');
% hold off;
% 
% % r=linspace(0,1,20);
% % hold on;
% % for i=1:200:1000
% %     plot(r,socr(i,:));
% % end
% % xlabel('r');
% % ylabel('SOC');
% % legend('1s', '200s','400s','600s','8000s');
% % hold off;
