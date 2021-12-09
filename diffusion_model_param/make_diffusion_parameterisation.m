%% Load experiment and OCV lookup 

clear;
load brOCV;
load CC_25;


k0=[0.01195; 50; 1000; 1.2]; %R0, I0,  tauD, kd
lowBound = [k0(1)/3;  1;    700;   k0(4)/2];
upBound = [k0(1)*3;   200;  2800;  k0(4)*1.5];

iniPopSpread=2;
initGaPopSize = 50;
randInitPopFactor = max(min(randn(length(k0),initGaPopSize)/2,1),-1);
initGaPop = k0./iniPopSpread.*randInitPopFactor + repmat(k0,1,initGaPopSize);
genAlgOpts.InitialPopulationMatrix = initGaPop';		
genAlgOpts.MaxTime = 60*10;

maxOptTime = 60*10; % Set optimisation time
maxHybridFminIters = 500; % Max number of IPA iterations


% IPA optimisation options
hybridopts = optimoptions('fmincon','Display','none','MaxIterations',maxHybridFminIters);
				
				% GA optimisation options
genAlgOpts = optimoptions('ga','Display','iter','InitialPopulationMatrix',...
				initGaPop','UseParallel', true,'MaxTime',maxOptTime,'MaxGenerations',1000.*length(k0),...
				'HybridFcn',{@fmincon,hybridopts});
 ocvData=BrOcv;

%% 5C

currData=currentSeries{7};
socData=socRefSeries{7};
voltageData=voltageSeries{7};
tempData=tempSeries{7};
timeData=1:length(currData)';

genAlgOpts.MaxTime = 60*3;
maxOptTime = 60*3;
globObj = @(k)Diffusion_Param_Optim_Function(k,currData,timeData,socData,voltageData,tempData,ocvData);
[xOptTmp,fOpt(1,1),~,~] = ga(globObj,length(k0),[],[],[],[],lowBound,upBound,[],genAlgOpts);

LUT25(:,7)=xOptTmp;
save param_5C.mat xOptTmp

%% 4C
currData=currentSeries{6}(1:450);
socData=socRefSeries{6}(1:450);
voltageData=voltageSeries{6}(1:450);
tempData=tempSeries{6}(1:450);
timeData=1:length(currData)';

k0=xOptTmp;
genAlgOpts.MaxTime = 60*3;
maxOptTime = 60*3;
globObj = @(k)Diffusion_Param_Optim_Function(k,currData,timeData,socData,voltageData,tempData,ocvData);
[xOptTmp,fOpt(1,1),~,~] = ga(globObj,length(k0),[],[],[],[],lowBound,upBound,[],genAlgOpts);

LUT25(:,6)=xOptTmp;
save param_4C.mat xOptTmp
%% 3C
currData=currentSeries{5}(1:700);
socData=socRefSeries{5}(1:700);
voltageData=voltageSeries{5}(1:700);
tempData=tempSeries{5}(1:700);
timeData=1:length(currData)';

k0=xOptTmp;
genAlgOpts.MaxTime = 60*5;
maxOptTime = 60*5;
globObj = @(k)Diffusion_Param_Optim_Function(k,currData,timeData,socData,voltageData,tempData,ocvData);
[xOptTmp,fOpt(1,1),~,~] = ga(globObj,length(k0),[],[],[],[],lowBound,upBound,[],genAlgOpts);

LUT25(:,5)=xOptTmp;
save param_3C.mat xOptTmp
%% 2C
currData=currentSeries{4}(1:1200);
socData=socRefSeries{4}(1:1200);
voltageData=voltageSeries{4}(1:1200);
tempData=tempSeries{4}(1:1200);
timeData=1:length(currData)';

k0=xOptTmp;
globObj = @(k)Diffusion_Param_Optim_Function(k,currData,timeData,socData,voltageData,tempData,ocvData);
[xOptTmp,fOpt(1,1),~,~] = ga(globObj,length(k0),[],[],[],[],lowBound,upBound,[],genAlgOpts);

LUT25(:,4)=xOptTmp;
save param_2C.mat xOptTmp
%% 1C

currData=currentSeries{3}(500:2840);
socData=socRefSeries{3}(500:2840);
voltageData=voltageSeries{3}(500:2840);
tempData=tempSeries{3}(500:2840);
timeData=1:length(currData)';

genAlgOpts.MaxTime = 60*6;
maxOptTime = 60*6;
k0=xOptTmp;
globObj = @(k)Diffusion_Param_Optim_Function(k,currData,timeData,socData,voltageData,tempData,ocvData);
[xOptTmp,fOpt(1,1),~,~] = ga(globObj,length(k0),[],[],[],[],lowBound,upBound,[],genAlgOpts);

LUT25(:,3)=xOptTmp;
save param_1C.mat xOptTmp

%% 0.5C
currData=currentSeries{2}(3000:6500);
socData=socRefSeries{2}(3000:6500);
voltageData=voltageSeries{2}(3000:6500);
tempData=tempSeries{2}(3000:6500);
timeData=1:length(currData)';


k0=xOptTmp;
globObj = @(k)Diffusion_Param_Optim_Function(k,currData,timeData,socData,voltageData,tempData,ocvData);
[xOptTmp,fOpt(1,1),~,~] = ga(globObj,length(k0),[],[],[],[],lowBound,upBound,[],genAlgOpts);

LUT25(:,2)=xOptTmp;
save param_05C.mat xOptTmp
%% 0.1 C


voltageData=voltageSeries{1}(22000:end);
currData=(currentSeries{1}(22000:end));
socData=socRefSeries{1}(22000:end);
tempData=tempSeries{1}(2200:end);
timeData=1:length(currData)';


genAlgOpts.MaxTime = 60*15;
maxOptTime = 60*15;
k0=xOptTmp;
globObj = @(k)Diffusion_Param_Optim_Function(k,currData,timeData,socData,voltageData,tempData,ocvData);
[xOptTmp,fOpt(1,1),~,~] = ga(globObj,length(k0),[],[],[],[],lowBound,upBound,[],genAlgOpts);

LUT25(:,1)=xOptTmp;
save param_01C.mat xOptTmp


%% Run and plot
params=xOptTmp;

[Vsim,Vdiff]=diffusion_model_run(params,currData,timeData,socData,tempData,ocvData);

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

save LUT25.mat LUT25;



 %% Test run
% % 
load LUT25.mat
load lincc_25.mat;
params=LUT25;
currData_t=lincc_25{3}(3:600,1);
socData_t=lincc_25{3}(3:600,3);
voltageData_t=lincc_25{3}(3:600,2);
tempData_t=35*ones(1391,1);
timeData=3:length(currData_t)+2;
% 

 Vsim=diffusion_model_run_lut(params,currData_t,timeData,socData_t,ocvData);
error_val=mean(abs(Vsim-voltageData_t));
% 
% 
% 
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



plot(socData_t,voltageData_t-Vsim);
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
