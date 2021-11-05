 % Confidential
 % 
 % -SYNOPSIS
 % Function to perform ECN parameterisation to all cycles of DCIR test data
 %
 % -NOTES
 %    Version:      1.0
 %    Author:       Breathe Battery Technologies
 %    E-mail:       christian.korte@breathe.technology
 %    File: 		ECN_Paramsation.m
 %    Copyright (c) 2021 Breathe Battery Technologies Ltd.
 % -LINK
 %    https://www.breathe.technology/
 % 
 % -LICENSE
 % Copyright Breathe Battery Technologies Ltd. May 2021 - All Rights Reserved
 % Provided under licence to Rimac Automobili d.o.o.

function[Results, socBp]  = ECN_Paramsation(cycleTabs,EcnSettings,StpInfo,currentFile,OcvLuts)
overallTimer = tic;

fastTestRun = EcnSettings.fastTestRun;
debugPlot = EcnSettings.debugPlot;

fastPulseStp = StpInfo.fastPulseStp;
slowPulseStp = StpInfo.slowPulseStp;
dvdtStp = StpInfo.dvdtStep;

optiParamnNames = {'R_0','R_1','tau_1','R_2','tau_2'};

%% 1 Fit process

if EcnSettings.cycles == 0 % Fit all cycles
	cycleIter = 1:length(cycleTabs);
else % Fit only desired cycles
	cycleIter = EcnSettings.cycles;
end

% Loop through all cycles
for i_c = 1:length(cycleIter)
	
	fprintf(['\n\n\t' datestr(now,'dd/mm hh:MM:ss') '\t Fitting cycle ' num2str(cycleIter(i_c)) ', Number ' num2str(i_c) ' of ' num2str(length(cycleIter)) ': \t'])
	
	% Get cycle index of current cycle
	cycIdx = cycleIter(i_c);
	
	timeData = cycleTabs{cycIdx}.timeTest - min(cycleTabs{cycIdx}.timeTest);
	currData = cycleTabs{cycIdx}.currCell;
	stpData = cycleTabs{cycIdx}.stepIdx;
	socData = cycleTabs{cycIdx}.socOut;
	ocvData = cycleTabs{cycIdx}.ocpotCell;
	voltData = cycleTabs{cycIdx}.voltCell;
	
	assert(all(voltData > 0), 'Voltage data not greater than 0')
	assert(all(ocvData > 0), 'OCV data not greater than 0')

	%% 1.1 Fit OCV data to resting voltage before pulses (the estimated OCV (ocvEst))
	ocvAdjusted = zeros(length(stpData),1);
	
	% Get indices before and after fast and slow pulses
	idxBeforeFastPulse = find(stpData == fastPulseStp,1,'first');
	idxEndFastPulse = find(stpData == fastPulseStp,1,'last');
	idxBeforeSlowPulse = find(stpData == slowPulseStp,1,'first');
	idxEndSlowPulse = find(stpData == slowPulseStp,1,'last');
	
	assert(	~isempty(idxBeforeFastPulse) && ...
			~isempty(idxEndFastPulse) && ...
			~isempty(idxBeforeSlowPulse) && ...
			~isempty(idxEndSlowPulse),'At least one pulse index for the OCV adjustment is empty')
	
	% Get estimated OCV before and after fast and slow pulses
	ocvEstBeforeFastPulse = voltData(idxBeforeFastPulse -1);
	ocvEstBeforeSlowPulse = voltData(idxBeforeSlowPulse - 1);
	ocvEstAfterSlowPulse = voltData(end);
	
	% Set OCV before the fast pulse to estimated OCV before the fast pulse
	ocvAdjusted(1:idxBeforeFastPulse -1) = ocvEstBeforeFastPulse;
	
	% Adjust Hf Pulse OCV -------------------------------------------------
	ocvLutTmpFastPulse = ocvData(idxBeforeFastPulse:idxEndFastPulse); % Get OCV for fast pulse according to OCV LUT
	ocvLutTmpFastPulseZero = ocvLutTmpFastPulse - ocvLutTmpFastPulse(1); % Set it to zero at the beginning to compare to estimated OCV
	ocvEstFastPulseRange =  ocvEstBeforeSlowPulse - ocvEstBeforeFastPulse; % Calculate estimated OCV range 
	
	% Transform the LUT OCV to the same range as the estimated OCV
	ocvLutFastPulseTmpAdjustedZero = ocvLutTmpFastPulseZero.* ocvEstFastPulseRange/(ocvLutTmpFastPulseZero(end) - ocvLutTmpFastPulseZero(1));
	% Add the OCV shape to estimated OCV before the fast pulse
	ocvAdjusted(idxBeforeFastPulse:idxEndFastPulse) = ocvLutFastPulseTmpAdjustedZero + ocvEstBeforeFastPulse;
	
	% Set OCV between fast and slow pulse to OCV after fast pulse
	ocvAdjusted(idxEndFastPulse : idxBeforeSlowPulse -1) = ocvAdjusted(idxEndFastPulse);
	
	% Adjust Lf Pulse OCV -------------------------------------------------
	ocvLutTmpSlowPulse = ocvData(idxBeforeSlowPulse:idxEndSlowPulse); % Get OCV for slow pulse according to OCV LUT
	ocvLutTmpSlowPulseZero = ocvLutTmpSlowPulse - ocvLutTmpSlowPulse(1); % Set it to zero at the beginning to compare to estimated OCV
	ocvEstSlowPulseRange =  ocvEstAfterSlowPulse - ocvEstBeforeSlowPulse; % Calculate estimated OCV range 
	
	% Transform the LUT OCV to the same range as the estimated OCV
	ocvLutSlowPulseTmpAdjustedZero = ocvLutTmpSlowPulseZero.* ocvEstSlowPulseRange/(ocvLutTmpSlowPulseZero(end) - ocvLutTmpSlowPulseZero(1));
	
	% Add the OCV shape to estimated OCV before the slow pulse
	ocvAdjusted(idxBeforeSlowPulse:idxEndSlowPulse) = ocvLutSlowPulseTmpAdjustedZero + ocvEstBeforeSlowPulse;
	
	% Set OCV after slow pulse to OCV after fast pulse	
	ocvAdjusted(idxEndSlowPulse :length(voltData) ) = ocvAdjusted(idxEndSlowPulse);
	
	% Overpotential data = measured voltage - OCV (adjusted)
	voltOpData = voltData - ocvAdjusted;
		
	%% 1.2 Begin Fitting
	
	socBp(cycIdx) = socData(1); % The soc at the beginning of the data is the SoC for this cycle/breakpoint
	
	if EcnSettings.skipLastIteration
		fitIterations = [1:3]; % Fitting iteration numbers	
	else
		fitIterations = [1:4]; % Fitting iteration numbers
	end
	% values i_f can take:
	% 1: estimate k0
	% 2: Fit HFP to fast pulse
	% 3: Fit LFP/AP to slow pulse
	% 4: Fit AP to all pulses (optional)

	i_f = 1;
	repFlag = false; % Flag that indicates an iteration is a repeat (true = repeat)
	
	while i_f <= max(fitIterations)
		
		%% 1.2.1 Select relevant data for each fitting iteration
		if ~repFlag % Data selection is not needed if this is a repetition
			
			fprintf(['\n\t\tFit iteration ' num2str(i_f) '/' num2str(max(fitIterations)) '\t']);
			switch i_f
				case {1,2} % Fast pulse
					
					if EcnSettings.fitFastRelax % Use data of relaxation
						dataIdx = stpData>=fastPulseStp & stpData<slowPulseStp;
					else % Don't use data of relaxation
						dataIdx = stpData==fastPulseStp;
					end
					
					dataIdx((find(dataIdx,1)-1):(find(dataIdx,1,'last'))) = 1; % extend one data point before where current is 0
					
					% Extract the data used for the the optimisation
					voltFitData{1}= voltOpData(dataIdx);
					currFitData{1}= currData(dataIdx);
					timeFitData{1}= timeData(dataIdx);
					stpFitData{1} = stpData(dataIdx);
					socFitData{1} = zeros(length(stpFitData{1}),1); % SoC is not fitted in this iteration
					
				case 3 % Slow pulse
					
					if EcnSettings.fitSlowRelax
						dataIdx = (stpData>=slowPulseStp) & (1:length(stpData) > 300)'; % Second part is to ignore the first 300 entries in data, which may include previous breakpoint data
						oneHourAfterPulseIdx = find(stpData == dvdtStp & (1:length(stpData) > 300), 1);
						if ~isempty(oneHourAfterPulseIdx)
							dataIdx(oneHourAfterPulseIdx:end) = 0; % Ignore data after first measurement of voltage change
						end
					else
						dataIdx = (stpData==slowPulseStp) & (1:length(stpData) > 300)'; % Second part is to ignore the first 300 entries in data, which may include previous breakpoint data
					end
					
					dataIdx((find(dataIdx,1)-1):(find(dataIdx,1,'last'))) = 1; % extend one data point before where current is 0
					
					% Extract the data used for the the optimisation
					voltFitData{1}= voltOpData(dataIdx);
					currFitData{1}= currData(dataIdx);
					timeFitData{1}= timeData(dataIdx);
% 					socFitData{1} = zeros(length(timeFitData{1}),1); % SoC is not fitted in this iteration
%                     socFitData2Tmp = socData(dataIdx);  % This data is only used in fit iteration 4 below
% 					socFitData1Tmp = socData(dataIdx);					
% 					socRange = (socFitData1Tmp(1) - socFitData2Tmp(end));					
% 					socFitData{1} = -(socFitData1Tmp - socFitData1Tmp(1))/socRange;
                    socFitData{1} = socData(dataIdx);	
					
					
				case 4 % Fast and slow pulse
					
					if EcnSettings.fitFastRelax % Use data of relaxation
						dataIdx = (stpData>=fastPulseStp & stpData<slowPulseStp);
					else % Don't use data of relaxation
						dataIdx = (stpData==fastPulseStp);
					end
					dataIdx((find(dataIdx,1)-1):(find(dataIdx,1,'last'))) = 1; % extend one data point before where current is 0
					
					% Extract the data used for the the optimisation
					
					% Use data from the previous slow pulse...
					voltFitData{2} = voltFitData{1};
					currFitData{2} = currFitData{1};
					timeFitData{2} = timeFitData{1};
					
					% ...and data from fast pulse
					voltFitData{1}= voltOpData(dataIdx);
					currFitData{1}= currData(dataIdx);
					timeFitData{1}= timeData(dataIdx);
					
					socFitData1Tmp = socData(dataIdx);
					
					socRange = (socFitData1Tmp(1) - socFitData2Tmp(end));
					
					socFitData{1} = -(socFitData1Tmp - socFitData1Tmp(1))/socRange;
					socFitData{2} = -(socFitData2Tmp - socFitData1Tmp(1))/socRange;

					% SoC variable should go from 0 to 1 for
					% interpolation
					assert(max(socFitData{2})==1,'SoC does not go to 1');
					assert(min(socFitData{1})==0,'SoC does not go to 0');

			end
		end
		
		assert((length(timeFitData{1}) == length(currFitData{1})) & (length(timeFitData{1}) == length(voltFitData{1})),'All input data must be the same length.');
		
		%% 1.2.2 Set initial guess and determine bounds
		
		%factor by which the fast parameters can change
		constrainFacHfEstimates = EcnSettings.constrainFacHfEstimates(i_f); 
		%factor by which the slow parameters can change
		constrainFacLfEstimates = EcnSettings.constrainFacLfEstimates(i_f); 
					
		% The variables are the global lower and upper bounds
		globalLb = EcnSettings.paramLowerBounds;
		globalUb = EcnSettings.paramUpperBounds;
		
		clear k0 %k0 is this initial guess
		switch i_f
			case 1 % Curve fit fast data, only for the initial guess of first full optimisation
				
				firstTransientIdx = find(stpFitData{1} == fastPulseStp,1); % get index of first voltage spike
				% Estimate R_0 from first voltage pulse
				initRGuess = voltFitData{1}(firstTransientIdx)/currFitData{1}(firstTransientIdx); % R = V/I
				
				if initRGuess < 0.0001 % Sometimes there is a datapoint before voltage increase, leading to a low or negative value of R_0
					initRGuess = voltFitData{1}(firstTransientIdx+1)/currFitData{1}(firstTransientIdx+1); % R = V/I
				end
				
				kCurveFit(1) = initRGuess;
				
				rc1FitIdx = (stpFitData{1} == fastPulseStp); % Get first pulse index
				xFit = timeFitData{1}(rc1FitIdx) - min(timeFitData{1}(rc1FitIdx)); %Time data of pulse
				currFit = mean(currFitData{1}(rc1FitIdx));	% Average value of the current pulse
				% Voltage data of pulse, without the contribution of R_0
				yFit = (voltFitData{1}(rc1FitIdx) - currFitData{1}(rc1FitIdx)*kCurveFit(1,1))/currFit; 
				
				% Settings for curve fit
				fo_exp = fitoptions('Method','NonlinearLeastSquares',...
					'Lower',[0.001,1],...
					'Upper',[1,15],...
					'StartPoint',[0.001 2]);
				
				% Define curve fit equation
				g_exp = fittype('r*(1-exp(-x/(tau)))','options',fo_exp);
				
				% Do curve fit
				[expFitModRc1, ~] = fit(xFit, yFit, g_exp);
				
				% Assign results to kCurveFit variable
				kCurveFit(2) = expFitModRc1.r;
				kCurveFit(3) = expFitModRc1.tau;
				
				xOpt{i_f,i_c} = kCurveFit'; % This is used as the initial guess for k in the optimisation where i_f == 2
				
			case 2 % GA fit fast parameters
				
				lastResults = xOpt{i_f-1,i_c};
				
				% Set bounds, proportional to results of first curve fit
				lowBound = lastResults/EcnSettings.constrainFacHfEstimates(i_f);
				upBound = lastResults*EcnSettings.constrainFacHfEstimates(i_f);
				
				k0 = lastResults; % initial guess
				
				fitRCs = 1; %1 = Fit only R1 C1, 0 = fit all parameters
				
				useSocDependence = false; % Assume SoC is constant
				
			case 3 % GA fit slow parameters to slow pulse
				
				lastResults = xOpt{i_f-1,i_c};
				% Because there is no guess for k(4:5), guess that k0(4) =
				% k(2) and k0(5) = k(3). (k(1:5) = R_0, R_1, C_1, tauD, kD)
				lastResults = [lastResults; 3000; 1.2];
				
				k0 = lastResults;
				% Use bounds based on constrainfactors for first 3 entries,
				% and global bounds for the rest
				lowBound = globalLb;
				lowBound([1 2 3]) = lastResults(1:3)'./EcnSettings.constrainFacHfEstimates(i_f);
				upBound = globalUb;
				upBound([1 2 3]) = lastResults(1:3)'.*EcnSettings.constrainFacHfEstimates(i_f);
				
				fitRCs = 0; %1 = Fit only R1 C1, 0 = fit all parameters
				
				useSocDependence = false;% Assume SoC is constant
				
			case 4
				
				lastResults = xOpt{i_f-1,i_c};
				k0 = lastResults;
				
				lowBound = ([k0(1:3)./constrainFacHfEstimates; k0(4:5)./constrainFacLfEstimates ]);
				upBound = ([k0(1:3).*constrainFacHfEstimates; k0(4:5).*constrainFacLfEstimates ]);
				
				useSocDependence = EcnSettings.fitSocDependence;
				
				% When interpolating SoC, the variable k must have 10
				% entries instead of 5. The entries 6:10 are the parameter
				% values for the second SoC breakpoint (SoC at the end of the slow pulse)
				if useSocDependence
					k0 = repmat(k0,2,1);
					lowBound = repmat(lowBound,2,1);
					upBound = repmat(upBound,2,1);
				end
				
				fitRCs = 0; %1 = Fit only R1 C1, 0 = fit all parameters
		end
			
		%% 1.2.3 Set optimisation options and do optimisation

		if ~repFlag
			iniPopSpread = EcnSettings.initPopSpread(i_f); % initial population spread factor
			optMaxTime = 60*EcnSettings.optiTimeMins(i_f); % Set maximum time for optimisation
		end
		
		if fastTestRun % Reduce population size during a fast test run
			initGaPopSize = 2;
		else
			initGaPopSize = 50;
		end
		
		if i_f >1
			% Generate a normally-distributed random initial population factor between -1 and 1
			randInitPopFactor = max(min(randn(length(k0),initGaPopSize)/2,1),-1);
			
			% Check initial value is the same size as the bounds
			assert(all(size(lowBound) == size(k0)),'Lower bounds not the same size as k0')
			assert(all(size(upBound) == size(k0)),'Upper bounds not the same size as k0')
		end
		
		clear xOptTmp
		switch i_f
			case 1
			case 2
				
				% Calculate initial solution population based on initial
				% guess k0 and the population spread.
				initGaPop = k0./iniPopSpread.*randInitPopFactor + repmat(k0,1,initGaPopSize);
				
				% Set times very low for fast test runs
				if fastTestRun
					maxOptTime = 0.1;
					maxHybridFminIters = 1;
					genAlgOpts.HybridFcn = {};
				else
					maxOptTime = optMaxTime; % Set optimisation time
					maxHybridFminIters = 500; % Max number of IPA iterations
				end
				
				% IPA optimisation options
				hybridopts = optimoptions('fmincon','Display','none','MaxIterations',maxHybridFminIters);
				
				% GA optimisation options
				genAlgOpts = optimoptions('ga','Display','none','InitialPopulationMatrix',...
					initGaPop','UseParallel', true,'MaxTime',maxOptTime,'MaxGenerations',1000.*length(k0),...
					'HybridFcn',{@fmincon,hybridopts});
				
				if debugPlot
					genAlgOpts.PlotFcn = {@gaplotbestf,@gaplotstopping}; % Plots progress of GA
				end
				
				% Define the optimisation function, as a function of
				% variable (k)
				globObj = @(k) ECN_Param_Optim_Function(k,currFitData{1},timeFitData{1},socFitData{1},voltFitData{1},useSocDependence,fitRCs,OcvLuts);
				
				% Perform the optimisation with the bounds
				[xOptTmp,fOpt(i_f,i_c),~,~] = ga(globObj,length(k0),[],[],[],[],lowBound,upBound,[],genAlgOpts);
				
				% Assign solution variable
				xOpt{i_f,i_c} = xOptTmp';
				
			case 3
				% Calculate initial solution population based on initial
				% guess k0 and the population spread.
				initGaPop = k0./iniPopSpread.*randInitPopFactor + repmat(k0,1,initGaPopSize);
				genAlgOpts.InitialPopulationMatrix = initGaPop';
				
				if fastTestRun
				else
					genAlgOpts.MaxTime = optMaxTime;
				end
				
				% Define the optimisation function, as a function of
				% variable (k)
				globObj = @(k) ECN_Param_Optim_Function(k,currFitData{1},timeFitData{1},socFitData{1},voltFitData{1},useSocDependence,fitRCs,OcvLuts);
				
				% Perform the optimisation with the bounds
				[xOptTmp,fOpt(i_f,i_c),~,~] = ga(globObj,length(k0),[],[],[],[],lowBound,upBound,[],genAlgOpts);
				
				% Assign solution variable			
				xOpt{i_f,i_c} = xOptTmp';
				
			case 4
				
				if fastTestRun
				else
					genAlgOpts.MaxTime = optMaxTime;
				end
				
				% Calculate initial solution population based on initial
				% guess k0 and the population spread.
				initGaPop = k0./iniPopSpread.*randInitPopFactor + repmat(k0,1,initGaPopSize);
				genAlgOpts.InitialPopulationMatrix = initGaPop';
				
				% In this case, the optimisation function is the sum of the
				% optimisation functions of the fast pulse and the slow pulse
				globObj = @(k) ECN_Param_Optim_Function(k,currFitData{1},timeFitData{1},socFitData{1},voltFitData{1},useSocDependence,fitRCs,OcvLuts) ...
					+ ECN_Param_Optim_Function(k,currFitData{2},timeFitData{2},socFitData{2},voltFitData{2},useSocDependence,fitRCs,OcvLuts);
				
				% Perform optimisation
				[xOptTmp,fOpt(i_f,i_c),~,~] = ga(globObj,length(k0),[],[],[],[],lowBound,upBound,[],genAlgOpts);
				
				if useSocDependence % Make entry of xOpt the same, whether SoC dependence is fitted or not
					xOpt{i_f,i_c} = [xOptTmp(1:5)' xOptTmp(6:10)'];
				else
					xOpt{i_f,i_c} = [xOptTmp(1:5)' xOptTmp(1:5)'];
				end	
		end
		
		if i_f > 1 % This in done for all optimisation steps, but not for the first curve fit
				
			% Warn if parameters are close to bounds
			if i_f == 4 && EcnSettings.fitSocDependence
				% Get indices of parameters that are close to bounds
				% If a parameter is less than 5% bigger than the lower
				% bound, or more than 5% smaller than the upper bound, a
				% warning is displayed. 
				% (This does not mean that the optimisation result is bad! It is only for information on setting the bounds)
				paramsCloseToBoundsLow = 	((xOptTmp(1:5)./lowBound(1:5)') < 1.05) + ((xOptTmp(6:10)./lowBound(6:10)') < 1.05);
				paramsCloseToBoundsUp = 	((xOptTmp(1:5)./upBound(1:5)') < 1.05) + ((xOptTmp(6:10)./upBound(6:10)') < 1.05);
			else
				% Get indices of parameters that are close to bounds
				paramsCloseToBoundsLow = 	(xOptTmp./lowBound') < 1.05;
				paramsCloseToBoundsUp = 	(xOptTmp./upBound') < 1.05;
			end
					
			% Repeat fit iteration if f is too low, otherwise move to next iteration
			if (fOpt(i_f,i_c) < EcnSettings.fThreshold(i_f)) || repFlag % the repFlag condition means there is at maximum one repeat
				repFlag = false;

				% Print warning if some parameters are close to lower bounds
				if any(xOptTmp'./lowBound < 1.05 )
					fprintf(['\n\t\tWARNING: Parameters ' sprintf('%s, ',optiParamnNames{logical(paramsCloseToBoundsLow)}) 'are close to lower bounds for fit iteration ' num2str(i_f)])
				end
				% Print warning if some parameters are close to upper bounds
				if any(xOptTmp'./upBound > 0.95 )
					fprintf(['\n\t\tWARNING: Parameters ' sprintf('%s, ',optiParamnNames{logical(paramsCloseToBoundsUp)}) 'are close to upper bounds for fit iteration ' num2str(i_f)])
				end

			else
				fprintf('!Repeat Iteration!');
				repFlag = true; % Repeat iteration
				iniPopSpread = iniPopSpread*2; % Increase spread of initial points
				optMaxTime = optMaxTime*1.5; % Increase optimisation time
			end
        end,

		%% 1.2.4 Plot figure showing fit of voltage model to estimated parameters
		if ~repFlag
			
			% For the final iteration, print the f value
			if i_f == max(fitIterations)
				ofjFcn = fOpt(i_f,i_c);
				fprintf(['\n\tf = ' num2str(ofjFcn,3) '.'])
			end
			
			if debugPlot
				
				% Calculate the voltage response with the model using the
				% optimal parameters of this iteration
				switch i_f
					case 2
						close all
						[~, voltModel] = ECN_Param_Optim_Function(xOptTmp,currFitData{1},timeFitData{1},socFitData{1},voltFitData{1},useSocDependence,fitRCs,OcvLuts);
					case 3
						[~, voltModel] = ECN_Param_Optim_Function(xOptTmp,currFitData{1},timeFitData{1},socFitData{1},voltFitData{1},useSocDependence,fitRCs,OcvLuts);
					case 4 %Adjust the data to plot both LF and HF
						[~, voltModelTmp1] = ECN_Param_Optim_Function(xOptTmp,currFitData{1},timeFitData{1},socFitData{1},voltFitData{1},useSocDependence,fitRCs,OcvLuts);
						[~, voltModelTmp2] = ECN_Param_Optim_Function(xOptTmp,currFitData{2},timeFitData{2},socFitData{2},voltFitData{2},useSocDependence,fitRCs,OcvLuts);
						voltModel = [voltModelTmp1; voltModelTmp2];
						timeFitData{1} = [timeFitData{1}; timeFitData{2} - (min(timeFitData{2}) - max(timeFitData{1}))];
						voltFitData{1} = [voltFitData{1}; voltFitData{2}];
				end
				
				if i_f >= 2
					
					% Plot voltage fit and esimated parameters
					figure			
					figTitle = (['Cycle ' num2str(cycleIter(i_c)) ' Fit Iteration ' num2str(i_f)]);
					
					sgtitle(figTitle) % Subfigure Title
					
					% Plot measured and model voltage on subplot 1
					subplot(2,1,1);
					hold on
					plot(timeFitData{1},[voltFitData{1} voltModel])
                    % plot(timeFitData{1},)% Plot meaused voltage vs model voltage
					legend('Data','Model')
					ylabel('Voltage (V)')
									
					if i_f == 4
						pltDataR = zeros(5,2);
						pltDataTau = zeros(5,2);
					else
						pltDataR = zeros(5,1);
						pltDataTau = zeros(5,1);						
					end

					% pltDataR contains the R values, pltDataTau contains
					% the tau values
					switch i_f
						case 2
							pltDataR([1 2],1) = xOpt{i_f,i_c}(1:2);
							pltDataTau([3],1) =  xOpt{i_f,i_c}(3);
						case 3
							pltDataR([1 2 4],1) =  xOpt{i_f,i_c}([1 2 4]);
							pltDataTau([3 5],1) =  xOpt{i_f,i_c}([3 5]);
						case 4
							pltDataR([1 2 4],:) =  xOpt{i_f,i_c}([1 2 4],:);
							pltDataTau([3 5],:) =  xOpt{i_f,i_c}([3 5],:);
					end
					
					% Plot Parameters on subplot 2
					subplot(2,1,2);
					figTitle = ['f = ' num2str(fOpt(i_f,i_c))];
					
					hold on
					bar(1:length(pltDataR),pltDataR) % Plot resistance values
					ylabel('Resistance (Ω)')
					ylim([0 0.1])
					
					yyaxis right		% plot taus on different axis
					bar(1:5,pltDataTau)
					set(gca,'yscale','log') % logarithmic y scale for tau
					ylim([1 500])		
					ylabel('τ (s)')
					
					% Label X xomponents
					xticks([1:5])
					xticklabels({'R_0','R_1','\tau_1','R_2','\tau_2'})
					
					pltTxt = pltDataR + pltDataTau; % Data used for text on the plot
					pltTxt(pltTxt == 0) = NaN; % Set 0 entries to NaN (to make it clear that these parameters were not fitted)
					
					yyaxis left
					% Plot text on bar chart
					text((1:length(pltTxt)) - 0.35,0.09*ones(1,length(pltTxt)),num2str(pltTxt(:,1),'%1.3g'),'vert','middle','horiz','center','rotation',90);
					if i_f == 4 &&  EcnSettings.fitSocDependence % Plot text for parameters at second SoC breakpoint
						text((1:length(pltTxt)) + 0.35,0.09*ones(1,length(pltTxt)),num2str(pltTxt(:,2),'%1.3g'),'vert','middle','horiz','center','rotation',90);
					end
					title(figTitle)
					
					
					% Save figure of last fit iteration
					if 	i_f == max(fitIterations)
						saveFolder = ['Intermediate_Results\Figures\File' num2str(currentFile)];
						
						if ~isfolder(saveFolder)
							mkdir(saveFolder)
						end
						
						if EcnSettings.saveDebugPlotAfterFit
							saveas(gcf,[saveFolder '\Cycle_' num2str(i_c) '.png']);
						end
					end
				end
			end
			i_f = i_f + 1;
		end
	end
end

%% 2 Store results in result struct
Results.f.optAf = fOpt;
Results.x.optAf = xOpt;

fprintf(['\nTotal time to fit file: ' num2str(round(toc(overallTimer)/3600,2)) ' hrs.\n']);

end






