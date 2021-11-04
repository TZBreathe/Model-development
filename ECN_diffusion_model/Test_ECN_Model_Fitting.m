 % Confidential
 % 
 % -SYNOPSIS
 % Test different ECN models by fitting to a single param run. DCIR data
 % already segmented. DCIR File No. 11 or 12.
 %
 

%% 1. Set up Paranmeterisation
clear variables
tic

% Use the diary function to record the output to the command line
diaryDateStr = datestr(now, 'yymmdd_HHMMSS');
diaryName = (['Diary\' diaryDateStr '.log']);

% Make folders for diary and intermediate results
if ~isfolder('Diary')
	mkdir('Diary');
end
diary (diaryName) % Start recording in diary

if ~isfolder('Intermediate_Results')
	mkdir('Intermediate_Results');
end

%% 2. Parameterisation Settings

EcnSettings.fastTestRun = false; % true: makes optimisation very fast to see if everything is working. false: normal optimisation
EcnSettings.debugPlot = true; % true: plot some debugging figures (optimisation progress and fit at every iteration) false; do not plot anything
EcnSettings.saveDebugPlotAfterFit = true; % true: save figure of fit at every iteration, false: do not save plot

EcnSettings.fitSocDependence = true; % true: interpolate SoC in fit for final fit iteration. false: do not interpolate SoC (use only SoC from beginning of cycle)
EcnSettings.skipLastIteration = true; % true: Skip final iteration (3 iterations total). false: do not skip final iteration (4 iterations total)

EcnSettings.fitFastRelax = false; %true: Fit relaxation of fast pulse. false: do not fit relaxation of fast pulse
EcnSettings.fitSlowRelax = false; %true: Fit relaxation of slow pulse. false: do not fit relaxation of slow pulse

% For the constrain settings, the upper bound is:
% initialGuess*constrainFac, 
% and the lower bound is:
% initialGuess/constrainFac
EcnSettings.constrainFacHfEstimates = [1 1.5 2.5 2]; % Factor by wich HF parameters are constrained for each fit iteration
EcnSettings.constrainFacLfEstimates = [NaN NaN 2 2.5]; % Factor by wich LF parameters are constrained for each fit iteration

% Expexted maximum value of objective function, if the true value is lower,
% an iteration is repeated
EcnSettings.fThreshold = [0 70 70 50];

EcnSettings.initPopSpread = [0 2 2 2]; % Number governing the initial solution population spread (higher = less spread, minimum = 1)
EcnSettings.optiTimeMins = [1 3 5 6]; % Time in minutes for each GA optimisation iteration

EcnSettings.paramLowerBounds = [0.001	2e-3	0.5		1e3      0.8]'; % R_0, R_1, tau_1, tauD, kD
EcnSettings.paramUpperBounds = [0.05	0.3	     50		2e3		1.5]';

% Breakpoints to use in the ECN LUT
Ecn.Dims.soc = 0:0.05:1;
Ecn.Dims.temperature = [-10 0 10 25 40 55 70];
Ecn.Dims.crate = [-50 -30 -10 10 20]/4.8;

tempBpTolerance = 0; %(°C) This value determines how far away the acutual temperature 
% (from the file name) can be from the SoC LUT dimension in <Ecn.Dims.temperature>

crateBpTolerance = 0.2; %(hr^-1) This value determines how far away the acutual C-rate 
% (in the measurement) can be from the C-rate LUT dimension in <Ecn.Dims.crate>
% This can vary a bit because the cell capacity varies

socBpTolerance = 5e-3; %(0 to 1) This value determines how far away the acutual SoC
% (in the measurement) can be from the SoC LUT dimension in <Ecn.Dims.soc>
% This can vary a bit because of cycler inaccuracy

% This variable determines which cycles are used for the fitting.
% Setting the value to 0 fits all cycles where the starting SoC is a member of Ecn.Dims.soc
fittingCycles = [0];

%% 3. Load data for parameterisation
% Load Data from File_List

fileTable = readtable('DCIR_File_List.xlsx');
fileIterator = 11; % Choose the single param test at 25 deg 2C rate.

% load OCV data from ECN_Parameters (the other data, such as ECN, is not used)
load('ECN_Parameters.mat','S48X')

OCV = S48X.OCV;

% Initialise Ecn results container
ecnComponentNames = {'R_0','R_1','tau_1','R_2','tau_2'};
for i_c = 1:length(ecnComponentNames)
	Ecn.Components.(ecnComponentNames{i_c}) = zeros(length(Ecn.Dims.soc),length(Ecn.Dims.temperature),length(Ecn.Dims.crate));
end

%% 4. Run Parameterisation


	currentFile = fileIterator;
	
	% Print some info on fit and expected duration
%  	fprintf(['\n--------------------------------------------------------------------------------------------------------------------------'...
% 		'\nFitting file ' num2str(currentFile) ', Number ' num2str(currentFile) ' of ' num2str(length(fileIterator)) ...
% 		'. Current: ' num2str(fileTable.Current(currentFile)) ' A, Temperature: ' num2str(fileTable.Temperature(currentFile)) ...
% 		' °C. Elapsed time: ' num2str(toc/3600,'%.3g') ' hrs, remaining Time: ' num2str((toc/(i_d/(1+length(fileIterator)))-toc)./3600,'%.3g') ' hrs.\n'...
% 		'--------------------------------------------------------------------------------------------------------------------------']);
	
	loadFromIntermediates = false; % Set to true to use data stored in \Intermedate_Results\ instead of running the optimisation
	if loadFromIntermediates
		load(['Intermediate_Results\File_No_' num2str(currentFile) '.mat'])
    end
		
		%% 4.1 Parameterisation Set-up
		
		% Load Data from DCIR_Data
		load(['DCIR_Data\File_No_' num2str(currentFile) '.mat'],'dataTable');
		
		% Store OCV soc LUT
		OcvLuts.soc = OCV.Dims.soc;
		
		tempBp = fileTable.Temperature(currentFile); % Get temperature breakpoint value
		tempIdx = (OCV.Dims.temp == tempBp);			% Get temperature breakpoint index for OCV
		
		% If no value of tempIdx is true, use the nearest one that exists
		% (for example because OCV data is not available at that temperature)
		if ~any(tempIdx)
			[~,tempIdx]=min(abs(OCV.Dims.temp - fileTable.Temperature(currentFile))); % Find nearest existing temperature breakpoint
			fprintf('\n') %This ensures the warning is on one line, and at the beginning
			warning(sprintf(['Temperature index ' num2str(fileTable.Temperature(currentFile)) ' °C does not exist in OCV table, using nearest: '...
				num2str(OCV.Dims.temp(tempIdx)) ' °C!'])) % Display warning that temperature breakpoint does not exist in the OCV
		end
		
		capCell = max(dataTable.capEst)/1000; % Get cell capacity (value is in mAh, convert to Ah)
		
		%% 4.2 Assign which step numbers correspond to which pulse (based on whether it is charging or discharging data)
		if fileTable.Current(currentFile) > 0 	% Charging 
			chargeDirection = 1;

			OcvLuts.ocv = OCV.Components.ocv(:,tempIdx) +  OCV.Components.hystAmp(:,tempIdx) + OCV.Components.hystInst(:,tempIdx);
			
			crateBp = fileTable.Current(currentFile)./capCell;
			
			StpInfo.fastPulseStp = 78;
			StpInfo.slowPulseStp = 80;
			StpInfo.dvdtStep = 84; %Index when voltage change is measured (used to only use the first hour of relaxation)
			lastPulseStep = 109; % index of last step. Data after this is ignored
			
		else			% Discharging 
			chargeDirection = -1;

			OcvLuts.ocv = OCV.Components.ocv(:,tempIdx) - OCV.Components.hystAmp(:,tempIdx) - OCV.Components.hystInst(:,tempIdx);

			crateBp = fileTable.Current(currentFile)./capCell;
			
			StpInfo.fastPulseStp = 71;
			StpInfo.slowPulseStp = 73;
			StpInfo.dvdtStep = 81; %Index when voltage change is measured (used to only use the first hour of relaxation)
			lastPulseStep = 102; % index of last step. Data after this is ignored
		end
		
		%% 4.3 Make the variable cycleCount correspond to one set of optimisation data (cycle) for the optimisation process
		% See documentation for further explanation
		
		% Take table entries as arrays, because this makes operations much
		% faster
		stepIdxArray = dataTable.stepIdx; 
		cycleCounterArray = dataTable.cycleCount;
		
		% Initialise output array
		cycleCounterArrayUpdate = zeros(length(cycleCounterArray),1);
		
		% Set counters to 0 at the beginning
		cycleCountNow = 0;
		cycleCountAdd = 0;
		
		% Loop through all time steps
		for i = 2:length(dataTable.cycleCount)-1
			idxNow = stepIdxArray(i);
			idxOld = stepIdxArray(i-1);
			
			if idxNow == StpInfo.fastPulseStp && idxOld ~= StpInfo.fastPulseStp %Find transition from rest to hf pulse
				cycleCountAdd = 1;
			end
			if cycleCounterArray(i) ~= 0
				cycleCountNow = cycleCounterArray(i);
			else
				cycleCounterArray(i) = cycleCountNow;
			end
			
			cycleCounterArrayUpdate(i) = cycleCounterArray(i) + cycleCountAdd;
			
			if cycleCounterArrayUpdate(i) == cycleCounterArray(i + 1)
				cycleCountAdd = 0;
			end
			
			if dataTable.stepIdx(i) >= lastPulseStep % For datapoints after the last pulse, set counter to 0
				cycleCounterArrayUpdate(i:end) = 0;
				break
			end
		end
		
		cycleCounterArrayUpdate = circshift(cycleCounterArrayUpdate,-50,1); % include some rest period (50 data points) from the previous step before the fast pulse
		dataTable.cycleCount = cycleCounterArrayUpdate; % Write updated counter array to dataTable
        
		%% 4.4 Run Parameterisation
		[cycleTabls, dataTable] = Prepare_ECN_Step_Data(dataTable,capCell,OcvLuts,chargeDirection); % This script divides the data of a test into cycles
		
		% Determine which cycles to fit (either all with an SoC in the ECN LUT, or user-defined cycles)
		if fittingCycles == 0 % This means that all cycles (that have an SoC in the ECN LUT) are fitted
			
			cyclesToFit = []; % Initialise array
			for i_c = 1:length(cycleTabls) % Loop over all cycles
				socBpNow = cycleTabls{i_c}.socOut(1);	% Get the SoC at the beginning of this cycle
				if any(abs(Ecn.Dims.soc - socBpNow) < socBpTolerance) % If this condition is true, the SoC at the beginning of a pulse corresponds to an SoC in the LUT...
					%... so the parameter fitting is done for this cycle...
					cyclesToFit(end+1) = i_c;
				else %... if it is not true, do not fit this cycle or any subsequent cycles (voltage limit was reached and the SoC values are no longer correct).
					break
				end
			end
			EcnSettings.cycles = cyclesToFit;
		else % If fittingCycles is not zero, the user can choose which cycles to fit
			EcnSettings.cycles = fittingCycles;			
		end
		
		[RunResults, socBpsTest] = ECN_Paramsation(cycleTabls,EcnSettings,StpInfo,currentFile,OcvLuts); % Function that runs the main parameterisation
			
 		save(['Intermediate_Results\File_No_' num2str(currentFile) '.mat'],'-v7.3')
		
	
	
	%% 4.5 Save parameters to ECN
% 	clear saveComponents
% 	
	if fittingCycles == 0
		
		socBpsSave = false(1,length(Ecn.Dims.soc)); % initialise array of indices for which parameters are saved
		socEcnIdx = false(1,length(Ecn.Dims.soc)); % initialise array of indices for which parameters are saved
		
		for i_bp = 1:length(RunResults.x.optAf)	% This loops over all SoC breakpoints that were parameterised
			
			% Choose which SoC breakpoint data is used to save the parameters
			% 1 means SoC at the beginning of a cycle
			% 2 means SoC at the end of a cycle
			socResultsIdx = 1;
			
			socEcnIdx = logical(abs(Ecn.Dims.soc - socBpsTest(i_bp)) < socBpTolerance);
			
			% Store optimisation results (x values) in saveComponents array,
			% which contains R, C, and tau values
			saveComponents.R_0(socEcnIdx) =  RunResults.x.optAf{end,i_bp}(1,socResultsIdx);
			saveComponents.R_1(socEcnIdx) = RunResults.x.optAf{end,i_bp}(2,socResultsIdx);
			saveComponents.C_1(socEcnIdx) = RunResults.x.optAf{end,i_bp}(3,socResultsIdx)./RunResults.x.optAf{end,i_bp}(2,socResultsIdx);
			saveComponents.R_2(socEcnIdx) = RunResults.x.optAf{end,i_bp}(4,socResultsIdx);
			saveComponents.C_2(socEcnIdx) = RunResults.x.optAf{end,i_bp}(5,socResultsIdx)./RunResults.x.optAf{end,i_bp}(4,socResultsIdx);
			saveComponents.tau_1(socEcnIdx) = RunResults.x.optAf{end,i_bp}(3,socResultsIdx);
			saveComponents.tau_2(socEcnIdx) = RunResults.x.optAf{end,i_bp}(5,socResultsIdx);
			
			socBpsSave(socEcnIdx) = 1;
		end
		
		tempBpSave = (abs(Ecn.Dims.temperature - tempBp) <= tempBpTolerance); % Find index of temperature breakpoint (with tolerance)
		crateBpSave = (abs(Ecn.Dims.crate - crateBp) <= crateBpTolerance);	% Find index of C-rate breakpoint (with tolerance)
		
		fprintf(['\nSoC breakpoints saved: ' sprintf('%d%%%%, ',round(Ecn.Dims.soc(socBpsSave)*100)) '\b\b.\n']) % Display which SoC breakpoints are saved
		
		ecnFieldNames = fieldnames(saveComponents);
		for i_f = 1:length(ecnFieldNames) % Loop over all Ecn.Components fieldnames
			Ecn.Components.(ecnFieldNames{i_f})(socBpsSave,tempBpSave,crateBpSave) = saveComponents.(ecnFieldNames{i_f})(socBpsSave);
		end
		
	end

%% 5 Save data and exit

% if fittingCycles == 0
% 	save('48X_ECN_Results.mat','Ecn') % Save ECN file at the end
% end
% 
% fprintf(['\n----------------------------------------------------------------\nOptimisation Finished!  Total time: '...
% 	num2str(round(toc/3600,2),'%.3g') ' hrs.\n----------------------------------------------------------------\n'])
% diary off % Turn off diary

