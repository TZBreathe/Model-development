 % Confidential
 % 
 % -SYNOPSIS
 % This function retrieves and sorts the DCIR test data from .csv files and
 % saves then as .mat files
 %
 % -NOTES
 %    Version:      1.0
 %    Author:       Breathe Battery Technologies
 %    E-mail:       christian.korte@breathe.technology
 %    File: 		Get_And_Sort_DCIR_Data.m
 %    Copyright (c) 2021 Breathe Battery Technologies Ltd.
 % -LINK
 %    https://www.breathe.technology/
 % 
 % -LICENSE
 % Copyright Breathe Battery Technologies Ltd. May 2021 - All Rights Reserved
 % Provided under licence to Rimac Automobili d.o.o.
 
function [] = Get_And_Sort_DCIR_Data()
tic

fprintf('Getting and sorting DCIR data...')

if ~isfolder('DCIR_Data')
	mkdir('DCIR_Data');
end

% Location of .xlsx data files as a cell vector
dcirFolders = {	'Tests\DCIR\DCIR_10°C\',...
	'Tests\DCIR\DCIR_25°C\',...
	'Tests\DCIR\DCIR_40°C\',...
	'Tests\DCIR\DCIR_55°C\',...
	'Tests\DCIR\DCIR_70°C\'};

fileListTable = array2table([]); % Make empty table to store overview of tests
i_t = 1; %Overall counter for DCIR test files

% Loop over all folders
for i_fol = 1:length(dcirFolders)
	
	fileNames = dir([dcirFolders{i_fol} '*.csv']); % Find all .csv files in the current folder
	
	for i_fil = 1:length(fileNames) % Loop over all files in the current folder
		
		matFileName = ['DCIR_Data\File_No_' num2str(i_t) '.mat']; % Assign at MATLAB file namt
		
		readDataTable = readtable([fileNames(i_fil).folder '\' fileNames(i_fil).name]);
		
		dataTable = array2table([]); %empty table of cycle data
		
		%% Convert data to readable names
		dataTable.timeTest = readDataTable.TotalTime_Seconds_;
		dataTable.stepIdx = readDataTable.Step;
		dataTable.currCell = readDataTable.Current_A_;
		dataTable.voltCell = readDataTable.Voltage_V_;
		dataTable.ahChrge = readDataTable.ChargeCapacity_mAh_/1000;
		dataTable.ahDchrge = readDataTable.DischargeCapacity_mAh_/1000;
		dataTable.tempCell = readDataTable.Cell1__C_;
		readDataTable.Real_capacity_Var13_(isnan(readDataTable.Real_capacity_Var13_)) = 0; % Set all NaN values to 0
		readDataTable.Number_of_pulses_Var4_(isnan(readDataTable.Number_of_pulses_Var4_)) = 0; % Set all NaN values to 0
		dataTable.capEst = readDataTable.Real_capacity_Var13_;
		dataTable.cycleCount = readDataTable.Number_of_pulses_Var4_;
		
		save(matFileName,'dataTable','-v7.3');
		
		fprintf(['\nElapsed time: ' num2str(round(toc/60,1)) ' min. Read and saved data ' num2str(i_fil) ' of ' num2str(length(fileNames)) ' in folder ' num2str(i_fol) ' of ' num2str(length(dcirFolders)) '.\n']);
		
		%% Make table with overview of file list and details of each file
		
		fileListTable.No(i_t) = i_t; % Number associated with file
		% Extract temperature from file name
		fileListTable.Temperature(i_t) = str2num(extractBetween(string(fileNames(i_fil).folder),'DCIR_','°C'));
		
		% Extract current from file name
		if contains(fileNames(i_fil).name,'CHG')
			fileListTable.Current(i_t) = str2num(extractBetween(string(fileNames(i_fil).name),'CHG ','A 48X'));		
		elseif contains(fileNames(i_fil).name,'DCH')
			fileListTable.Current(i_t) = -str2num(extractBetween(string(fileNames(i_fil).name),'DCH ','A 48X'));					
		else
			error(['File ' fileNames(i_fil).name ' is not named as expected. The name should begin with: DCIR XXX YYA 48X..., where XXX is either DCH or CHG and YY is the current amplitude.'])
		end
		
		%Save file details
		fileListTable.Excel_Name(i_t) = {fileNames(i_fil).name}; % Original Excel file name
		fileListTable.Path(i_t) = {fileNames(i_fil).folder};	% Path to original file
		fileListTable.MATLAB_File_Name(i_t) = {matFileName};	% Matlab file name
		
		i_t = i_t + 1;
	end
end

% Save overview of data to an Excel file. This is used as a reference in
% Make_ECN to obtain the C-rates and temperatures of each data files
writetable(fileListTable,'DCIR_File_List.xlsx','Sheet',1,'Range','A1','WriteVariableNames',true);

end
