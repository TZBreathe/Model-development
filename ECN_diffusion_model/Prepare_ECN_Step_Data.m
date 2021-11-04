 % Confidential
 % 
 % -SYNOPSIS
 % This function prepares the DCIR test data by splitting it into cycles (each cycle is an SoC breakpoint)
 %
 % -NOTES
 %    Version:      1.0
 %    Author:       Breathe Battery Technologies
 %    E-mail:       christian.korte@breathe.technology
 %    File: 		Prepare_ECN_Step_Data.m
 %    Copyright (c) 2021 Breathe Battery Technologies Ltd.
 % -LINK
 %    https://www.breathe.technology/
 % 
 % -LICENSE
 % Copyright Breathe Battery Technologies Ltd. May 2021 - All Rights Reserved
 % Provided under licence to Rimac Automobili d.o.o.

function [cycleTabs, dataTable] = Prepare_ECN_Step_Data(dataTable,capCell,OcvLuts,chargeDirection)
%% Extract Cycle Data

if chargeDirection == 1 % Charging dataset
   	dataTable.ahTot = dataTable.ahChrge;
elseif chargeDirection == -1 % Discharging dataset
   	dataTable.ahTot = capCell - dataTable.ahDchrge;
end

% Calculate SoC at each point in dataTable
dataTable.socOut = dataTable.ahTot/capCell;

% Check if SoC is within a reasonable range (sometimes it is slightly lower than 0 or higher than 1)
if any(dataTable.socOut > 1.03) || any(dataTable.socOut < -0.04)
	error('SoC out of range!!')
end

% Generate OCV data by interpolating OCV LUT
dataTable.ocpotCell = interp1(OcvLuts.soc,OcvLuts.ocv,dataTable.socOut,'linear','extrap');

% Split data into cycles (using cycleCount)
% Each cycle corresponds to one SoC breakpoint in the ECN LUT
for i_c = 1:(max(dataTable.cycleCount))
	% Each cycle has its own table of data for only this cycle
	cycleTabs{i_c} = dataTable(dataTable.cycleCount == (i_c),:);
end

end

