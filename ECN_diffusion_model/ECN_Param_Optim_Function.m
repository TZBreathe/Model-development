 % Confidential
 % 
 % -SYNOPSIS
 % This function provides the optimisation function for the ECN
 % parameterisation process
 %
 % -NOTES
 %    Version:      1.0
 %    Author:       Breathe Battery Technologies
 %    E-mail:       christian.korte@breathe.technology
 %    File: 		ECN_Param_Optim_Function.m
 %    Copyright (c) 2021 Breathe Battery Technologies Ltd.
 % -LINK
 %    https://www.breathe.technology/
 % 
 % -LICENSE
 % Copyright Breathe Battery Technologies Ltd. May 2021 - All Rights Reserved
 % Provided under licence to Rimac Automobili d.o.o.
 

function [obj, voltModel] = ECN_Param_Optim_Function(k,currData,timeData,socData,voltData,useSocDependence,fitRCs,OcvLuts)
% This allows the user to choose which degrees of freedom to use (dependence on soc)
% k is the optimisation variable

calcOnlyHf = false; %when false, the voltage model only calculates the response of the fast RC pair

if fitRCs == 1 %only fit R0, R1 and C1
	lengthParamVec = 3;
	calcOnlyHf = true;
else % Fit everything else
	lengthParamVec = 5;
end

expectedLengthK = lengthParamVec + useSocDependence*lengthParamVec;

% These asserts can be removed to increase code speed
assert(length(k) == expectedLengthK,['Length of k (' num2str(length(k)) ') is not consistent with input fitRCs. k should have length ' num2str(expectedLengthK)])
assert(size(k,1) == 1,'k is not a row vector')

% when fitting parameters as a function of SoC (useSocDependence = true)
% the optimisation variable has 10 parameters instead of 5. The first 5 are
% for the SoC breakpoint at the beginning of the data, the second five are
% for the SoC breakpoint at the end of the data
if useSocDependence
	kWide = [k(1:lengthParamVec)' k((lengthParamVec+1):end)'];
else
	kWide(:,1) = k'; % double for values at 2 SoCs
	kWide(:,2) = k'; % double for values at 2 SoCs
end

% Calculate the voltage resopnse
voltModel = diffusion_discrete(kWide,currData,timeData,socData,calcOnlyHf,OcvLuts);

obj = 1e6.*sum(sum((voltModel - voltData).^2))./length(timeData); 
% Factor of 1e6 results in higher numbers, which means the optimisation does not stop prematurely
% Divide by length of timeData to normalise the objective value and compare
% across different data sets

end

