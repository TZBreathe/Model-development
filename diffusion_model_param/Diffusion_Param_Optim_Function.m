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
 

function [obj, voltModel] = Diffusion_Param_Optim_Function(k,currData,timeData,socData,voltData,tempData,OcvLuts)
% This allows the user to choose which degrees of freedom to use (dependence on soc)
% k is the optimisation variable

% lengthParamVec = 5;



% Calculate the voltage resopnse
voltModel = diffusion_model(k,currData,timeData,socData,tempData,OcvLuts);

obj = 1e3.*sum(sum((voltModel - voltData).^2))./length(timeData); 
% Factor of 1e6 results in higher numbers, which means the optimisation does not stop prematurely
% Divide by length of timeData to normalise the objective value and compare
% across different data sets

end

