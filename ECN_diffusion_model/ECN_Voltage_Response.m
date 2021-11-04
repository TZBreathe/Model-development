 % Confidential
 % 
 % -SYNOPSIS
 % This function calculates the voltage response of the ECN
 %
 % -NOTES
 %    Version:      1.0
 %    Author:       Breathe Battery Technologies
 %    E-mail:       christian.korte@breathe.technology
 %    File: 		ECN_Voltage_Response.m
 %    Copyright (c) 2021 Breathe Battery Technologies Ltd.
 % -LINK
 %    https://www.breathe.technology/
 % 
 % -LICENSE
 % Copyright Breathe Battery Technologies Ltd. May 2021 - All Rights Reserved
 % Provided under licence to Rimac Automobili d.o.o.

function voltOut = ECN_Voltage_Response(k,currData,timeData,socData,calcOnlyHf,OCVLuts)
% This function calculates the voltage response of the ECN
% The input parameters are:
% k: array or R/C components
% currData: current data to simulate
% time: time data to simulate
% socData: soc data to simulate
OCVLuts=[];
v0 = [0 0];

% Set the settings for the ODE solver
odeOpts = odeset('RelTol',1e-4,'AbsTol',1e-8,'MaxStep',5);

% Call the ODE function of the two RC pairs
[~,vOde] = ode23t(@(t,y)TWO_RC_ODE(t,y,k,currData,timeData,socData), timeData, v0,odeOpts);

% Initialise r_0
r_0 = k(1,1)*ones(length(socData),1);

% if the SoC changes, calculate the SoC at each point...
if length(unique(socData)) >1
	% Calulate r_0 at each point in time
	for i_t = 1:length(timeData)
		r_0(i_t) = k(1,1).*(1-socData(i_t)) + k(1,2).*socData(i_t);
	end
else %... otherwise use the initialised value of r_0
	
end

voltOp = vOde(:,1) + vOde(:,2) + r_0.*currData;

voltOut = voltOp;

	function dVoltdt = TWO_RC_ODE(t,y,k,currentData,tData,soc)
		
		dVoltdt = zeros(2,1);

		currNow = interp1(tData,currentData,t);
		socNow = interp1(tData,soc,t);

		r_1 = k(2,1).*(1-socNow) + k(2,2).*socNow;
		tau_1 = k(3,1).*(1-socNow) + k(3,2).*socNow;
		
		if ~calcOnlyHf
			r_2 = k(4,1).*(1-socNow) + k(4,2).*socNow ;
			tau_2 = k(5,1).*(1-socNow) + k(5,2).*socNow;
		end
		
		dVoltdt(1,:) = -y(1)./(tau_1) + currNow/(tau_1/r_1);
		if ~calcOnlyHf
			dVoltdt(2,:) = -y(2)./(tau_2) + currNow/(tau_2/r_2);
		end
    end


% 	function dVoltdt = TWO_RC_ODE_sqrt(t,y,k,currentData,tData,soc)
% 		
% 		dVoltdt = zeros(2,1);
% 
% 		currNow = interp1(tData,currentData,t);
% 		socNow = interp1(tData,soc,t);
% 
% 		r_1 = k(2,1).*(1-socNow) + k(2,2).*socNow;
% 		tau_1 = k(3,1).*(1-socNow) + k(3,2).*socNow;
% 		
% 		if ~calcOnlyHf
% 			r_2 = k(4,1).*(1-socNow) + k(4,2).*socNow ;
% 			tau_2 = k(5,1).*(1-socNow) + k(5,2).*socNow;
% 		end
% 		
% 		dVoltdt(1,:) = -y(1)./((tau_1)) + currNow./(tau_1/r_1);
% 		if ~calcOnlyHf
% 			dVoltdt(2,:) = -y(2)./(2*sqrt(t)*tau_2) + currNow/(2*sqrt(t)*tau_2/r_2);
% 		end
% 	end

end

