%% A simple diffusion+ECN model

function [voltOut soc Vdiff]=diffusion_model_run(k,currData,timeData,socData,OcvLuts)
%params

N=1; % controls timestep
dt=1/N;

SoC0=socData(1);
voltOut=ones(length(timeData),1);
Vrc=0;

R_0=k(1);
I_0 = k(2);
alpha=k(3);
tau0=k(4);
beta=k(5);
kd=k(6);

% diffusion settings
Q=4.7*3600; %capacity, should be an argument but lazy for the moment
R = 1; % particle radius [m]
Nr = 20; % number of "shells" radially
dR = R/Nr; % width of each "shell"
Sa = 4*pi*(R*(1:Nr)/Nr).^2; % outer surface area of each shell
dV = (4/3)*pi*((R*(1:Nr)/Nr).^3-(R*(0:Nr-1)/Nr).^3); % vol. of ea. shell
SoC = SoC0*ones(1,(Nr)); % concentration profile versus "r" dimension
SoCs = zeros(size(timeData)); % concentration at surface
SoCavg=SoC0*ones(size(timeData));
SoCs(1) = SoC0;

h(1)=0;
k_hyst=10;
hyst=OcvLuts.Components.hystAmp(:,5); %Hyst data parameters, fifth column for 25 deg.
hyst_0=OcvLuts.Components.hystInst(:,5);
% SoCr=ones(length(timeData),Nr); %internal SoC, maybe useful for gradient]

%if using finer timestep 'interploate' input current array to larger size
times=N*length(timeData);
curr_Data=currData;
% for ii=2:length(currData)
%   curr_Data(N*ii+1:N*ii+10)=currData(ii);
% end
% curr_Data(1:N)=currData(1);

% calc diffusion 
for timestep = 1:times
  

IR0=R_0.*curr_Data(timestep);     
Vbv=0.0256*2*asinh(currData(timestep)/(I_0*(socData(timestep)+alpha)*(1-socData(timestep)+alpha))); %BV-like overpotential, larger at high & low soc

M_hyst=interp1(OcvLuts.Dims.soc,hyst,socData(timestep));
M0=interp1(OcvLuts.Dims.soc,hyst_0,socData(timestep));
h=exp(-dt*abs(curr_Data(timestep))*k_hyst/(Q))*h+sign(curr_Data(timestep))*(1-exp(-dt*abs(curr_Data(timestep)*k_hyst/(Q))));
U_hyst=M_hyst.*h+sign(curr_Data(timestep)).*M0;

tau=tau0/(socData(timestep)+beta); 
flux = -1/tau*diff(SoC)/dR; % flux at surfaces between "bins"
M = flux.*Sa(1:end-1); % total SoC crossing surface between bins
SoC= SoC+ ([0 M] - [M 0])*dt./dV; % conc. change via diffusion
SoC(end) = SoC(end) + (curr_Data(timestep)/3/Q)*Sa(end)*dt/dV(end); % at boundary
SoCs(timestep) = min(1,SoC(end)); % surface soc
% SoCavg(timestep)=(SoC*dV')/(4/3*pi*R^3); % average soc
SoCavg(timestep)= socData(timestep);
OCVcell(timestep)=interp1(OcvLuts.Dims.soc,OcvLuts.Components.ocv(:,5),SoCavg(timestep));% ocv at avearge soc
OCVcell_surf(timestep)=interp1(OcvLuts.Dims.soc,OcvLuts.Components.ocv(:,5),SoCs(timestep)); %ocv at surface soc
Vdiff(timestep)=-(kd*(OCVcell(timestep)-OCVcell_surf(timestep)));
v_Out(timestep)=IR0+Vbv+Vdiff(timestep)+OCVcell(timestep)+U_hyst;
% v_Out(timestep)=U_hyst;

soc(timestep)=SoCs(timestep);
end

voltOut=v_Out';

% for i=2:length(timeData)
%     
% voltOut(i)=v_Out(N*(i-1)+1);
% end
% voltOut(1)=v_Out(1);


end

