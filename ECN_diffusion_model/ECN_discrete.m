%% A discrete version of ECN model

function voltOut=ECN_discrete(k,currData,timeData,socData,calcOnlyHf,OCV)

N=10; % controls timestep
dt=1/N;
R_0=k(1);
R_1 = k(2);
C_1=k(3)/k(2);
if ~calcOnlyHf
R_2=k(4);
C_2=k(5)/k(4);
end

Vrc_1=0;
Vrc_2=0;
voltOut=0;

%if using finer timestep 'interploate' input current array to larger size
times=N*length(timeData);
for ii=2:length(currData)
  curr_Data(N*ii+1:N*ii+10)=currData(ii);
end
curr_Data(1:N)=currData(1);

for timestep = 1:times
    
    Vrc_1=R_1*(exp(-dt/R_1/C_1).*(Vrc_1/R_1)+(1-exp(-dt/R_1/C_1)).*curr_Data(timestep));

    IR0=R_0.*curr_Data(timestep);  
    v_Out(timestep)=Vrc_1+IR0; 
    
    
    
    if ~calcOnlyHf
    Vrc_2=R_2*(exp(-dt/R_2/C_2).*(Vrc_2/R_2)+(1-exp(-dt/R_2/C_2)).*curr_Data(timestep));
    v_Out(timestep)=Vrc_1+Vrc_2+IR0; 
    end
end
 
for i=2:length(timeData)
voltOut(i)=v_Out(N*(i-1)+1);
end
voltOut(1)=v_Out(1);
voltOut=voltOut';
% voltOut=v_Out;

end

