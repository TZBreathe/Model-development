
currData=table2array(cycleTabls{1, 10}(:,'currCell'));
timeData=table2array(cycleTabls{1, 10}(:,'timeTest'));
timeData=timeData-timeData(1);
socData=table2array(cycleTabls{1, 10}(:,'socOut'));


% currData=[zeros(1,100) 20*ones(1,50) zeros(1,500)]';
% timeData=linspace(1,length(currData),length(currData))';
% socData=ones(1,length(currData))';


calcOnlyHf=0;
k=[0.015 0.003 2 0.0054 40]';
k2=[k k];

V_disc=ECN_discrete(k,currData,timeData,socData,calcOnlyHf);
figure();
plot(timeData,V_disc);hold on;

V_ode=ECN_Voltage_Response(k2,currData,timeData,socData,calcOnlyHf);
% figure();
plot(timeData,V_ode);
% plot(timeData,(V_ode'-V_disc)/0.4);
%% random test

dt=0.1;
R_1=0.003; C_1=2/0.003;
Vrc_1=zeros(1,100);
Vrc_2=zeros(1,10);
  for ii=2:100
    Vrc_1(ii)=R_1*(exp(-dt/R_1/C_1)*(Vrc_1(ii-1)/R_1)+(1-exp(-dt/R_1/C_1))*1);
  end
 
  
  for i=1:1:10
  V1(i)=Vrc_1(10*i);
  end
 V1(1)=0;
 
   hold on; plot(V1,'o');
  
  dt2=1;
  for ii=2:10
  Vrc_2(ii)=R_1*(exp(-dt2/R_1/C_1).*(Vrc_2(ii-1)/R_1)+(1-exp(-dt2/R_1/C_1)).*1);
  end
 plot(Vrc_2,'s'); 