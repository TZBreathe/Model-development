%Get experimental CC charge data
clear;

currentSeries={};
voltageSeries={};
socRefSeries={};
%25 degrees
Br{1}=gdFun.Load_Breathe_Run(1225);
capacity=max(Br{1}.RunData.dataTable{:,'ahTotal'});
%% 
% Get previous disch capacity for soc estimation
for i=1:7
   Idx=find((Br{1}.RunData.dataTable{:,'stepIdx'}==(2*i-1))); 
   cap_disch(i)=Br{1}.RunData.dataTable{Idx(end),'ahTotal'};
end
cap_disch(1)=4.734; %value from last conditioning cycle
  
for i=1:7
 Idx=find((Br{1}.RunData.dataTable{:,'stepIdx'}==2*i) & (Br{1}.RunData.dataTable{:,'currCell'}>0));
 currentSeries{i}= Br{1}.RunData.dataTable{Idx,'currCell'};
 voltageSeries{i}= Br{1}.RunData.dataTable{Idx,'voltCell'};
 socRefSeries{i}= Br{1}.RunData.dataTable{Idx,'ahTotal'}./capacity+(capacity-cap_disch(i))./capacity;
end


