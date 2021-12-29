%Get experimental CC charge data
clear;
load brOCV;

currentSeries={};
voltageSeries={};
socRefSeries={};
tempSeries=[];
%25 degrees
Br{1}=gdFun.Load_Breathe_Run(1225);
capacity=max(Br{1}.RunData.dataTable{:,'ahTotal'});
%% 
% Get previous disch capacity for soc estimation
for i=1:7
   Idx=find((Br{1}.RunData.dataTable{:,'stepIdx'}==(2*i-1))); 
   cap_disch(i)=Br{1}.RunData.dataTable{Idx(end),'ahTotal'};
end
cap_disch(1)=capacity; %assume the first cycle starts with 0 soc (indicated by inverse ocv lookup)

for i=1:7
 Idx=find((Br{1}.RunData.dataTable{:,'stepIdx'}==2*i) & (Br{1}.RunData.dataTable{:,'currCell'}>0));
 currentSeries{i}= Br{1}.RunData.dataTable{Idx,'currCell'};
 voltageSeries{i}= Br{1}.RunData.dataTable{Idx,'voltCell'};
 socRefSeries{i}= Br{1}.RunData.dataTable{Idx,'ahTotal'}./capacity+(capacity-cap_disch(i))./capacity;
 tempSeries{i}=Br{1}.RunData.dataTable{Idx,'tempCell'}-4.5;
end


%% temp
hold on
plot(socRefSeries{5},voltageSeries{5});
plot(BrOcv.Dims.soc,BrOcv.Components.ocv (:,5));
hold off
figure

OCVS=interp1(BrOcv.Dims.soc,BrOcv.Components.ocv (:,5),socRefSeries{5});
plot(socRefSeries{5},voltageSeries{5}-OCVS);
