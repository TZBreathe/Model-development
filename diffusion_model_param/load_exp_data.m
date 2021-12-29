%Get experimental Lincc charge data
clear;
currentSeries={};
voltageSeries={};
socRefSeries={};
tempSeries=[];
%25 degrees
%max 4 cycles
Br{1}=gdFun.Load_Breathe_Run(456); % 12.5A, 20mV
Br{2}=gdFun.Load_Breathe_Run(459); % 12.5 A, 5mV
Br{3}=gdFun.Load_Breathe_Run(476); % 15A, 10mV
Br{4}=gdFun.Load_Breathe_Run(477); % 15A, 0mV
Br{5}=gdFun.Load_Breathe_Run(871); % 15A, -10mV
%min cycle 5
Br{6}=gdFun.Load_Breathe_Run(1067); % 12.5A, 10mV
Br{7}=gdFun.Load_Breathe_Run(1102); % 12.5A, -20mV
Br{8}=gdFun.Load_Breathe_Run(1071); % 20 A, -20mV
Br{9}=gdFun.Load_Breathe_Run(1103); % 20A, -20mV
Br{10}=gdFun.Load_Breathe_Run(1112); % 26A, -20mV

%35 degrees
Br{11}=gdFun.Load_Breathe_Run(1073); % 12.5A, 10mV
Br{12}=gdFun.Load_Breathe_Run(891); % 15A, 0mV, only 4 cycles
Br{13}=gdFun.Load_Breathe_Run(1072); % 20A, -20mV
Br{14}=gdFun.Load_Breathe_Run(1206); % 30A, -20mV


%45 degrees
Br{15}=gdFun.Load_Breathe_Run(1069); % 12.5A, 10mV, issues
Br{16}=gdFun.Load_Breathe_Run(1109); % 26A, 10mV 
Br{17}=gdFun.Load_Breathe_Run(1068); % 20A, -20mV, issues
Br{18}=gdFun.Load_Breathe_Run(1205); % 30A, -20mV


%% Get charge voltage 

for i=1:5
 Idx=find((Br{i}.RunData.dataTable{:,'cycleNumber'}==2) & (Br{i}.RunData.dataTable{:,'currCell'}>0));
 currentSeries{i}= Br{i}.RunData.dataTable{Idx,'currCell'};
 voltageSeries{i}= Br{i}.RunData.dataTable{Idx,'voltCell'};
 socRefSeries{i}= Br{i}.RunData.dataTable{Idx,'socOut'};
 tempSeries{i}=Br{i}.RunData.dataTable{Idx,'tempCell'};
end

for i=6:18
 Idx=find((Br{i}.RunData.dataTable{:,'cycleNumber'}==6) & (Br{i}.RunData.dataTable{:,'currCell'}>0));
 currentSeries{i}= Br{i}.RunData.dataTable{Idx,'currCell'};
 voltageSeries{i}= Br{i}.RunData.dataTable{Idx,'voltCell'};
 socRefSeries{i}= Br{i}.RunData.dataTable{Idx,'socOut'};
  tempSeries{i}=Br{i}.RunData.dataTable{Idx,'tempCell'};
end

for i=12
 Idx=find((Br{i}.RunData.dataTable{:,'cycleNumber'}==2) & (Br{i}.RunData.dataTable{:,'currCell'}>0));
 currentSeries{i}= Br{i}.RunData.dataTable{Idx,'currCell'};
 voltageSeries{i}= Br{i}.RunData.dataTable{Idx,'voltCell'};
 socRefSeries{i}= Br{i}.RunData.dataTable{Idx,'socOut'};
  tempSeries{i}=Br{i}.RunData.dataTable{Idx,'tempCell'};
end

%%
for i=1:10
lincc_25{i}=[currentSeries{i} voltageSeries{i} socRefSeries{i} tempSeries{i}];
end
save 'lincc_25.mat' lincc_25;

for i=1:4
lincc_35{i}=[currentSeries{i+10} voltageSeries{i+10} socRefSeries{i+10} tempSeries{i+10}];
end
save 'lincc_35.mat' lincc_35;

for i=1:4
lincc_45{i}=[currentSeries{i+14} voltageSeries{i+14} socRefSeries{i+14} tempSeries{i+14}];
end
save 'lincc_45.mat' lincc_45;


%% temp compare discharge V

clear;
currentSeries={};
voltageSeries={};
socRefSeries={};
tempSeries=[];

Br{1}=gdFun.Load_Breathe_Run(456); % 12.5A, 20mV
Br{2}=gdFun.Load_Breathe_Run(459); % 12.5 A, 5mV
Br{3}=gdFun.Load_Breathe_Run(476); % 15A, 10mV
Br{4}=gdFun.Load_Breathe_Run(477); % 15A, 0mV
Br{5}=gdFun.Load_Breathe_Run(1225); % 12.5A, 20mV
for i=1:4
 Idx=find((Br{i}.RunData.dataTable{:,'cycleNumber'}==2) & (Br{i}.RunData.dataTable{:,'currCell'}>0));
 currentSeries{i}= Br{i}.RunData.dataTable{Idx,'currCell'};
 voltageSeries{i}= Br{i}.RunData.dataTable{Idx,'voltCell'};
%  socRefSeries{i}= Br{i}.RunData.dataTable{Idx,'socOut'};
%  tempSeries{i}=Br{i}.RunData.dataTable{Idx,'tempCell'};
end





hold on
for i=1:2
    plot(voltageSeries{i},'DisplayName',strcat('No=',num2str(i)));
   
end
legend('show');
hold off

figure
hold on
for i=3:4
    plot(voltageSeries{i},'DisplayName',strcat('No=',num2str(i)));
   
end
legend('show');
hold off