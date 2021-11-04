

for i=1:5
Ecn.Components.R_0(:,:,i)=Ecn.Components.R_0(:,:,4);   
Ecn.Components.R_1(:,:,i)=Ecn.Components.R_1(:,:,4);
Ecn.Components.R_2(:,:,i)=Ecn.Components.R_2(:,:,4); 
Ecn.Components.tau_1(:,:,i)=Ecn.Components.tau_1(:,:,4);
Ecn.Components.tau_2(:,:,i)=Ecn.Components.tau_2(:,:,4); 


end

for j=1:7
   Ecn.Components.R_0(:,j,:)=Ecn.Components.R_0(:,4,:);   
   Ecn.Components.R_1(:,j,:)=Ecn.Components.R_1(:,4,:);   
   Ecn.Components.R_2(:,j,:)=Ecn.Components.R_2(:,4,:);   
   Ecn.Components.tau_1(:,j,:)=Ecn.Components.tau_1(:,4,:);   
   Ecn.Components.tau_2(:,j,:)=Ecn.Components.tau_2(:,4,:);   
 
   
end

Ecn.Components.R_0(Ecn.Components.R_0==0)=1e-5;
Ecn.Components.R_1(Ecn.Components.R_1==0)=1e-5;
Ecn.Components.tau_1(Ecn.Components.tau_1==0)=2;
Ecn.Components.R_2(Ecn.Components.R_2==0)=1000;
Ecn.Components.tau_2(Ecn.Components.tau_2==0)=1;
