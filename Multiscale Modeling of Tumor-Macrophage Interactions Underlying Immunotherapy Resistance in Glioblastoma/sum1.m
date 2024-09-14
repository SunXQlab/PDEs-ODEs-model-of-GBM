function sum_1

% Written by Xiaoqiang Sun, Haofeng Lin   Apr 05, 2021

addpath('funfun'); 

global eta1    %Degradation Rate of Drugs
global D_d1    %Diffusion Rate of Drugs
global delt_t  %Time Interval in Simulation
global L  L2    %Size of the x-axis during simulation
global B_3  
global A_0   
global d_max         
global period        %This is for boundary drug administration with a cycle
global inter         %Transition interval of the step function
global period1 period_2      %Real-time period
global inter1        %Actual time interval
global TimeLength    %Total duration in days
global m0        %Amount of drug designed for boundary condition three
global Ptest1
global EGF IGF1
global EGFR_I IGF1R_I
global C_T_0 C_M1_0 C_M2_0 CSF1_0 EGF_0 IGF1_0
global ERK AKT
global x_0
global V7 a V81 V82 V9 V10 d7 d8 d9 d10 K71 K72 K81 K82 K83 K91 K92 K101 K102 K73 K93 n
global A_T_EGF A_T_IGF1 A_M1 A_M2 epsilon r_T a_ERK K1 K2 K31 K32 K4 a_AKT  d_T d_TM1 a_M2M1 a_CSF1  a_A K3 B_1 B_2 b_EGF
global CSF1_max EGF_max IGF1_max dmax1 dmax2 dmax3  
global select_CSF1R_I select_EGFR_I select_IGF1R_I 
Parameters_nondimensional_5_gai_jin_fangcheng_5
load('canshu1-24hours-ode15s.mat','p');
Ptest1=p(1,1:22);

TimeLength=400;
L2=100;
L=1/x_0 ;  %Total length of the axis during simulation
x = 0:L/L2:L;            %Spatial growth of the tumor
delt_t=1/91.5926;      %Time interval during simulation
t=0:delt_t:delt_t*TimeLength;
period_2=delt_t*14; %sin
period=7;           %This is for boundary drug administration with a cycle
period1=period*delt_t;  %Actual half-cycle
T=period1*2; %Actual cycle
inter=1;             
inter1=inter*delt_t;   
eta1=eta_d;
D_d1=D_d;


select_CSF1R_I=1;%When selecting adaptive case 5, PDEPE3 needs to be used for solving
select_EGFR_I=0;
select_IGF1R_I=0;

dmax1=1;
dmax2=1;
dmax3=1;

global b t5 b00 c00 fin1 fin2 tumor_mean_min
b=[1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];
b=b*1;
t50=10000*ones(40,1);
t5=t50*delt_t;
t5(1)=0;
b00=1;
c00=1;
fin1=0.3;     %Adaptive drug dosage threshold
fin2=0.05;

tumor_mean_min=0;
options=odeset('reltol',1e-20);%Precision function cannot be omitted
m = 2;%spherical
sol = pdepe2(m,@pdex4pde,@odex4ode,@pdex4ic,@pdex4bc,x,t,options);
% sol = pdepe3(m,@pdex4pde,@odex4ode,@pdex4ic,@pdex4bc,x,t,options);%%adaptive therapy
u1 = sol(:,:,1);
u2 = sol(:,:,2);
u3 = sol(:,:,3);
u4 = sol(:,:,4);
u5 = sol(:,:,5);
u6 = sol(:,:,6);
u7 = sol(:,:,7);
u8 = sol(:,:,8);
u9 = sol(:,:,9);
u10 = sol(:,:,10);
C_T_mean=sum(u1(:,:),2)/100; 
M1_mean=sum(u2(:,:),2)/100; 
M2_mean=sum(u3(:,:),2)/100; 

%   TimeLength=400
%   A1=zeros(TimeLength+1,L2+1);
%   A2=zeros(TimeLength+1,L2+1);
%   A3=zeros(TimeLength+1,L2+1);
%   A4=zeros(TimeLength+1,L2+1);
%   A5=zeros(TimeLength+1,L2+1);
%   A6=zeros(TimeLength+1,L2+1);
%   A7=zeros(TimeLength+1,L2+1);
%   A8=zeros(TimeLength+1,L2+1);
%   A9=zeros(TimeLength+1,L2+1);
%   A10=zeros(TimeLength+1,L2+1);
% 
% %   
%   x00=0:1:L2;
%   t00=0:1:TimeLength;
% 
% %   TimeLength=200;
% % %   
% %   L=L*x_0
% global b1 b2 b3 b4
% global t51 t52 t53 t54 t55 t56 t57 t58
% b1=1;b2=0;b3=1;b4=0;
% t51=0*delt_t;t52=99*delt_t;t53=100*delt_t;t54=300*delt_t;t55=100*delt_t;t56=149*delt_t;t57=150*delt_t;t58=201*delt_t;
% global b t5
% % b=[1 0.5 1 0.5 1 0.5 1 0.5 1 0.5 1 0.5 1 0.5 1 0.5 1 0.5 1 0.5 1 0.5 1 0.5];
% b=[1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];
% t50=[0  25  26  60  61  100 101 150 151 200 201 300];
% t5=t50*delt_t;
% Drug2(0.1,49.5*delt_t,4,1)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
% % d_max=1;
% m0=1;
% period=7; 
% period1=period*delt_t;
% period_2=delt_t*14;
% period=30;           
% period1=period*delt_t;  
% T=period1*2; 
%   for i= 0:1:TimeLength
%       for j=0:1:L2
%           i1=i*delt_t;
%           j1=L/L2*j;
%           if j1==0
%               j1=0.0000001;
%           end
% %           drugt1=Drug2(j1,i1,1,1);
%           drugt2=Drug2(j1,i1,2,1);
% %           drugt3=Drug2(j1,i1,3,1);
% %           drugt4=Drug2(j1,i1,4,1);
% %           drugt5=Drug2(j1,i1,5,1);
% %           IL411=A_Drugsimulation2(j1,i1,1,1);
% %           IL412=A_Drugsimulation2(j1,i1,2,1);
% %           IL413=A_Drugsimulation2(j1,i1,3,1);
% %           IL414=A_Drugsimulation2(j1,i1,4,1);
% %           IL415=A_Drugsimulation2(j1,i1,5,1);
% %           A1(i+1,j+1)=drugt1;
%           A2(i+1,j+1)=drugt2;
% %           A3(i+1,j+1)=drugt3;
% %           A4(i+1,j+1)=drugt4;
% %           A5(i+1,j+1)=drugt5;
% %           A6(i+1,j+1)=IL411;
% %           A7(i+1,j+1)=IL412;
% %           A8(i+1,j+1)=IL413;
% %           A9(i+1,j+1)=IL414;
% %           A10(i+1,j+1)=IL415;
%       end
%   end

figure
pcolor(x,t,u1)
% title('Tumor')
% title('Cancer Cell Density','FontWeight','Bold','FontSize',18);
title('Cancer cell density','FontWeight','Bold','FontSize',18,'FontName','Arial');
% title('Continuous Drug','FontWeight','Bold','FontSize',18);
% colorbar('ylim',[0.3,0.7],'ytick',(0.3:0.1:0.7));
% colorbar('Ticks',0:0.1:1);
colorbar;
shading interp;
   set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength,'FontName','Arial');
   set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
   set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
   set(gca,'FontSize',16,'FontName','Arial');



figure
pcolor(x,t,u2)
title('M1 Macrophage density','FontWeight','Bold','FontSize',18,'FontName','Arial');
% colorbar('ylim',[0.3,1],'ytick',(0.3:0.1:1));
colorbar;
shading interp;
   set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength,'FontName','Arial');
   set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
   set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
   set(gca,'FontSize',16,'FontName','Arial');



figure
pcolor(x,t,u3)
title('M2 Macrophage density','FontWeight','Bold','FontSize',18,'FontName','Arial');
colorbar;
shading interp;
   set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength,'FontName','Arial');
   set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
   set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
   set(gca,'FontSize',16,'FontName','Arial');



figure
pcolor(x,t,u4)
title('CSF1 Concentration','FontWeight','Bold','FontSize',18,'FontName','Arial');
% colorbar('ylim',[0.1,0.4],'ytick',(0.1:0.05:0.4)); 
colorbar;
shading interp;
   set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength,'FontName','Arial');
   set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
   set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
   set(gca,'FontSize',16,'FontName','Arial');



 figure
 pcolor(x,t,u5)
 title('EGF Concentration','FontWeight','Bold','FontSize',18,'FontName','Arial');
% colorbar('ylim',[0.3,0.6],'ytick',(0.3:0.05:0.6)); 
colorbar;
shading interp;
   set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength,'FontName','Arial');
   set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
   set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
   set(gca,'FontSize',16,'FontName','Arial');

 figure
 pcolor(x,t,u6)
 title('IGF1 Concentration','FontWeight','Bold','FontSize',18,'FontName','Arial');
% colorbar('ylim',[0,0.9],'ytick',(0:0.1:0.9));
colorbar;
shading interp;
   set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength,'FontName','Arial');
   set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
   set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
   set(gca,'FontSize',16,'FontName','Arial');

    
figure
 pcolor(x,t,u7)
 title('EGFR Concentration','FontWeight','Bold','FontSize',18,'FontName','Arial');
% colorbar('ylim',[0.7,0.78],'ytick',(0.7:0.01:0.78)); 
colorbar;
shading interp;
   set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength,'FontName','Arial');
   set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
   set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
   set(gca,'FontSize',16,'FontName','Arial');

    
figure
 % surf(x,t,u5)
 pcolor(x,t,u8)
 title('ERK Concentration','FontWeight','Bold','FontSize',18,'FontName','Arial');
% colorbar('ylim',[0.55,0.75],'ytick',(0.55:0.05:0.75));
colorbar;
shading interp;
   set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength,'FontName','Arial');
   set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
   set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
   set(gca,'FontSize',16,'FontName','Arial');

    
figure
 pcolor(x,t,u9)
 title('IGF1R Concentration','FontWeight','Bold','FontSize',18,'FontName','Arial');
% colorbar('ylim',[0,0.025],'ytick',(0:0.005:0.025)); 
colorbar;
shading interp;
   set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength,'FontName','Arial');
   set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
   set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
   set(gca,'FontSize',16,'FontName','Arial');

figure
 pcolor(x,t,u10)
 title('AKT Concentration','FontWeight','Bold','FontSize',18,'FontName','Arial');
% colorbar('ylim',[0,0.45],'ytick',(0:0.05:0.45)); 
colorbar;
shading interp;
   set(gca,'ytick',0:delt_t*TimeLength/5:delt_t*TimeLength); set(gca,'yticklabel',0:TimeLength/5:TimeLength,'FontName','Arial');
   set(gca,'xtick',0:L/5:L); set(gca,'xticklabel',0:20:100,'FontName','Arial');
%    xlabel('Position','FontWeight','Bold','FontSize',18);ylabel('Time (Days)','FontWeight','Bold','FontSize',18);
   xlabel('Position','FontName','Arial');ylabel('Time (Days)','FontName','Arial');
   set(gca,'FontWeight','Bold','FontSize',14,'FontName','Arial');
   set(gca,'FontSize',16,'FontName','Arial');


   figure
   plot(0:delt_t:delt_t*TimeLength,C_T_mean,'-','LineWidth',3);
   title('Average cancer cell density','FontWeight','Bold','FontSize',18,'FontName','Arial')
   set(gca,'xtick',(0:delt_t*TimeLength/5:delt_t*TimeLength));
   set(gca,'xticklabel',0:TimeLength/5:TimeLength,'FontName','Arial');
   set(gca,'xlim',[0 delt_t*TimeLength]);
   set(gca,'ylim',[0 1]);
   xlabel('Time (Days)','FontWeight','Bold','FontSize',18,'FontName','Arial');ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18,'FontName','Arial');
   set(gca,'FontWeight','Bold','FontSize',16,'FontName','Arial');  


figure
   plot(0:delt_t:delt_t*TimeLength,M1_mean,'-','LineWidth',3);
   title('Average macrophage density','FontWeight','Bold','FontSize',18)
   set(gca,'xtick',(0:delt_t*TimeLength/5:delt_t*TimeLength));
   set(gca,'xticklabel',0:TimeLength/5:TimeLength);
%    set(gca,'ytick',(0:0.1:0.8));
%    set(gca,'yticklabel',0:0.1:0.8)
   set(gca,'xlim',[0 delt_t*TimeLength]);
   set(gca,'ylim',[0 1]);
   xlabel('Time (Days)','FontWeight','Bold','FontSize',18);ylabel('Density','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',16);  
   hold on
   plot(0:delt_t:delt_t*TimeLength,M2_mean,'-','LineWidth',3);
%    title('Average cell density','FontWeight','Bold','FontSize',18)
   set(gca,'xtick',(0:delt_t*TimeLength/5:delt_t*TimeLength));
   set(gca,'xticklabel',0:TimeLength/5:TimeLength);
%    set(gca,'ytick',(0:0.1:0.8));
%    set(gca,'yticklabel',0:0.1:0.8)
   set(gca,'xlim',[0 delt_t*TimeLength]);
   set(gca,'ylim',[0 1]);
   xlabel('Time (Days)','FontWeight','Bold','FontSize',18);ylabel('Density','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',16);




end


