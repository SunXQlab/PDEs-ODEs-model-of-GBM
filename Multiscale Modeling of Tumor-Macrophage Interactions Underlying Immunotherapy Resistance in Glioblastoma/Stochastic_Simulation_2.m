    %%%%%% Stochastic Simulatin
%%%%%% Xiaoqiang Sun,Haofeng Lin @ SYSU
%%%%%% 2021,5
clear all
% close all
% clc

tic
 
%%% LHS_PRCC
clear all
% close all
% clc

tic


addpath('funfun'); 

global eta1    
global D_d1    
global delt_t  
global L  L2    
global B_3  
global A_0   
global d_max         
global period        
global inter         
global period1       
global inter1        
global TimeLength    
global m0        
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

Parameters_nondimensional_6_gai_jin_fangcheng_6
load('canshu1-24hours-ode15s.mat','p');
Ptest1=p(1,1:22);


TimeLength=400;
L2=100;
L=1/x_0 ;  %Total length of the axis during simulation
x = 0:L/L2:L;            %Spatial growth of the tumor
delt_t=1/91.5926;      %Time interval during simulation
t=0:delt_t:delt_t*TimeLength;

period=20;           %This is for boundary drug administration with a cycle
period1=period*delt_t;  %Actual half-cycle
inter=1;             %Interval
inter1=inter*delt_t;   %Actual Interval
eta1=eta_d;
D_d1=D_d;


dmax1=1;                  
dmax2=1;
dmax3=1;

select_CSF1R_I=1;
select_EGFR_I=0;
select_IGF1R_I=0;

global A_T_EGF A_T_IGF1 r_T a_ERK a_AKT K1 K2 K310 K320 K4 A_M1 epsilon d_T d_TM1 a_M2M1 a_CSF1 a_A B1 B2 B3 D_d1 eta1
Parameter=[ A_T_EGF A_T_IGF1 r_T a_ERK a_AKT K1 K2 K310 K320 K4 A_M1 epsilon d_T d_TM1 a_M2M1 a_CSF1 a_A B1 B2 B3 D_d1 eta1];

p = 100;   % Number of points
N = 1;   % Number of dimensions

mu=[1.05*a_ERK 0.05*a_A];
sigma=[0.005*a_ERK 0;0 0.03*a_A];
% sigma=[1 0 0;0 1 0;0 0 1];
for i=1:1:1000
bm41=mvnrnd(mu,sigma,p);
aaaaa=(bm41<=0);
bbbbb=sum(aaaaa,1);
ccccc1=sum(bbbbb,2);
if ccccc1==0
    break;
end
i=i+1;
end

mu=[1.05*a_ERK 2.47*a_A];
sigma=[0.005*a_ERK 0;0 32*a_A];
% sigma=[1 0 0;0 1 0;0 0 1];
for i=1:1:1000
bm42=mvnrnd(mu,sigma,p);
aaaaa=(bm42<=0);
bbbbb=sum(aaaaa,1);
ccccc2=sum(bbbbb,2);
if ccccc2==0
    break;
end
i=i+1;
end

bm4=zeros(p,2);
bm4=bm41;
rrr=0.6;
for i=1:1:p
rrr1=rand;

bm4(i,:)=bm41(i,:)*heaviside(rrr1-rrr)+bm42(i,:)*heaviside(rrr-rrr1);
end

aaaaaaaaaa=0;
figure
plot(1:p,bm4(:,1));
figure
plot(1:p,bm4(:,2));

r_RG=zeros(1,p);
r_K=zeros(1,p);

C=zeros(p,TimeLength+1);



t6=zeros(100,500);
global iii
for i=1:p
    iii=i
    a_ERK=bm4(i,1);
    a_A=bm4(i,2);


options=odeset('reltol',1e-3);%accuracy of solution
m = 2;%dimensionality of the solution
sol = pdepe2(m,@pdex4pde,@odex4ode,@pdex4ic,@pdex4bc,x,t,options);
u1 = sol(:,:,1);
% u2 = sol(:,:,2);
% u3 = sol(:,:,3);
% u4 = sol(:,:,4);
% u5 = sol(:,:,5);
% u6 = sol(:,:,6);
% u7 = sol(:,:,7);
% u8 = sol(:,:,8);
% u9 = sol(:,:,9);
% u10 = sol(:,:,10);
t6(i,:)=t5/delt_t;
C_T_mean=sum(u1(:,:),2)/100; 
C(i,:)=C_T_mean;

    if C(i,TimeLength+1)~=min(C(i,:))
        r_RG(i)=(C(i,TimeLength+1)-min(C(i,:)))/(TimeLength+1-find(C(i,:)==min(C(i,:)),1));
    else
        r_RG(i)=0;
    end
    if C(i,1)~=min(C(i,:))
        r_K(i)=(C(i,1)-min(C(i,:)))/(find(C(i,:)==min(C(i,:)),1));
    else
        r_K(i)=0;
    end
cc=(C<0)
end

save('lianxuyongdrug_CSF1RI_1_400.mat','C','r_RG','r_K','bm4','t6');  

threshold=0.697;
[p p12]=size(C);
C0=C;
C1=zeros(p,1);
C2=C;
C3=zeros(p,1);

for person=1:1:p
    if isempty(find(C2(person,:)>=threshold,1))
        aaaaaa=max(C2(person,:));
        aaaaaa1=find(C2(person,:)==aaaaaa);
        if aaaaaa1~=p12
            if aaaaaa1+50<=p12
                if aaaaaa>=0.67
                    if (aaaaaa-sum(C2(person,aaaaaa1+1:aaaaaa1+50))/50)/aaaaaa<0.005
                        C3(person,1)=aaaaaa1;
                    end                        
                end
            end
        end
    else
        C3(person,1)=find(C2(person,:)>=threshold,1); 
    end
end
C4=zeros(TimeLength+1,1);
for i=1:1:TimeLength+1  
    if isempty(find(C3==i,1))
        aaaaa=i;
    else
       aaa=find(C3==i);
       C4(i,1)=length(aaa);        
    end
end
Survival_Rate=100*ones(TimeLength+1,1);
qqq=p;
for t1=2:1:TimeLength+1    
    qqq=qqq-C4(t1,1);
    Survival_Rate(t1,1)=qqq*100/p;
end
C5=zeros(TimeLength/20+1,1);
for i=1:1:TimeLength/20
    C5(i+1,1)=sum(C4((i-1)*20+1:i*20,1));
end

ST1=[0 66 34 0 0 0 0 0];%data
PST1=C5(1:8,1);
figure,plot([0:20:140],ST1,'b o-- ','Linewidth',2.5,'MarkeredgeColor','b','MarkerfaceColor','b');axis([0,140,0,100])
hold on, plot([0:20:140],PST1,'r o-- ','Linewidth',2.5,'MarkeredgeColor','r','MarkerfaceColor','r');axis([0,140,0,100])
set(gca,'xtick',(0:20:140));
xlim([0 140]);
legend('Experimental data','Prediction')
title('No treatment','FontWeight','Bold','FontSize',18);
set(gca,'xticklabel',{'0','20','40','60','80','100','120','140'})
xlabel('Time (Days)','FontWeight','Bold','FontSize',18);
ylabel('Frequency','FontWeight','Bold','FontSize',18);
set(gca,'FontSize',18);     
figure_FontSize=18;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',18),'FontSize',figure_FontSize);

ST2=[0 0 0 9 25 21 1 0];%data
PST2=C5(1:8,1);
figure,plot(0:20:140,ST2,'o--','Color',[64/255,116/255,52/255],'Linewidth',2.5,'MarkeredgeColor',[64/255,116/255,52/255],'MarkerfaceColor',[64/255,116/255,52/255])
hold on, plot(0:20:140,PST2,'o-- ','Color','r','Linewidth',2.5,'MarkeredgeColor','r','MarkerfaceColor','r')
xlim([0 140]);
ylim([0 100]);
legend('Experimental data','Prediction')
title('CSF1R Inhibitor','FontWeight','Bold','FontSize',18);
set(gca,'xtick',0:20:140);
set(gca,'xticklabel',{'0','20','40','60','80','100','120','140'})
xlabel('Survival time (days)','FontWeight','Bold','FontSize',18);
ylabel('Frequency','FontWeight','Bold','FontSize',18);
set(gca,'FontSize',18);     
figure_FontSize=18;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',18),'FontSize',figure_FontSize);


%
load('canshucanshu2yue24xin4_nodrug_CSFRI_1_400day_lianxu_111.mat','Survival_Rate','C5','C3','C','bm4');  
ST1=[0 66 34 0 0 0 0 0];%data
PST1=C5(1:8,1);
figure,plot([0:20:140],ST1,'b o-- ','Linewidth',2.5,'MarkeredgeColor','b','MarkerfaceColor','b');axis([0,140,0,100])
hold on, plot([0:20:140],PST1,'r o-- ','Linewidth',2.5,'MarkeredgeColor','r','MarkerfaceColor','r');axis([0,140,0,100])
set(gca,'xtick',(0:20:140));
xlim([0 140]);
legend('Experimental data','Prediction')
title('No treatment','FontWeight','Bold','FontSize',18);
set(gca,'xticklabel',{'0','20','40','60','80','100','120','140'})
xlabel('Time (Days)','FontWeight','Bold','FontSize',18);
ylabel('Frequency','FontWeight','Bold','FontSize',18);
set(gca,'FontSize',18);     
figure_FontSize=18;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',18),'FontSize',figure_FontSize);
% 
load('canshucanshu2yue24xin4_yongdrug_CSFRI_1_400day_lianxu_111.mat','Survival_Rate','C5','C3','C','bm4');  
ST2=[0 0 0 9 25 21 1 0];%data
PST2=C5(1:8,1);
figure,plot(0:20:140,ST2,'o--','Color',[64/255,116/255,52/255],'Linewidth',2.5,'MarkeredgeColor',[64/255,116/255,52/255],'MarkerfaceColor',[64/255,116/255,52/255])
hold on, plot(0:20:140,PST2,'o-- ','Color','r','Linewidth',2.5,'MarkeredgeColor','r','MarkerfaceColor','r')
xlim([0 140]);
ylim([0 100]);
legend('Experimental data','Prediction')
title('CSF1R Inhibitor','FontWeight','Bold','FontSize',18);
set(gca,'xtick',0:20:140);
set(gca,'xticklabel',{'0','20','40','60','80','100','120','140'})
xlabel('Survival time (days)','FontWeight','Bold','FontSize',18);
ylabel('Frequency','FontWeight','Bold','FontSize',18);
set(gca,'FontSize',18);     
figure_FontSize=18;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',18),'FontSize',figure_FontSize);

figure,
for person=1:p
    hold on;
plot(0:140,C(person,1:140+1));
end
set(gca,'xtick',0:20:140); 
set(gca,'xticklabel',0:20:140);
xlim([0 140]);
% title('d_T_M_1 ')
xlabel('Time (days)','FontWeight','Bold','FontSize',18);
ylabel('Cancer cell density','FontWeight','Bold','FontSize',18);
set(gca,'FontSize',18);     
figure_FontSize=18;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',18),'FontSize',figure_FontSize);
title('Stochastic simulation','FontWeight','Bold','FontSize',18);


load('canshucanshu2yue24xin4_nodrug_CSFRI_1_400day_lianxu_111.mat','Survival_Rate','C5','C3','C','bm4'); 
Time_1=140;
figure,
plot(0:1:Time_1,Survival_Rate(1:Time_1+1,1),'color',[150/255,150/255,150/255],'Linewidth',2.5,'MarkeredgeColor','b','MarkerfaceColor','y');
set(gca,'xtick',0:20:Time_1); 
set(gca,'xticklabel',0:20:Time_1);
% xlim([0 TimeLength+1]);
xlim([0 Time_1+1]);
ylim([0 100]);
xlabel('Time(days)','FontWeight','Bold','FontSize',18);
ylabel('Percent survival','FontWeight','Bold','FontSize',18);
set(gca,'FontSize',18);     
figure_FontSize=18;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',18),'FontSize',figure_FontSize);
hold on

load('canshucanshu2yue24xin4_yongdrug_CSFRI_1_400day_lianxu_111.mat','Survival_Rate','C5','C3','C','bm4');  
plot(0:1:Time_1,Survival_Rate(1:Time_1+1,1),'color',[220/255,20/255,60/255],'Linewidth',2.5,'MarkeredgeColor','b','MarkerfaceColor','y');
set(gca,'xtick',0:20:Time_1); 
set(gca,'xticklabel',0:20:Time_1);
% xlim([0 TimeLength+1]);
xlim([0 Time_1+1]);
ylim([0 100]);
xlabel('Time(days)','FontWeight','Bold','FontSize',18);
ylabel('Percent survival','FontWeight','Bold','FontSize',18);
set(gca,'FontSize',18);     
figure_FontSize=18;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',18),'FontSize',figure_FontSize);
legend('No treatment','CSF1R Inhibitor')








load('canshucanshu2yue24xin4_yongdrug_CSFRI_1_400day_5yue8special_CSF1RI_lianxu_EGFRI_400_111.mat','Survival_Rate','C5');  

Time_1=400;
figure,
plot(0:1:Time_1,Survival_Rate(1:Time_1+1,1),'color',[77/255 175/255 74/255],'Linewidth',2.5,'MarkeredgeColor','b','MarkerfaceColor','y');
set(gca,'xtick',0:Time_1/5:Time_1); 
set(gca,'xticklabel',0:Time_1/5:Time_1);
% xlim([0 TimeLength+1]);
xlim([0 Time_1+1]);
ylim([0 100]);
xlabel('Time(days)','FontWeight','Bold','FontSize',18);
ylabel('Percent Survival','FontWeight','Bold','FontSize',18);
set(gca,'FontSize',18);     
figure_FontSize=18;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',18),'FontSize',figure_FontSize);
hold on

 
load('canshucanshu2yue24xin4_yongdrug_CSFRI_1_400day_5yue8special_CSFRI_1_EGFRI_1_400day_zishiying_total_111.mat','Survival_Rate','C5');  
plot(0:1:Time_1,Survival_Rate(1:Time_1+1,1),'color',[231/255 41/255 138/255],'Linewidth',2.5,'MarkeredgeColor','b','MarkerfaceColor','y');
set(gca,'xtick',0:Time_1/5:Time_1); 
set(gca,'xticklabel',0:Time_1/5:Time_1);
% xlim([0 TimeLength+1]);
xlim([0 Time_1+1]);
ylim([0 100]);
xlabel('Time(days)','FontWeight','Bold','FontSize',18);
ylabel('Percent Survival','FontWeight','Bold','FontSize',18);
set(gca,'FontSize',18);     
figure_FontSize=18;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',18),'FontSize',figure_FontSize);
hold on


load('canshucanshu2yue24xin4_yongdrug_CSFRI_1_400day_5yue8special_CSF1RI_lianxu_IGF1RI_400_111.mat','Survival_Rate','C5');  
plot(0:1:Time_1,Survival_Rate(1:Time_1+1,1),'Color',[152/255 78/255 163/255],'Linewidth',2.5,'MarkeredgeColor','b','MarkerfaceColor','y');
set(gca,'xtick',0:Time_1/5:Time_1); 
set(gca,'xticklabel',0:Time_1/5:Time_1);
% xlim([0 TimeLength+1]);
xlim([0 Time_1+1]);
ylim([0 100]);
xlabel('Time(days)','FontWeight','Bold','FontSize',18);
ylabel('Percent Survival','FontWeight','Bold','FontSize',18);
set(gca,'FontSize',18);     
figure_FontSize=18;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',18),'FontSize',figure_FontSize);
hold on
% % 

load('canshucanshu2yue24xin4_yongdrug_CSFRI_1_400day_5yue8special_CSFRI_1_IGF1RI_1_400day_zishiying_total_111.mat','Survival_Rate','C5');  
plot(0:1:Time_1,Survival_Rate(1:Time_1+1,1),'Color',[255/255 127/255 0/255],'Linewidth',2.5,'MarkeredgeColor','b','MarkerfaceColor','y');
set(gca,'xtick',0:Time_1/5:Time_1); 
set(gca,'xticklabel',0:Time_1/5:Time_1);
% xlim([0 TimeLength+1]);
xlim([0 Time_1+1]);
ylim([0 100]);
xlabel('Time(days)','FontWeight','Bold','FontSize',18);
ylabel('Percent Survival','FontWeight','Bold','FontSize',18);
set(gca,'FontSize',18);     
figure_FontSize=18;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',18),'FontSize',figure_FontSize);


% legend('No Drug','Continuous  CSF1RI','Continuous  CSF1RI & Continuous EGFRI', 'Continuous  CSF1RI & Continuous IGF1RI')
    legend('Continuous  CSF1RI & Continuous EGFRI','Adaptive  CSF1RI & Continuous EGFRI', ...
    'Continuous  CSF1RI & Continuous IGF1RI',...
    'Adaptive  CSF1RI & Continuous IGF1RI' );


