%%% LHS_PRCC
%%%%%% Xiaoqiang Sun,Haofeng Lin @ SYSU
%%%%%% 2021,5
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


TimeLength=200;
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

select_CSF1R_I=0;
select_EGFR_I=0;
select_IGF1R_I=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters
% global D_IGF1 d_IGF1 D_T r1 K0 Alpha Kp d1 Beta Kd
% global A1 A2 R1 D1 B1 B2 e D_d eta  r1 d1
% global C_T_mean
global A_T_EGF A_T_IGF1 r_T a_ERK a_AKT K1 K2 K31 K32 K4 A_M1 epsilon d_T d_TM1 a_M2M1 a_CSF1 a_A B1 B2 B3 D_d1 eta1
Parameter=[ A_T_EGF A_T_IGF1 r_T a_ERK a_AKT K1 K2 K31 K32 K4 A_M1 epsilon d_T d_TM1 a_M2M1 a_CSF1 a_A B1 B2 B3 D_d1 eta1];

%%%%  LHS

p = 1000;   % Number of points
N = 22;   % Number of dimensions
% lb = zeros(1,9); % lower bounds for parameters
lb = [A_T_EGF A_T_IGF1 r_T a_ERK a_AKT K1 K2 K31 K32 K4 A_M1 epsilon d_T d_TM1 a_M2M1 a_CSF1 a_A B_1 B_2 B_3 D_d1 eta1]*0.9; % lower bounds for parameters
ub = [A_T_EGF A_T_IGF1 r_T a_ERK a_AKT K1 K2 K31 K32 K4 A_M1 epsilon d_T d_TM1 a_M2M1 a_CSF1 a_A B_1 B_2 B_3 D_d1 eta1]*1.1; % upper bounds for parameters
X = lhsdesign(p,N,'criterion','correlation');
D = bsxfun(@plus,lb,bsxfun(@times,X,(ub-lb)));

r_RG=zeros(1,p);
r_K=zeros(1,p);

C=zeros(p,TimeLength+1);
%%%%%%%  Model Simulation
for i=1:1:p
    i
    A_T_EGF=D(i,1);
    A_T_IGF1=D(i,2);
    r_T=D(i,3);
    a_ERK=D(i,4);
    a_AKT=D(i,5);
    K1=D(i,6);
    K2=D(i,7);
    K31=D(i,8);
    K32=D(i,9);
    K4=D(i,10);
    A_M1=D(i,11);
    epsilon=D(i,12);
    d_T=D(i,13);
    d_TM1=D(i,14);
    a_M2M1=D(i,15);
    a_CSF1=D(i,16);
    a_A=D(i,17);
    B_1=D(i,18);
    B_2=D(i,19);
    B_3=D(i,20);
    D_d1=D(i,21);
    eta1=D(i,22);
    
    
options=odeset('reltol',1e-3);
m = 2;
sol = pdepe2(m,@pdex4pde,@odex4ode,@pdex4ic,@pdex4bc,x,t,options);
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
C_T_mean=sum(u1(:,:),2)/(L2+1); 

C(i,:)=C_T_mean;

[asizeC,bsizeC]=size(C);
TimeLength=bsizeC-1;
threshold=0.697;



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

[asizeC,bsizeC]=size(C);
TimeLength=bsizeC-1;
threshold=0.697;
[p p12]=size(C);
C0=C;
C1=zeros(p,1);
C2=C;
C3=zeros(p,1);
%C3 represents the number of days since death, C4 represents the number of deaths per day, 
% and C5 represents the number of deaths every 20 days.
for person=1:1:p
    if isempty(find(C2(person,:)>=threshold,1))
        aaaaaa=max(C2(person,:));
        aaaaaa1=find(C2(person,:)==aaaaaa);
        if aaaaaa1~=p12
            if aaaaaa1+50<=p12
                if aaaaaa>=0.6
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

save('1000personcanshuminganxingfenxi_nodrug_C3_0.9-1.1.mat','C','r_RG','r_K','D','C3','C4','C5');   

load('1000person_nodrug_C3_11para_0.9_1.1_6.mat','C','r_RG','r_K','D','C3','C4','C5');
load('1000person_CSF1RI_1_C3_11para_0.9_1.1_6.mat','C','r_RG','r_K','D','C3','C4','C5');   

[aa1,aa2]=size(D);
N=aa2;
if aa2==22
%%% Ranking
    A_T_EGF_1=D(:,1);
    A_T_IGF1_1=D(:,2);
    r_T_1=D(:,3);
    a_ERK_1=D(:,4);
    a_AKT_1=D(:,5);
    K1_1=D(:,6);
    K2_1=D(:,7);
    K31_1=D(:,8);
    K32_1=D(:,9);
    K4_1=D(:,10);
    A_M1_1=D(:,11);
    epsilon_1=D(:,12);
    d_T_1=D(:,13);
    d_TM1_1=D(:,14);
    a_M2M1_1=D(:,15);
    a_CSF1_1=D(:,16);
    a_A_1=D(:,17);
    B_1_1=D(:,18);
    B_2_1=D(:,19);
    B_3_1=D(:,20);
    D_d1_1=D(:,21);
    eta1_1=D(:,22);
    P_O=[A_T_EGF_1,A_T_IGF1_1,r_T_1,a_ERK_1,a_AKT_1,K1_1,K2_1,K31_1,K32_1,K4_1,A_M1_1,epsilon_1,d_T_1,d_TM1_1,a_M2M1_1,a_CSF1_1,a_A_1,B_1_1,B_2_1,B_3_1,D_d1_1,eta1_1];
    Para_label={'A_{TEGF}','A_{TIGF1}','r_T','{\alpha}_{ERK}','{\alpha}_{AKT}','K1','K2','K31','K32','K4','A_{M1}','{\epsilon}','d_T','d_{TM1}','{\alpha}_{M2M1}','{\alpha}_{CSF1}','{\alpha}_A','B_1','B_2','B_3','D_d','{\eta}'}; 
elseif aa2==11
    a_ERK_1=D(:,1);
    a_AKT_1=D(:,2);
    K1_1=D(:,3);
    K2_1=D(:,4);
    K31_1=D(:,5);
    K32_1=D(:,6);
    K4_1=D(:,7);
    d_TM1_1=D(:,8);
    a_M2M1_1=D(:,9);
    a_CSF1_1=D(:,10);
    a_A_1=D(:,11);
    P_O=[a_ERK_1,a_AKT_1,K1_1,K2_1,K31_1,K32_1,K4_1,d_TM1_1,a_M2M1_1,a_CSF1_1,a_A_1];
    Para_label={'{\alpha}_{ERK}','{\alpha}_{AKT}','K1','K2','K31','K32','K4','d_{TM1}','{\alpha}_{M2M1}','{\alpha}_{CSF1}','{\alpha}_A'};  
else 
    error('Cannot calculate with given values')
end

    r_RG_O=r_RG';
    r_K_O=r_K';

    C3_O=C3;

    [P_R,ind_P]=sort(P_O);
    [r_RG_R,ind_RG]=sort(r_RG_O);
    [r_K_R,ind_K]=sort(r_K_O);

    [P_R,ind_P]=sort(P_O);
    [r_RG_R,ind_RG]=sort(C3_O);

        CC=zeros(N,1);
        Pvalue=zeros(N,1);
%%%%  Pearson CC for raw data
       for j=1:N
         a=1:N;
         a(a==j)=[];
%         [rho,pval] = corr(r_RG_O,P_O(:,j));
        [rho,pval] = corr(C3_O,P_O(:,j));
        CC(j)=rho;
        Pvalue(j)=pval;
       end

b1=[110/255,165/255,211/255];b2=[241/255,153/255,157/255];
figure,
bar(1:N,CC)  %
set(gca,'xtick',1:1:N)
set(gca,'xticklabel',Para_label)
rotateticklabel(gca,315);
% xlabel('Parameters','FontWeight','Bold','FontSize',18);
ylabel('Sensitivity Coefficient','FontWeight','Bold','FontSize',18);
title('Correlation Coefficient','FontWeight','Bold','FontSize',18);

set(gca,'FontSize',14);     
figure_FontSize=18;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',14),'FontSize',figure_FontSize);

figure,
bar(1:N,Pvalue)  %
set(gca,'xtick',1:1:N)
set(gca,'xticklabel',Para_label)
rotateticklabel(gca,315);
% xlabel('Parameters','FontWeight','Bold','FontSize',18);
ylabel('Sensitivity Coefficient','FontWeight','Bold','FontSize',18);
title('Correlation Coefficient Pvalue','FontWeight','Bold','FontSize',18);

set(gca,'FontSize',14);     
figure_FontSize=18;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',14),'FontSize',figure_FontSize);


    toc
    
    









