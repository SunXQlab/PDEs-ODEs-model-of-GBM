function [c,f1,s1] = pdex4pde(x,t,u,DuDx)
global A_T_EGF A_T_IGF1 A_M1 A_M2 epsilon r_T a_ERK K1 K2 K31 K32 K4 a_AKT  d_T d_TM1 a_M2M1 a_CSF1 a_A B_1 B_2 b_EGF
global CSF1_max EGF_max IGF1_max
global ERK AKT 
global x_0
global Ptest1
global delt_t
global select_CSF1R_I dmax1

c = [1; 1; 1; 1; 1; 1];  
t%write time
f1 = [1,0,0,0,-A_T_EGF*u(1),-A_T_IGF1*u(1);
      0,1,0,-A_M1*u(2),0,0;
      0,0,1,-A_M2*u(3),0,0;
      0,0,0,1/epsilon,0,0;
      0,0,0,0,1/epsilon,0;
      0,0,0,0,0,1/epsilon]*DuDx; 

     drug_CSF1R_I=Drug2(x,t,select_CSF1R_I,dmax1);
     A=A_Drugsimulation2(x,t,select_CSF1R_I,dmax1);
     R1=r_T*(1+a_ERK*u(8)/(K1+u(8))+a_AKT*u(10)/(K2+u(10)));
     D1=d_T*(1+d_TM1*u(2));  % death rate    
     R3=a_M2M1*u(3)-u(2)*(a_CSF1*u(4)/(K31+K32*drug_CSF1R_I+u(4))+a_A*A/(K4+A));

s1 = [R1*u(1)*(1-u(1))-D1*u(1); R3; -R3;B_1/epsilon*(u(1)*(1-u(4)/CSF1_max)-u(4));1/epsilon*B_2*(u(3)*(1-u(5)/EGF_max)-u(5));1/epsilon*(u(3)*A*(1-u(6)/IGF1_max)-u(6))]; 

end