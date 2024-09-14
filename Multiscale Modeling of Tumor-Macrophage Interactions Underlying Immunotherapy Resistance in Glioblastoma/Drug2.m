function aaa = Drug2(x,t,select,dmax,dmax_1)

global d_max m0 b1 b2
% if isempty(dmax)
% else
    d_max=dmax;
    m0=dmax;
%     b1=dmax;
% end
% if nargin >=5
%     b2=dmax_1;
% end
if x==0
    x=0.0000001;
end
switch select
    case 0
        aaa=0;
    case 1
        aaa=Drug_1(x,t);
    case 2
        aaa=Drug_2(x,t);
    case 3
        aaa=Drug_3(x,t);
    case 4
        aaa=Drug_4(x,t);
    case 5
        aaa=Drug_5(x,t) ;
end
end

%Constant function
function q = Drug_1(x,t)
global d_max
global eta1
global D_d1
global L   %Maximum value of x
   C_d=d_max;
   drugt1=0;   
for n=1:1:100
    fi_n=(2*C_d*L/n/pi*(-1)^(n+1)+4*C_d*L/(n*pi)^3*((-1)^n-1)+C_d*(2*L)/n/pi*(-1)^n);
    Tn1=(fi_n*exp(-(D_d1*(n*pi/L)^2+eta1)*t)+2*eta1*C_d*L/(n*pi)*(-1)^n/(D_d1*(n*pi/L)^2+eta1)*(1-exp(-(D_d1*(n*pi/L)^2+eta1)*t)));
    Tn11=Tn1*sin(n*pi*x/L);
    drugt1=drugt1+Tn11; 
end
q=1/x*drugt1+C_d;
end

%Sine function
function q = Drug_2(x,t)
global d_max
global eta1
global D_d1
global L   %Maximum value of x
global  period_2
% d_initial=d_max/L;
w=2*pi/period_2;
% w=pi/(TimeLength*delt_t/(period1/delt_t));
   C_d=d_max;
   C_d=C_d*2/3;
   d22=C_d/2*(sin(w*t)+2);
drugt2=0;
for n=1:1:100
    fi_n=(2*C_d*L/n/pi*(-1)^(n+1)+4*C_d*L/(n*pi)^3*((-1)^n-1)+C_d*(2*L)/n/pi*(-1)^n);
    Tn2_1=fi_n*exp(-(D_d1*(n*pi/L)^2+eta1)*t);
    Tn2_2=2*eta1*C_d*L/(n*pi)*(-1)^n/(D_d1*(n*pi/L)^2+eta1)*(1-exp(-(D_d1*(n*pi/L)^2+eta1)*t));
    Tn2_3=eta1*C_d*(D_d1*(n*pi/L)^2+eta1)*L/(n*pi)*(-1)^n/((D_d1*(n*pi/L)^2+eta1)^2+w^2)*sin(w*t);
    Tn2_4=-eta1*C_d*w*L/(n*pi)*(-1)^n/((D_d1*(n*pi/L)^2+eta1)^2+w^2)*cos(w*t);
    Tn2_5=eta1*C_d*w*L/(n*pi)*(-1)^n/((D_d1*(n*pi/L)^2+eta1)^2+w^2)*exp(-(D_d1*(n*pi/L)^2+eta1)*t);
    Tn2_6=C_d*w*(D_d1*(n*pi/L)^2+eta1)*L/(n*pi)*(-1)^n/((D_d1*(n*pi/L)^2+eta1)^2+w^2)*cos(w*t);
    Tn2_7=-C_d*w*(D_d1*(n*pi/L)^2+eta1)*L/(n*pi)*(-1)^n/((D_d1*(n*pi/L)^2+eta1)^2+w^2)*exp(-(D_d1*(n*pi/L)^2+eta1)*t);
    Tn2_8=C_d*w^2*L/(n*pi)*(-1)^n/((D_d1*(n*pi/L)^2+eta1)^2+w^2)*sin(w*t);
    Tn2=Tn2_1+Tn2_2+Tn2_3+Tn2_4+Tn2_5+Tn2_6+Tn2_7+Tn2_8;
    drugt2=drugt2+Tn2*sin(n*pi*x/L);
end
q=1/x*drugt2+d22;
end


%%This involves the periodic medication function, 
% with Hermite interpolation used to connect the functions for calculating Tn, where integration is involved, it's done through 
% cumulative summation of the integral.
%The first part
%The integral function one, corresponding to h1(t, t1), where t1 is equal to to.
function q=h_1_ceshiair2_1(k,t1,t,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
q1=m0/(D_d1*(n*pi/L)^2+eta1)*(exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t1))-exp(-(D_d1*(n*pi/L)^2+eta1)*(t-k*T)));
q=2*eta1*L*(-1)^n/n/pi*q1;
end
%The second integral function, corresponding to h2(t, t1).
function q=h_2_ceshiair2_2(k,t1,t,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q11=m0/d^2*(1+2/d*(t1-(k+1/2)*T+d))*(t1-(k+1/2)*T)^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t1))-m0*exp(-(D_d1*(n*pi/L)^2+eta1)*(t-(k+1/2)*T+d));
q1=q11/(D_d1*(n*pi/L)^2+eta1);
q21=(2*m0/d^3*(t1-(k+1/2)*T)^2+2*m0/d^2*(t1-(k+1/2)*T)*(1+2/d*(t1-(k+1/2)*T+d)))*exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t1));
q2=-q21/(D_d1*(n*pi/L)^2+eta1)^2;
q31=(8*m0/d^3*(t1-(k+1/2)*T)+2*m0/d^2*(1+2/d*(t1-(k+1/2)*T+d)))*exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t1))+6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(t-(k+1/2)*T+d));
q3=q31/(D_d1*(n*pi/L)^2+eta1)^3;
q41=12*m0/d^3*(exp(-(D_d1*(n*pi/L)^2+eta1)*(t-(k+1/2)*T+d))-exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t1)));
q4=q41/(D_d1*(n*pi/L)^2+eta1)^4;
q=(q1+q2+q3+q4)*2*eta1*L*(-1)^n/n/pi;
end
%The second integral function, corresponding to h3(t, t1).
function q=h_3_ceshiair2_3(k,t1,t,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q11=m0/d^2*(1-2/d*(t1-(k+1)*T))*(t1-(k+1)*T+d)^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t1));
q1=q11/(D_d1*(n*pi/L)^2+eta1);
q21=-(-2*m0/d^3*(t1-(k+1)*T+d)^2+2*m0/d^2*(t1-(k+1)*T+d)*(1-2/d*(t1-(k+1)*T)))*exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t1));
q2=q21/(D_d1*(n*pi/L)^2+eta1)^2;
q31=(-8*m0/d^3*(t1-(k+1)*T+d)+2*m0/d^2*(1-2/d*(t1-(k+1)*T)))*exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t1))-6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(t-(k+1)*T+d));
q3=q31/(D_d1*(n*pi/L)^2+eta1)^3;
q41=12*m0/d^3*(exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t1))-exp(-(D_d1*(n*pi/L)^2+eta1)*(t-(k+1)*T+d)));
q4=q41/(D_d1*(n*pi/L)^2+eta1)^4;
q=(q1+q2+q3+q4)*2*eta1*L*(-1)^n/n/pi;
end
%The second integral function, corresponding to h4(t, t1).
function q=h_4_ceshiair2_4(k,t1,t,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q11=(2*m0/d^3*(t1-(k+1/2)*T)^2+2*m0/d^2*(1+2/d*(t1-(k+1/2)*T+d))*(t1-(k+1/2)*T))*exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t1));
q1=q11/(D_d1*(n*pi/L)^2+eta1);
q21=-((8*m0/d^3*(t1-(k+1/2)*T)+2*m0/d^2*(1+2/d*(t1-(k+1/2)*T+d)))*exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t1))+6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1))*(t-(k+1/2)*T+d));
q2=q21/(D_d1*(n*pi/L)^2+eta1)^2;
q31=12*m0/d^3*(exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t1))-exp(-(D_d1*(n*pi/L)^2+eta1)*(t-(k+1/2)*T+d)));
q3=q31/(D_d1*(n*pi/L)^2+eta1)^3;
q=(q1+q2+q3)*2*L*(-1)^n/n/pi;
end
%The second integral function, corresponding to h5(t, t1).
function q=h_5_ceshiair2_5(k,t1,t,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q11=(-2*m0/d^3*(t1-(k+1)*T+d)^2+2*m0/d^2*(1-2/d*(t1-(k+1)*T))*(t1-(k+1)*T+d))*exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t1));
q1=q11/(D_d1*(n*pi/L)^2+eta1);
q21=(8*m0/d^3*(t1-(k+1)*T+d)-2*m0/d^2*(1-2/d*(t1-(k+1)*T)))*exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t1))+6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1))*(t-(k+1)*T+d);
q2=q21/(D_d1*(n*pi/L)^2+eta1)^2;
q31=-12*m0/d^3*(exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t1))-exp(-(D_d1*(n*pi/L)^2+eta1)*(t-(k+1)*T+d)));
q3=q31/(D_d1*(n*pi/L)^2+eta1)^3;
q=(q1+q2+q3)*2*L*(-1)^n/n/pi;
end

%The periodic medication function.
function q = Tn34(t,n)
global D_d1;
global eta1;
global L
global period1
global inter1
global m0
T=2*period1;
d=inter1;
q=0;
q1=0;
q2=0;
q3=0;
q4=0;
d_0=m0;
d_initial=d_0/L;
d1=floor(t/period1);
d2=floor(d1/2);
d3=rem(d1,2);
t0=t-d1*period1;
if d2~=0 
for i=0:1:(d2-1)
    t1=i*T;
    q1=h_1_ceshiair2_1(i,t1,t,n);
    t1=(i+1/2)*T-d;
    q2=h_1_ceshiair2_1(i,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
    t1=(i+1/2)*T-d;
    q1=h_2_ceshiair2_2(i,t1,t,n);
    t1=(i+1/2)*T;
    q2=h_2_ceshiair2_2(i,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
    t1=(i+1)*T-d;
    q1=h_3_ceshiair2_3(i,t1,t,n);
    t1=(i+1)*T;
    q2=h_3_ceshiair2_3(i,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
    t1=(i+1/2)*T-d;
    q1=h_4_ceshiair2_4(i,t1,t,n);
    t1=(i+1/2)*T;
    q2=h_4_ceshiair2_4(i,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
    t1=(i+1)*T-d;
    q1=h_5_ceshiair2_5(i,t1,t,n);
    t1=(i+1)*T;
    q2=h_5_ceshiair2_5(i,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
end
end
if d3==0 && t0<(period1-inter1)
    t1=d2*T;
    q1=h_1_ceshiair2_1(d2,t1,t,n);
    t1=t;
    q2=h_1_ceshiair2_1(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3; 
end
if d3==0 && t0>=(period1-inter1)
    t1=d2*T;
    q1=h_1_ceshiair2_1(d2,t1,t,n);
    t1=(d2+1/2)*T-d;
    q2=h_1_ceshiair2_1(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3; 
    t1=(d2+1/2)*T-d;
    q1=h_2_ceshiair2_2(d2,t1,t,n);
    t1=t;
    q2=h_2_ceshiair2_2(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3; 
    t1=(d2+1/2)*T-d;
    q1=h_4_ceshiair2_4(d2,t1,t,n);
    t1=t;
    q2=h_4_ceshiair2_4(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3; 
end
if d3==1 && t0<(period1-inter1)
    t1=d2*T;
    q1=h_1_ceshiair2_1(d2,t1,t,n);
    t1=(d2+1/2)*T-d;
    q2=h_1_ceshiair2_1(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3; 
    t1=(d2+1/2)*T-d;
    q1=h_2_ceshiair2_2(d2,t1,t,n);
    t1=(d2+1/2)*T;
    q2=h_2_ceshiair2_2(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3; 
    t1=(d2+1/2)*T-d;
    q1=h_4_ceshiair2_4(d2,t1,t,n);
    t1=(d2+1/2)*T;
    q2=h_4_ceshiair2_4(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3; 
end 
if d3==1 && t0>=(period1-inter1)
    t1=d2*T;
    q1=h_1_ceshiair2_1(d2,t1,t,n);
    t1=(d2+1/2)*T-d;
    q2=h_1_ceshiair2_1(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3; 
    t1=(d2+1/2)*T-d;
    q1=h_2_ceshiair2_2(d2,t1,t,n);
    t1=(d2+1/2)*T;
    q2=h_2_ceshiair2_2(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
    t1=(d2+1)*T-d;
    q1=h_3_ceshiair2_3(d2,t1,t,n);
    t1=t;
    q2=h_3_ceshiair2_3(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
    t1=(d2+1/2)*T-d;
    q1=h_4_ceshiair2_4(d2,t1,t,n);
    t1=(d2+1/2)*T;
    q2=h_4_ceshiair2_4(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
    t1=(d2+1)*T-d;
    q1=h_5_ceshiair2_5(d2,t1,t,n);
    t1=t;
    q2=h_5_ceshiair2_5(d2,t1,t,n);
    q3=q2-q1;
    q4=q4+q3;
end
fi_n=(2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n);
q=fi_n*exp(-(D_d1*(n*pi/L)^2+eta1)*t)+q4;
end


function q = Drug_3(x,t) 
global L
global period1
global inter1
global m0
q0=0;
d1=floor(t/period1);
d2=floor(d1/2);
d3=rem(d1,2);
t0=t-d1*period1;
xk=(d1+1)*period1-inter1;
xk1=(d1+1)*period1;
if d3==0 && t0<=(period1-inter1)
    dnt=m0;
end
if d3==0 && t0>(period1-inter1)
    yk=m0;
    %yk1=0;
    dnt=(1+2*(t-xk)/(xk1-xk))*((t-xk1)/(xk-xk1))^2*yk+0+0+0;   
end
if d3==1 && t0<=(period1-inter1)
    dnt=0;
end
if d3==1 && t0>(period1-inter1)
    %yk=0;
    yk1=m0;
    dnt=0+(1+2*(t-xk1)/(xk-xk1))*((t-xk)/(xk1-xk))^2*yk1+0+0;    
end
    for n=1:1:20
    q0=q0+Tn34(t,n)*sin(n*pi*x/L);  
    end
q=q0/x+dnt;
end


%The following is the function for medication cessation, situation 4 of drug administration.
function q=h_0_air(t,n)
global D_d1
global L
global inter1
global eta1
global b1
d=inter1;
q11=b1/(D_d1*(n*pi/L)^2+eta1);
q12=1-exp(-(D_d1*(n*pi/L)^2+eta1)*t);
q=q11*q12;
end

function q=h_1_air(t52,t53,t,n)
global D_d1
global L
global inter1
global eta1
global b1 b2
d=inter1;
q01=b1/(D_d1*(n*pi/L)^2+eta1);
q02=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))-exp(-(D_d1*(n*pi/L)^2+eta1)*t);
q0=q01*q02;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*b1/d^2*(t52-t53)^2;
q13=(b1/d^2*(1+2/d*(t-t52))*(t-t53)^2+b2/d^2*(1-2/d*(t-t53))*(t-t52)^2);
q1=q11*(q13-q12);
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(2*b1/d^3*(t52-t53)^2+2*b1/d^2*(t52-t53));
q23=(2*b1/d^3*(t-t53)^2+2*b1/d^2*(1+2/d*(t-t52))*(t-t53)-2*b2/d^3*(t-t52)^2+2*b2/d^2*(1-2/d*(t-t53))*(t-t52));
q2=q21*(q23-q22);
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(8*b1/d^3*(t52-t53)+2*b1/d^2+2*b2/d^2*(1-2/d*(t52-t53)));
q33=(8*b1/d^3*(t-t53)+2*b1/d^2*(1+2/d*(t-t52))-8*b2/d^3*(t-t52)+2*b2/d^2*(1-2/d*(t-t53)));
q3=q31*(q33-q32);
q41=-1/(D_d1*(n*pi/L)^2+eta1)^4;
q42=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(12*b1/d^3-12*b2/d^3);
q43=(12*b1/d^3-12*b2/d^3);
q4=q41*(q43-q42);
q=q0+q1+q2+q3+q4;
end

function q=h_2_air(t52,t53,t,n)
global D_d1
global L
global inter1
global eta1
global b1 b2
d=inter1;
q01=b1/(D_d1*(n*pi/L)^2+eta1);
q02=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))-exp(-(D_d1*(n*pi/L)^2+eta1)*t);
q0=q01*q02;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*b1/d^2*(t52-t53)^2;
q13=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t53))*(b2/d^2*(t53-t52)^2);
q1=q11*(q13-q12);
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(2*b1/d^3*(t52-t53)^2+2*b1/d^2*(t52-t53));
q23=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t53))*(-2*b2/d^3*(t53-t52)^2+2*b2/d^2*(t53-t52));
q2=q21*(q23-q22);
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(8*b1/d^3*(t52-t53)+2*b1/d^2+2*b2/d^2*(1-2/d*(t52-t53)));
q33=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t53))*(2*b1/d^2*(1+2/d*(t53-t52))-8*b2/d^3*(t53-t52)+2*b2/d^2);
q3=q31*(q33-q32);
q41=-1/(D_d1*(n*pi/L)^2+eta1)^4;
q42=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(12*b1/d^3-12*b2/d^3);
q43=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t53))*(12*b1/d^3-12*b2/d^3);
q4=q41*(q43-q42);
q51=b2/(D_d1*(n*pi/L)^2+eta1);
q52=1-exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t53));
q5=q51*q52;
q=q0+q1+q2+q3+q4+q5;
end


function q=h_3_air(t52,t53,t,n)
global D_d1
global L
global inter1
global period1
global eta1
global b1 b2
d=inter1;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(2*b1/d^3*(t52-t53)^2+2*b1/d^2*(t52-t53));
q13=(2*b1/d^3*(t-t53)^2+2*b1/d^2*(1+2/d*(t-t52))*(t-t53)-2*b2/d^3*(t-t52)^2+2*b2/d^2*(1-2/d*(t-t53))*(t-t52));
q1=q11*(q13-q12);
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(8*b1/d^3*(t52-t53)+2*b1/d^2+2*b2/d^2*(1-2/d*(t52-t53)));
q23=(8*b1/d^3*(t-t53)+2*b1/d^2*(1+2/d*(t-t52))-8*b2/d^3*(t-t52)+2*b2/d^2*(1-2/d*(t-t53)));
q2=q21*(q23-q22);
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(12*b1/d^3-12*b2/d^3);
q33=(12*b1/d^3-12*b2/d^3);
q3=q31*(q33-q32);
q=q1+q2+q3;
end

function q=h_4_air(t52,t53,t,n)
global D_d1
global L
global inter1
global period1
global eta1
global b1 b2
d=inter1;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(2*b1/d^3*(t52-t53)^2+2*b1/d^2*(t52-t53));
q13=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t53))*(-2*b2/d^3*(t53-t52)^2+2*b2/d^2*(t53-t52));
q1=q11*(q13-q12);
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(8*b1/d^3*(t52-t53)+2*b1/d^2+2*b2/d^2*(1-2/d*(t52-t53)));
q23=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t53))*(2*b1/d^2*(1+2/d*(t53-t52))-8*b2/d^3*(t53-t52)+2*b2/d^2);
q2=q21*(q23-q22);
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(12*b1/d^3-12*b2/d^3);
q33=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t53))*(12*b1/d^3-12*b2/d^3);
q3=q31*(q33-q32);
q=q1+q2+q3;
end

function q = Tn4(t,n)
global D_d1;
global eta1;
global L
global period1
global inter1                                        
global b1 b2
global t51 t52 t53
d_0=b1;
d_initial=d_0/L;
if t>=t51 && t<t52
    q0=h_0_air(t,n)*eta1*2*L/(n*pi)*(-1)^n;
    q1=0;
elseif t>=t52 && t<t53
    q0=h_1_air(t52,t53,t,n)*eta1*2*L/(n*pi)*(-1)^n;
    q1=h_3_air(t52,t53,t,n)*2*L/(n*pi)*(-1)^n;
elseif t>=t53
    q0=h_2_air(t52,t53,t,n)*eta1*2*L/(n*pi)*(-1)^n;
    q1=h_4_air(t52,t53,t,n)*2*L/(n*pi)*(-1)^n;
else
    msg='Error occured';
    error(msg)
end
q2=q0+q1;
q=(2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n)*exp(-(D_d1*(n*pi/L)^2+eta1)*t)+q2;
end

%The following is the function for medication cessation, situation 4 of drug administration.
function q = Drug_4(x,t) 
global L
global t51 t52 t53
global inter1
global b1 b2
d=inter1;
if t>=t51 && t<t52
    dnt=b1;
elseif t>=t52 && t<t53
    dnt=b1/d^2*(1+2/d*(t-t52))*(t-t53)^2+b2/d^2*(1-2/d*(t-t53))*(t-t52)^2;
elseif t>=t53
    dnt=b2;
else
        msg='Error occured';
    error(msg)
end
q0=0;
    for n=1:1:20
    q0=q0+Tn4(t,n)*sin(n*pi*x/L);  
    end
q=q0/x+dnt;
end



%Below is the adaptive medication function, situation 5 of drug administration.
function q=p1i_air(t51,t,n,bi)
global D_d1
global L
global inter1
global eta1
d=inter1;
q11=bi/(D_d1*(n*pi/L)^2+eta1);
q12=1-exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t51));
q=q11*q12;
end

function q=p2i_air(t52,t53,t,n,bi1,bi2)
global D_d1
global L
global inter1
global period1
global eta1
d=inter1;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*bi1/d^2*(t52-t53)^2;
q13=(bi1/d^2*(1+2/d*(t-t52))*(t-t53)^2+bi2/d^2*(1-2/d*(t-t53))*(t-t52)^2);
q1=q11*(q13-q12);
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(2*bi1/d^3*(t52-t53)^2+2*bi1/d^2*(t52-t53));
q23=(2*bi1/d^3*(t-t53)^2+2*bi1/d^2*(1+2/d*(t-t52))*(t-t53)-2*bi2/d^3*(t-t52)^2+2*bi2/d^2*(1-2/d*(t-t53))*(t-t52));
q2=q21*(q23-q22);
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(8*bi1/d^3*(t52-t53)+2*bi1/d^2+2*bi2/d^2*(1-2/d*(t52-t53)));
q33=(8*bi1/d^3*(t-t53)+2*bi1/d^2*(1+2/d*(t-t52))-8*bi2/d^3*(t-t52)+2*bi2/d^2*(1-2/d*(t-t53)));
q3=q31*(q33-q32);
q41=-1/(D_d1*(n*pi/L)^2+eta1)^4;
q42=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(12*bi1/d^3-12*bi2/d^3);
q43=(12*bi1/d^3-12*bi2/d^3);
q4=q41*(q43-q42);
q=q1+q2+q3+q4;
end
function q=p3i_air(t52,t53,t,n,bi1,bi2)
global D_d1
global L
global inter1
global period1
global eta1
d=inter1;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(2*bi1/d^3*(t52-t53)^2+2*bi1/d^2*(t52-t53));
q13=(2*bi1/d^3*(t-t53)^2+2*bi1/d^2*(1+2/d*(t-t52))*(t-t53)-2*bi2/d^3*(t-t52)^2+2*bi2/d^2*(1-2/d*(t-t53))*(t-t52));
q1=q11*(q13-q12);
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(8*bi1/d^3*(t52-t53)+2*bi1/d^2+2*bi2/d^2*(1-2/d*(t52-t53)));
q23=(8*bi1/d^3*(t-t53)+2*bi1/d^2*(1+2/d*(t-t52))-8*bi2/d^3*(t-t52)+2*bi2/d^2*(1-2/d*(t-t53)));
q2=q21*(q23-q22);
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(12*bi1/d^3-12*bi2/d^3);
q33=(12*bi1/d^3-12*bi2/d^3);
q3=q31*(q33-q32);
q=q1+q2+q3;
end
function q=q1i_air(t51,t52,t,n,bi)
global D_d1
global L
global inter1
global eta1
d=inter1;
q11=bi/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))-exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t51));
q=q11*q12;
end

function q=q2i_air(t52,t53,t,n,bi1,bi2)
global D_d1
global L
global inter1
global period1
global eta1
d=inter1;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*bi1/d^2*(t52-t53)^2;
q13=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t53))*(bi2/d^2*(t53-t52)^2);
q1=q11*(q13-q12);
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(2*bi1/d^3*(t52-t53)^2+2*bi1/d^2*(t52-t53));
q23=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t53))*(-2*bi2/d^3*(t53-t52)^2+2*bi2/d^2*(t53-t52));
q2=q21*(q23-q22);
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(8*bi1/d^3*(t52-t53)+2*bi1/d^2+2*bi2/d^2*(1-2/d*(t52-t53)));
q33=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t53))*(2*bi1/d^2*(1+2/d*(t53-t52))-8*bi2/d^3*(t53-t52)+2*bi2/d^2);
q3=q31*(q33-q32);
q41=-1/(D_d1*(n*pi/L)^2+eta1)^4;
q42=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(12*bi1/d^3-12*bi2/d^3);
q43=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t53))*(12*bi1/d^3-12*bi2/d^3);
q4=q41*(q43-q42);
q=q1+q2+q3+q4;
end
function q=q3i_air(t52,t53,t,n,bi1,bi2)
global D_d1
global L
global inter1
global period1
global eta1
d=inter1;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(2*bi1/d^3*(t52-t53)^2+2*bi1/d^2*(t52-t53));
q13=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t53))*(-2*bi2/d^3*(t53-t52)^2+2*bi2/d^2*(t53-t52));
q1=q11*(q13-q12);
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(8*bi1/d^3*(t52-t53)+2*bi1/d^2+2*bi2/d^2*(1-2/d*(t52-t53)));
q23=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t53))*(2*bi1/d^2*(1+2/d*(t53-t52))-8*bi2/d^3*(t53-t52)+2*bi2/d^2);
q2=q21*(q23-q22);
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t52))*(12*bi1/d^3-12*bi2/d^3);
q33=exp(-(D_d1*(n*pi/L)^2+eta1)*(t-t53))*(12*bi1/d^3-12*bi2/d^3);
q3=q31*(q33-q32);
q=q1+q2+q3;
end

function q = Tn5(t,n)
global D_d1;
global eta1;
global L
global period1
global inter1
global b1 b2 b3 b4
global t51 t52 t53 t54 t55 t56 t57 t58
global b t5
d_0=b(1);
d_initial=d_0/L;


aa1=find(t>=t5);
aa2=aa1(end);
aa3=floor(aa2/2);
aa4=aa2-aa3*2;
q11=0;q12=0;q13=0;
if aa3==0
        q1=eta1*p1i_air(t5(1),t,n,b(1));
elseif aa3==1
    if aa4==0
        q11=q11+eta1*(q1i_air(t5(1),t5(2),t,n,b(1)));
        q1=q11+eta1*p2i_air(t5(2),t5(3),t,n,b(1),b(2))+p3i_air(t5(2),t5(3),t,n,b(1),b(2));
    else
        q11=eta1*(q1i_air(t5(1),t5(2),t,n,b(1)));
        q12=eta1*q2i_air(t5(2),t5(3),t,n,b(1),b(2));
        q13=q3i_air(t5(2),t5(3),t,n,b(1),b(2));
        q1=q11+q12+q13+eta1*p1i_air(t5(3),t,n,b(2));
    end
else
    if aa4==0
        q11=q11+eta1*(q1i_air(t5(1),t5(2),t,n,b(1)));
        for i=1:1:aa3-1
        q12=q12+eta1*q2i_air(t5(i*2),t5(i*2+1),t,n,b(i),b(i+1));
        q13=q13+q3i_air(t5(i*2),t5(i*2+1),t,n,b(i),b(i+1));
        q11=q11+eta1*(q1i_air(t5(i*2+1),t5(i*2+2),t,n,b(i+1)));
        end
        q1=q11+q12+q13+eta1*p2i_air(t5(aa2),t5(aa2+1),t,n,b(aa3),b(aa3+1))+p3i_air(t5(aa2),t5(aa2+1),t,n,b(aa3),b(aa3+1));
    else
        q11=q11+eta1*(q1i_air(t5(1),t5(2),t,n,b(1)));
        for i=1:1:aa3-1
        q12=q12+eta1*q2i_air(t5(i*2),t5(i*2+1),t,n,b(i),b(i+1));
        q13=q13+q3i_air(t5(i*2),t5(i*2+1),t,n,b(i),b(i+1));
        q11=q11+eta1*(q1i_air(t5(i*2+1),t5(i*2+2),t,n,b(i+1)));
        end
        q12=q12+eta1*q2i_air(t5(aa2-1),t5(aa2),t,n,b(aa3),b(aa3+1));
        q13=q13+q3i_air(t5(aa2-1),t5(aa2),t,n,b(aa3),b(aa3+1));
        q1=q11+q12+q13+eta1*p1i_air(t5(aa2),t,n,b(aa3+1));
    end
end

fi_n=(2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n);
q4=2*L/(n*pi)*(-1)^n*q1;
q=fi_n*exp(-(D_d1*(n*pi/L)^2+eta1)*t)+q4;
end

%Below is the adaptive medication function, situation 5 of drug administration.
function q = Drug_5(x,t) 
global L
global t51 t52 t53 t54 t55 t56 t57 t58 
global inter1
global b1 b2 b3 b4
global b t5
d=inter1;
aa1=find(t>=t5);
aa2=aa1(end);
aa3=floor(aa2/2);
aa4=aa2-aa3*2;

if aa3==0
    dnt=b(1);
else
    if aa4==0
        dnt=b(aa3)/d^2*(1+2/d*(t-t5(aa2)))*(t-t5(aa2+1))^2+b(aa3+1)/d^2*(1-2/d*(t-t5(aa2+1)))*(t-t5(aa2))^2;
    else
        dnt=b(aa3+1);
    end
end

q0=0;
    for n=1:1:20
    q0=q0+Tn5(t,n)*sin(n*pi*x/L);  
    end
q=q0/x+dnt;
end