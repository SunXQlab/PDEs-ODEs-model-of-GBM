function aaa=A_Drugsimulation2(x,t,select,dmax,dmax_1)

global d_max m0 b1 b2

    d_max=dmax;
    m0=dmax;

if x==0
    x=0.0000001;
end
switch select
    case 0
        aaa=0;
    case 1
        aaa=A_1(x,t);
    case 2
        aaa=A_2(x,t);
    case 3
        aaa=A_3(x,t);
    case 4
        aaa=A_4(x,t);
    case 5
        aaa=A_5(x,t);
end
end
%The following is the calculation of the concentration of A produced by the drug CSF1R_I at position x and time t.

% constant function
function q = TnA_1(t,n)
global d_max
global D_d1;
global eta1;
global L
C_d=d_max;
% d_initial=d_max/L;
fi_n=2*C_d*L/n/pi*(-1)^(n+1)+4*C_d*L/(n*pi)^3*((-1)^n-1)+C_d*(2*L)/n/pi*(-1)^n;
q=fi_n*(1-exp(-(D_d1*(n*pi/L)^2+eta1)*t))/(D_d1*(n*pi/L)^2+eta1)+2*eta1*C_d*(L/n/pi)*(-1)^n/(D_d1*(n*pi/L)^2+eta1)*(t+(exp(-(D_d1*(n*pi/L)^2+eta1)*t)-1)/(D_d1*(n*pi/L)^2+eta1));
end
% constant function
function q = A_1(x,t)
global d_max
global A_0
global L
global B_3  
q0=0;
C_d=d_max;
for n=1:1:100
    q0=q0+sin(n*pi/L*x)*TnA_1(t,n);
end
q=A_0+B_3/x*q0+B_3*C_d*t;
end

% Sine function
function q = TnA_2(t,n)
global d_max
global D_d1;
global eta1;
global L
global period_2
w=2*pi/period_2;
C_d=d_max;
C_d=C_d*2/3;
fi_n=(2*C_d*L/n/pi*(-1)^(n+1)+4*C_d*L/(n*pi)^3*((-1)^n-1)+C_d*(2*L)/n/pi*(-1)^n);
q1=fi_n*(1-exp(-(D_d1*(n*pi/L)^2+eta1)*t))/(D_d1*(n*pi/L)^2+eta1);
q2=2*eta1*C_d*(L/n/pi)*(-1)^n/(D_d1*(n*pi/L)^2+eta1)*(t+(exp(-(D_d1*(n*pi/L)^2+eta1)*t)-1)/(D_d1*(n*pi/L)^2+eta1));
q3=eta1*C_d*(D_d1*(n*pi/L)^2+eta1)*L/(n*pi)*(-1)^n/((D_d1*(n*pi/L)^2+eta1)^2+w^2)*(1-cos(w*t))/w;
q4=-eta1*C_d*L/(n*pi)*(-1)^n/((D_d1*(n*pi/L)^2+eta1)^2+w^2)*sin(w*t);
q5=eta1*C_d*w*L/(n*pi)*(-1)^n/((D_d1*(n*pi/L)^2+eta1)^2+w^2)*(1-exp(-(D_d1*(n*pi/L)^2+eta1)*t))/(D_d1*(n*pi/L)^2+eta1);
q6=C_d*(D_d1*(n*pi/L)^2+eta1)*L/(n*pi)*(-1)^n/((D_d1*(n*pi/L)^2+eta1)^2+w^2)*sin(w*t);
q7=-C_d*w*(D_d1*(n*pi/L)^2+eta1)*L/(n*pi)*(-1)^n/((D_d1*(n*pi/L)^2+eta1)^2+w^2)*(1-exp(-(D_d1*(n*pi/L)^2+eta1)*t))/(D_d1*(n*pi/L)^2+eta1);
q8=C_d*w*L/(n*pi)*(-1)^n/((D_d1*(n*pi/L)^2+eta1)^2+w^2)*(1-cos(w*t));
q=q1+q2+q3+q4+q5+q6+q7+q8;
end

% Sine function
function q = A_2(x,t)
global d_max
global A_0
global L
global B_3  
global period_2;
C_d=d_max;
C_d=C_d*2/3;
q0=0;
w=2*pi/period_2;
for n=1:1:100
    q0=q0+sin(n*pi*x/L)*TnA_2(t,n);
end
q=A_0+B_3/x*q0+B_3*C_d/2/w*(1-cos(w*t))+B_3*C_d*t;
end


% The periodic medication function.Situation 3.
function q=h_11_ceshiair2_11(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q1=2*eta1*L*(-1)^n/(n*pi)*m0/(D_d1*(n*pi/L)^2+eta1);
q2=s+exp(-(D_d1*(n*pi/L)^2+eta1)*(s-k*T))/(D_d1*(n*pi/L)^2+eta1);
q=q1*q2;
end
function q=h_12_ceshiair2_12(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q0=2*eta1*L*(-1)^n/(n*pi);
q11=m0/d^2*((s-(k+1/2)*T)^3+1/(2*d)*(s-(k+1/2)*T)^4);
q12=m0*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d))/(D_d1*(n*pi/L)^2+eta1);
q1=(q11+q12)/(D_d1*(n*pi/L)^2+eta1);
q21=2*m0/(3*d^3)*(s-(k+1/2)*T)^3+2*m0/d^2*(3/2*(s-(k+1/2)*T)^2+2/(3*d)*(s-(k+1/2)*T)^3);
q2=-q21/(D_d1*(n*pi/L)^2+eta1)^2;
q31=4*m0/d^3*(s-(k+1/2)*T)^2+2*m0/d^2*(s+1/d*(s-(k+1/2)*T+d)^2);
q32=-6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d))/(D_d1*(n*pi/L)^2+eta1);
q3=(q31+q32)/(D_d1*(n*pi/L)^2+eta1)^3;
q41=-12*m0/d^3*(exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d))/(D_d1*(n*pi/L)^2+eta1)+s);
q4=q41/(D_d1*(n*pi/L)^2+eta1)^4;
q=q0*(q1+q2+q3+q4);
end
function q=h_13_ceshiair2_13(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q0=2*eta1*L*(-1)^n/(n*pi);
q11=m0/d^2*((s-(k+1)*T+d)^3-1/(2*d)*(s-(k+1)*T+d)^4);
q1=q11/(D_d1*(n*pi/L)^2+eta1);
q21=-2*m0/(3*d^3)*(s-(k+1)*T+d)^3+2*m0/d^2*(3/2*(s-(k+1)*T+d)^2-2/(3*d)*(s-(k+1)*T+d)^3);
q2=-q21/(D_d1*(n*pi/L)^2+eta1)^2;
q31=-4*m0/d^3*(s-(k+1)*T+d)^2+2*m0/d^2*(s-1/d*(s-(k+1)*T)^2);
q32=6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T+d))/(D_d1*(n*pi/L)^2+eta1);
q3=(q31+q32)/(D_d1*(n*pi/L)^2+eta1)^3;
q41=12*m0/d^3*(s+exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T+d))/(D_d1*(n*pi/L)^2+eta1));
q4=q41/(D_d1*(n*pi/L)^2+eta1)^4;
q=q0*(q1+q2+q3+q4);
end
function q=h_14_ceshiair2_14(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q0=2*L*(-1)^n/(n*pi);
q11=2*m0/(3*d^3)*(s-(k+1/2)*T)^3+2*m0/d^2*(3/2*(s-(k+1/2)*T)^2+2/(3*d)*(s-(k+1/2)*T)^3);
q1=q11/(D_d1*(n*pi/L)^2+eta1);
q21=4*m0/d^3*(s-(k+1/2)*T)^2+2*m0/d^2*(s+1/d*(s-(k+1/2)*T+d)^2);
q22=-6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d))/(D_d1*(n*pi/L)^2+eta1);
q2=-(q21+q22)/(D_d1*(n*pi/L)^2+eta1)^2;
q31=12*m0/d^3*(s+exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d))/(D_d1*(n*pi/L)^2+eta1));
q3=q31/(D_d1*(n*pi/L)^2+eta1)^3;
q=q0*(q1+q2+q3);
end
function q=h_15_ceshiair2_15(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q0=2*L*(-1)^n/(n*pi);
q11=-2*m0/(3*d^3)*(s-(k+1)*T+d)^3+2*m0/d^2*(3/2*(s-(k+1)*T+d)^2-2/(3*d)*(s-(k+1)*T+d)^3);
q1=q11/(D_d1*(n*pi/L)^2+eta1);
q21=4*m0/d^3*(s-(k+1)*T+d)^2-2*m0/d^2*(s-1/d*(s-(k+1)*T)^2);
q22=-6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T+d))/(D_d1*(n*pi/L)^2+eta1);
q2=(q21+q22)/(D_d1*(n*pi/L)^2+eta1)^2;
q31=12*m0/d^3*(s+exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T+d))/(D_d1*(n*pi/L)^2+eta1));
q3=-q31/(D_d1*(n*pi/L)^2+eta1)^3;
q=q0*(q1+q2+q3);
end
% g1
function q=g_11_ceshiair2_11(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q1=2*eta1*L*(-1)^n/(n*pi)*m0/(D_d1*(n*pi/L)^2+eta1)^2;
q2=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-k*T))-exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d));
q=q1*q2;
end
% g2
function q=g_12_ceshiair2_12(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q0=2*eta1*L*(-1)^n/(n*pi);
q11=m0*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d));
q1=q11/(D_d1*(n*pi/L)^2+eta1)^2;
q21=6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T));
q22=6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d));
q2=-(q21+q22)/(D_d1*(n*pi/L)^2+eta1)^4;
q31=12*m0/d^3*(exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d))-exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T)));
q3=-q31/(D_d1*(n*pi/L)^2+eta1)^5;
q=q0*(q1+q2+q3);
end
% g3
function q=g_13_ceshiair2_13(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q0=2*eta1*L*(-1)^n/(n*pi);
q11=m0*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T));
q1=-q11/(D_d1*(n*pi/L)^2+eta1)^2;
q21=6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T+d));
q22=6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T));
q2=(q21+q22)/(D_d1*(n*pi/L)^2+eta1)^4;
q31=12*m0/d^3*(exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T))-exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T+d)));
q3=-q31/(D_d1*(n*pi/L)^2+eta1)^5;
q=q0*(q1+q2+q3);
end
% g4
function q=g_14_ceshiair2_14(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q0=2*L*(-1)^n/(n*pi);
q11=6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T));
q12=6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d));
q1=(q11+q12)/(D_d1*(n*pi/L)^2+eta1)^3;
q21=12*m0/d^3*(exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T))-exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1/2)*T+d)));
q2=-q21/(D_d1*(n*pi/L)^2+eta1)^4;
q=q0*(q1+q2);
end
% g5
function q=g_15_ceshiair2_15(k,s,n)
global D_d1
global L
global inter1
global period1
global m0
global eta1
T=2*period1;
d=inter1;
q0=2*L*(-1)^n/(n*pi);
q11=6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T));
q12=6*m0/d^2*exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T+d));
q1=-(q11+q12)/(D_d1*(n*pi/L)^2+eta1)^3;
q21=12*m0/d^3*(exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T))-exp(-(D_d1*(n*pi/L)^2+eta1)*(s-(k+1)*T+d)));
q2=q21/(D_d1*(n*pi/L)^2+eta1)^4;
q=q0*(q1+q2);
end

function q1=h_7_ceshiair7(s)
global m0
q1=m0*s;
end

function q5=h_8_ceshiair8(k,s)
global inter1
global period1
global m0
T=2*period1;
d=inter1;
q1=m0/d^2/2/d*s^4;
q2=m0/d^2/3*(3-2/d*(k+1/2)*T-4/d*(k+1/2)*T)*s^3;
q3=m0/d^2/2*(-6+6/d*(k+1/2)*T)*(k+1/2)*T*s^2;
q4=m0/d^2*(3-2/d*(k+1/2)*T)*(k+1/2)^2*T^2*s;
q5=q1+q2+q3+q4;
end

function q5=h_9_ceshiair9(k,s)
global inter1
global period1
global m0
T=2*period1;
d=inter1;
q1=-m0/d^2/2/d*s^4;
q2=m0/d^2/3*(1+2/d*(k+1)*T-4/d*(d-(k+1)*T))*s^3;
q3=m0/d^2/2*(-2/d*(d-(k+1)*T)^2+2*(1+2/d*(k+1)*T)*(d-(k+1)*T))*s^2;
q4=m0/d^2*(1+2/d*(k+1)*T)*(d-(k+1)*T)^2*s;
q5=q1+q2+q3+q4;
end

% The periodic medication function.Situation 3.
function q = TnA_3(t,n)
global D_d1;
global eta1;
global L
global period1  
global inter1
global m0
d=inter1;
T=2*period1;
d1=floor(t/period1);
d2=floor(d1/2);
d3=rem(d1,2);
t0=t-d1*period1;
q=0;
q4=0;
d_0=m0;
d_initial=d_0/L;
if d2~=0
for i=0:1:(d2-1)
    s=i*T;
    q1=h_11_ceshiair2_11(i,s,n);
    s=(i+1/2)*T-d;
    q2=h_11_ceshiair2_11(i,s,n);
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1/2)*T-d;
    q1=h_12_ceshiair2_12(i,s,n);
    s=(i+1/2)*T;
    q2=h_12_ceshiair2_12(i,s,n);
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1)*T-d;
    q1=h_13_ceshiair2_13(i,s,n);
    s=(i+1)*T;
    q2=h_13_ceshiair2_13(i,s,n);
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1/2)*T-d;
    q1=h_14_ceshiair2_14(i,s,n);
    s=(i+1/2)*T;
    q2=h_14_ceshiair2_14(i,s,n);
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1)*T-d;
    q1=h_15_ceshiair2_15(i,s,n);
    s=(i+1)*T;
    q2=h_15_ceshiair2_15(i,s,n);
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1/2)*T-d;
    q1=g_11_ceshiair2_11(i,s,n);
    s=t;
    q2=g_11_ceshiair2_11(i,s,n);
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1/2)*T;
    q1=g_12_ceshiair2_12(i,s,n);
    s=t;
    q2=g_12_ceshiair2_12(i,s,n);
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1)*T;
    q1=g_13_ceshiair2_13(i,s,n);
    s=t;
    q2=g_13_ceshiair2_13(i,s,n);
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1/2)*T;
    q1=g_14_ceshiair2_14(i,s,n);
    s=t;
    q2=g_14_ceshiair2_14(i,s,n);
    q3=q2-q1;
    q4=q4+q3;
    s=(i+1)*T;
    q1=g_15_ceshiair2_15(i,s,n);
    s=t;
    q2=g_15_ceshiair2_15(i,s,n);
    q3=q2-q1;
    q4=q4+q3;
end
end

if d3==0 && t0<(period1-inter1)
    s=d2*T;
    q1=h_11_ceshiair2_11(d2,s,n);
    s=t;
    q2=h_11_ceshiair2_11(d2,s,n);
    q3=q2-q1;
    q4=q4+q3; 
end
if d3==0 && t0>=(period1-inter1)
    s=d2*T;
    q1=h_11_ceshiair2_11(d2,s,n);
    s=(d2+1/2)*T-d;
    q2=h_11_ceshiair2_11(d2,s,n);
    q3=q2-q1;
    q4=q4+q3;  
    s=(d2+1/2)*T-d;
    q1=g_11_ceshiair2_11(d2,s,n);
    s=t;
    q2=g_11_ceshiair2_11(d2,s,n);
    q3=q2-q1;
    q4=q4+q3; 
    s=(d2+1/2)*T-d;
    q1=h_12_ceshiair2_12(d2,s,n);
    s=t;
    q2=h_12_ceshiair2_12(d2,s,n);
    q3=q2-q1;
    q4=q4+q3;
    s=(d2+1/2)*T-d;
    q1=h_14_ceshiair2_14(d2,s,n);
    s=t;
    q2=h_14_ceshiair2_14(d2,s,n);
    q3=q2-q1;
    q4=q4+q3;
end
if d3==1 && t0<(period1-inter1)
    s=d2*T;
    q1=h_11_ceshiair2_11(d2,s,n);
    s=(d2+1/2)*T-d;
    q2=h_11_ceshiair2_11(d2,s,n);
    q3=q2-q1;
    q4=q4+q3; 
    s=(d2+1/2)*T-d;
    q1=h_12_ceshiair2_12(d2,s,n);
    s=(d2+1/2)*T;
    q2=h_12_ceshiair2_12(d2,s,n);
    q3=q2-q1;
    q4=q4+q3;
    s=(d2+1/2)*T-d;
    q1=h_14_ceshiair2_14(d2,s,n);
    s=(d2+1/2)*T;
    q2=h_14_ceshiair2_14(d2,s,n);
    q3=q2-q1;
    q4=q4+q3;
    s=(d2+1/2)*T-d;
    q1=g_11_ceshiair2_11(d2,s,n);
    s=t;
    q2=g_11_ceshiair2_11(d2,s,n);
    q3=q2-q1;
    q4=q4+q3; 
    s=(d2+1/2)*T;
    q1=g_12_ceshiair2_12(d2,s,n);
    s=t;
    q2=g_12_ceshiair2_12(d2,s,n);
    q3=q2-q1;
    q4=q4+q3;
    s=(d2+1/2)*T;
    q1=g_14_ceshiair2_14(d2,s,n);
    s=t;
    q2=g_14_ceshiair2_14(d2,s,n);
    q3=q2-q1;
    q4=q4+q3;
end 
if d3==1 && t0>=(period1-inter1)
    s=d2*T;
    q1=h_11_ceshiair2_11(d2,s,n);
    s=(d2+1/2)*T-d;
    q2=h_11_ceshiair2_11(d2,s,n);
    q3=q2-q1;
    q4=q4+q3; 
    s=(d2+1/2)*T-d;
    q1=h_12_ceshiair2_12(d2,s,n);
    s=(d2+1/2)*T;
    q2=h_12_ceshiair2_12(d2,s,n);
    q3=q2-q1;
    q4=q4+q3;
    s=(d2+1/2)*T-d;
    q1=h_14_ceshiair2_14(d2,s,n);
    s=(d2+1/2)*T;
    q2=h_14_ceshiair2_14(d2,s,n);
    q3=q2-q1;
    q4=q4+q3;
    s=(d2+1/2)*T-d;
    q1=g_11_ceshiair2_11(d2,s,n);
    s=t;
    q2=g_11_ceshiair2_11(d2,s,n);
    q3=q2-q1;
    q4=q4+q3; 
    s=(d2+1/2)*T;
    q1=g_12_ceshiair2_12(d2,s,n);
    s=t;
    q2=g_12_ceshiair2_12(d2,s,n);
    q3=q2-q1;
    q4=q4+q3;
    s=(d2+1/2)*T;
    q1=g_14_ceshiair2_14(d2,s,n);
    s=t;
    q2=g_14_ceshiair2_14(d2,s,n);
    q3=q2-q1;
    q4=q4+q3;
    s=(d2+1)*T-d;
    q1=h_13_ceshiair2_13(d2,s,n);
    s=t;
    q2=h_13_ceshiair2_13(d2,s,n);
    q3=q2-q1;
    q4=q4+q3;
    s=(d2+1)*T-d;
    q1=h_15_ceshiair2_15(d2,s,n);
    s=t;
    q2=h_15_ceshiair2_15(d2,s,n);
    q3=q2-q1;
    q4=q4+q3;
end 

fi_n=2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n;
q5=fi_n*(1-exp(-(D_d1*(n*pi/L)^2+eta1)*t))/(D_d1*(n*pi/L)^2+eta1);
q=q4+q5;
end
% The periodic medication function.Situation 3.
function q = A_3(x,t)
global A_0
global L
global B_3   
global period1  
global inter1
global m0
d_0=m0;
d=inter1;
T=2*period1;
d1=floor(t/period1);
d2=floor(d1/2);
d3=rem(d1,2);
t0=t-d1*period1;
q0=0;
q4=0;
for k=0:1:d2-1
    s=k*T;
    q1=h_7_ceshiair7(s);
    s=(k+1/2)*T-d;
    q2=h_7_ceshiair7(s);
    q3=q2-q1;
    q4=q4+q3;
    s=(k+1/2)*T-d;
    q1=h_8_ceshiair8(k,s);
    s=(k+1/2)*T;
    q2=h_8_ceshiair8(k,s);
    q3=q2-q1;
    q4=q4+q3;
    s=(k+1)*T-d;
    q1=h_9_ceshiair9(k,s);
    s=(k+1)*T;
    q2=h_9_ceshiair9(k,s);
    q3=q2-q1;
    q4=q4+q3;    
end
if d3==0 && t0<=(period1-inter1)
    s=d2*T;
    q1=h_7_ceshiair7(s);
    s=t;
    q2=h_7_ceshiair7(s);
    q3=q2-q1;
    q4=q4+q3; 
end
if d3==0 && t0>(period1-inter1)
    s=d2*T;
    q1=h_7_ceshiair7(s);
    s=(d2+1/2)*T-d;
    q2=h_7_ceshiair7(s);
    q3=q2-q1;
    q4=q4+q3;  
    s=(d2+1/2)*T-d;
    q1=h_8_ceshiair8(d2,s);
    s=t;
    q2=h_8_ceshiair8(d2,s);
    q3=q2-q1;
    q4=q4+q3;    
end 
if d3==1 && t0<=(period1-inter1)
    s=d2*T;
    q1=h_7_ceshiair7(s);
    s=(d2+1/2)*T-d;
    q2=h_7_ceshiair7(s);
    q3=q2-q1;
    q4=q4+q3;  
    s=(d2+1/2)*T-d;
    q1=h_8_ceshiair8(d2,s);
    s=(d2+1/2)*T;
    q2=h_8_ceshiair8(d2,s);
    q3=q2-q1;
    q4=q4+q3; 
end 
if d3==1 && t0>(period1-inter1)
    s=d2*T;
    q1=h_7_ceshiair7(s);
    s=(d2+1/2)*T-d;
    q2=h_7_ceshiair7(s);
    q3=q2-q1;
    q4=q4+q3;  
    s=(d2+1/2)*T-d;
    q1=h_8_ceshiair8(d2,s);
    s=(d2+1/2)*T;
    q2=h_8_ceshiair8(d2,s);
    q3=q2-q1;
    q4=q4+q3; 
    s=(d2+1)*T-d;
    q1=h_9_ceshiair9(d2,s);
    s=t;
    q2=h_9_ceshiair9(d2,s);
    q3=q2-q1;
    q4=q4+q3; 
end 
for n=1:1:20
    q0=q0+TnA_3(t,n)*sin(n*pi*x/L);
end
q=A_0+B_3/x*q0+B_3*q4;
end


%Drug cessation function, situation 4.
function q=h_0_airjifen(s,n)
global D_d1
global L
global inter1
global eta1
global b1
d=inter1;
q11=b1/(D_d1*(n*pi/L)^2+eta1);
q12=s+exp(-(D_d1*(n*pi/L)^2+eta1)*s)/(D_d1*(n*pi/L)^2+eta1);
q=q11*q12;
end

function q=h_1_airjifen(t52,t53,s,n)
global D_d1
global L
global inter1
global period1
global eta1
global b1 b2
d=inter1;
q01=b1/(D_d1*(n*pi/L)^2+eta1);
q02=exp(-(D_d1*(n*pi/L)^2+eta1)*s)-exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52));
q03=1/(D_d1*(n*pi/L)^2+eta1);
q0=q01*q02*q03;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*b1/d^2*(t52-t53)^2;
q131=b1/d^2*(s-t53)^3/3;
q132=2*b1/d^3*(s^4/4-2*t53*s^3/3+t53^2*s^2/2-t52*s^3/3+t52*t53*s^2-t52*t53^2*s);
q133=b2/d^2*(s-t52)^3/3;
q134=-2*b2/d^3*(s^4/4-2*t52*s^3/3+t52^2*s^2/2-t53*s^3/3+t53*t52*s^2-t53*t52^2*s);
q1=q11*(q131+q132+q133+q134+q12/(D_d1*(n*pi/L)^2+eta1));
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(2*b1/d^3*(t52-t53)^2+2*b1/d^2*(t52-t53));
q231=2*b1/(3*d^3)*(s-t53)^3+b1/d^2*(s-t53)^2+4*b1/d^3*(s^3/3-(t52+t53)*s^2/2+t52*t53*s);
q232=-2*b2/(3*d^3)*(s-t52)^3+b2/d^2*(s-t52)^2-4*b2/d^3*(s^3/3-(t53+t52)*s^2/2+t53*t52*s);
q2=q21*(q231+q232+q22/(D_d1*(n*pi/L)^2+eta1));
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(8*b1/d^3*(t52-t53)+2*b1/d^2+2*b2/d^2*(1-2/d*(t52-t53)));
q33=4*b1/d^3*(s-t53)^2+2*b1/d^2*(s+(s-t52)^2/d)-4*b2/d^3*(s-t52)^2+2*b2/d^2*(s-(s-t53)^2/d);
q3=q31*(q33+q32/(D_d1*(n*pi/L)^2+eta1));
q41=-1/(D_d1*(n*pi/L)^2+eta1)^4;
q42=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(12*b1/d^3-12*b2/d^3);
q43=(12*b1/d^3-12*b2/d^3)*s;
q4=q41*(q43+q42/(D_d1*(n*pi/L)^2+eta1));
q=q0+q1+q2+q3+q4;
end

function q=h_2_airjifen(t52,t53,s,n)
global D_d1
global L
global inter1
global period1
global eta1
global b1 b2
d=inter1;
q01=b1/(D_d1*(n*pi/L)^2+eta1);
q02=exp(-(D_d1*(n*pi/L)^2+eta1)*s)-exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52));
q03=1/(D_d1*(n*pi/L)^2+eta1);
q0=q01*q02*q03;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*b1/d^2*(t52-t53)^2;
q13=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(b2/d^2*(t53-t52)^2);
q14=1/(D_d1*(n*pi/L)^2+eta1);
q1=q11*(q12-q13)*q14;
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(2*b1/d^3*(t52-t53)^2+2*b1/d^2*(t52-t53));
q23=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(-2*b2/d^3*(t53-t52)^2+2*b2/d^2*(t53-t52));
q24=1/(D_d1*(n*pi/L)^2+eta1);
q2=q21*(q22-q23)*q24;
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(8*b1/d^3*(t52-t53)+2*b1/d^2+2*b2/d^2*(1-2/d*(t52-t53)));
q33=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(2*b1/d^2*(1+2/d*(t53-t52))-8*b2/d^3*(t53-t52)+2*b2/d^2);
q34=1/(D_d1*(n*pi/L)^2+eta1);
q3=q31*(q32-q33)*q34;
q41=-1/(D_d1*(n*pi/L)^2+eta1)^4;
q42=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(12*b1/d^3-12*b2/d^3);
q43=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(12*b1/d^3-12*b2/d^3);
q44=1/(D_d1*(n*pi/L)^2+eta1);
q4=q41*(q42-q43)*q44;
q51=b2/(D_d1*(n*pi/L)^2+eta1);
q52=s+exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))/(D_d1*(n*pi/L)^2+eta1);
q5=q51*q52;
q=q0+q1+q2+q3+q4+q5;
end

function q=h_3_airjifen(t52,t53,s,n)
global D_d1
global L
global inter1
global period1
global eta1
global b1 b2
d=inter1;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(2*b1/d^3*(t52-t53)^2+2*b1/d^2*(t52-t53));
q131=2*b1/(3*d^3)*(s-t53)^3+b1/d^2*(s-t53)^2+4*b1/d^3*(s^3/3-(t52+t53)*s^2/2+t52*t53*s);
q132=-b2/(3*d^3)*(s-t52)^3+b2*(s-t52)^2-4*b2/d^3*(s^3/3-(t53+t52)*s^2/2+t53*t52*s);
q1=q11*(q131+q132+q12/(D_d1*(n*pi/L)^2+eta1));
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(8*b1/d^3*(t52-t53)+2*b1/d^2+2*b2/d^2*(1-2/d*(t52-t53)));
q23=4*b1/d^3*(s-t53)^2+2*b1/d^2*(s+(s-t52)^2/d)-4*b2/d^3*(s-t52)^2+2*b2/d^2*(s-(s-t53)^2/d);
q2=q21*(q23+q22/(D_d1*(n*pi/L)^2+eta1));
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(12*b1/d^3-12*b2/d^3);
q33=(12*b1/d^3-12*b2/d^3)*s;
q3=q31*(q33+q32/(D_d1*(n*pi/L)^2+eta1));
q=q1+q2+q3;
end
function q=h_4_airjifen(t52,t53,s,n)
global D_d1
global L
global inter1
global period1
global eta1
global b1 b2
d=inter1;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(2*b1/d^3*(t52-t53)^2+2*b1/d^2*(t52-t53));
q13=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(-2*b2/d^3*(t53-t52)^2+2*b2/d^2*(t53-t52));
q14=1/(D_d1*(n*pi/L)^2+eta1);
q1=q11*(q12-q13)*q14;
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(8*b1/d^3*(t52-t53)+2*b1/d^2+2*b2/d^2*(1-2/d*(t52-t53)));
q23=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(2*b1/d^2*(1+2/d*(t53-t52))-8*b2/d^3*(t53-t52)+2*b2/d^2);
q24=1/(D_d1*(n*pi/L)^2+eta1);
q2=q21*(q22-q23)*q24;
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(12*b1/d^3-12*b2/d^3);
q33=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(12*b1/d^3-12*b2/d^3);
q34=1/(D_d1*(n*pi/L)^2+eta1);
q3=q31*(q32-q33)*q34;
q=q1+q2+q3;
end
function q=h_5_airjifen(s)
global b1
q=b1*s;
end
function q=h_6_airjifen(t52,t53,s)
global inter1
global b1 b2
d=inter1;
q1=b1/d^2*((s-t53)^3/3+2/d*(s^4/4-2*t53*s^3/3+t53^2*s^2/2-t52*s^3/3+t52*t53*s^2-t52*t53^2*s));
q2=b2/d^2*((s-t52)^3/3-2/d*(s^4/4-2*t52*s^3/3+t52^2*s^2/2-t53*s^3/3+t53*t52*s^2-t53*t52^2*s));
q=q1+q2;
end
function q=h_7_airjifen(s)
global b2
q=b2*s;
end
%Drug cessation function, situation 4.
function q = TnA_4(t,n)
global D_d1;
global eta1;
global L
global b1
global t51 t52 t53
d_0=b1;
d_initial=d_0/L;
if t>=t51 && t<t52
    q4=eta1*(h_0_airjifen(t,n)-h_0_airjifen(0,n));
elseif t>=t52 && t<t53
    q41=eta1*(h_0_airjifen(t52,n)-h_0_airjifen(0,n));
    q42=eta1*(h_1_airjifen(t52,t53,t,n)-h_1_airjifen(t52,t53,t52,n));
    q43=h_3_airjifen(t52,t53,t,n)-h_3_airjifen(t52,t53,t52,n);
    q4=q41+q42+q43;
elseif t>=t53
    q41=eta1*(h_0_airjifen(t52,n)-h_0_airjifen(0,n));
    q42=eta1*(h_1_airjifen(t52,t53,t53,n)-h_1_airjifen(t52,t53,t52,n));
    q43=h_3_airjifen(t52,t53,t53,n)-h_3_airjifen(t52,t53,t52,n);
    q44=eta1*(h_2_airjifen(t52,t53,t,n)-h_2_airjifen(t52,t53,t53,n));
    q45=h_4_airjifen(t52,t53,t,n)-h_4_airjifen(t52,t53,t53,n);
    q4=q41+q42+q43+q44+q45;
else
    msg='Error occured';
    error(msg)    
end
q5=q4*2*L/n/pi*(-1)^n;
fi_n=2*d_0*L/n/pi*(-1)^(n+1)+4*d_0*L/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n;
q=fi_n*(1-exp(-(D_d1*(n*pi/L)^2+eta1)*t))/(D_d1*(n*pi/L)^2+eta1)+q5;
end
%Drug cessation function, situation 4.
function q = A_4(x,t)
global A_0
global L
global B_3  
global period1  
global inter1
global b1
global t51 t52 t53
d_0=b1;
d=inter1;
q0=0;
if t>=t51 && t<t52
    q41=h_5_airjifen(t51);
    q42=h_5_airjifen(t);
    q4=q42-q41;
elseif t>=t52 && t<t53
    q41=h_5_airjifen(t51);
    q42=h_5_airjifen(t52);
    q43=h_6_airjifen(t52,t53,t52);
    q44=h_6_airjifen(t52,t53,t);
    q4=q42-q41+q44-q43;
elseif t>=t53
    q41=h_5_airjifen(t51);
    q42=h_5_airjifen(t52);
    q43=h_6_airjifen(t52,t53,t52);
    q44=h_6_airjifen(t52,t53,t53);
    q45=h_7_airjifen(t53);
    q46=h_7_airjifen(t);
    q4=q42-q41+q44-q43+q46-q45;
else
    msg='Error occured';
    error(msg)    
end
for n=1:1:20
    q0=q0+TnA_4(t,n)*sin(n*pi*x/L);
end
q=A_0+B_3/x*q0+B_3*q4;
end




%Below is the adaptive medication function, situation 5 of drug administration.
function q=jifen_p1i_air(t51,s,n,bi)
global D_d1
global L
global inter1
global eta1
d=inter1;
q11=bi/(D_d1*(n*pi/L)^2+eta1);
q12=s+exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t51))/(D_d1*(n*pi/L)^2+eta1);
q=q11*q12;
end

function q=jifen_p2i_air(t52,t53,s,n,bi1,bi2)
global D_d1
global L
global inter1
global period1
global eta1
d=inter1;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*bi1/d^2*(t52-t53)^2;
q131=bi1/d^2*(s-t53)^3/3;
q132=2*bi1/d^3*(s^4/4-2*t53*s^3/3+t53^2*s^2/2-t52*s^3/3+t52*t53*s^2-t52*t53^2*s);
q133=bi2/d^2*(s-t52)^3/3;
q134=-2*bi2/d^3*(s^4/4-2*t52*s^3/3+t52^2*s^2/2-t53*s^3/3+t53*t52*s^2-t53*t52^2*s);
q1=q11*(q131+q132+q133+q134+q12/(D_d1*(n*pi/L)^2+eta1));
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(2*bi1/d^3*(t52-t53)^2+2*bi1/d^2*(t52-t53));
q231=2*bi1/(3*d^3)*(s-t53)^3+bi1/d^2*(s-t53)^2+4*bi1/d^3*(s^3/3-(t52+t53)*s^2/2+t52*t53*s);
q232=-2*bi2/(3*d^3)*(s-t52)^3+bi2/d^2*(s-t52)^2-4*bi2/d^3*(s^3/3-(t53+t52)*s^2/2+t53*t52*s);
q2=q21*(q231+q232+q22/(D_d1*(n*pi/L)^2+eta1));
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(8*bi1/d^3*(t52-t53)+2*bi1/d^2+2*bi2/d^2*(1-2/d*(t52-t53)));
q33=4*bi1/d^3*(s-t53)^2+2*bi1/d^2*(s+(s-t52)^2/d)-4*bi2/d^3*(s-t52)^2+2*bi2/d^2*(s-(s-t53)^2/d);
q3=q31*(q33+q32/(D_d1*(n*pi/L)^2+eta1));
q41=-1/(D_d1*(n*pi/L)^2+eta1)^4;
q42=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(12*bi1/d^3-12*bi2/d^3);
q43=(12*bi1/d^3-12*bi2/d^3)*s;
q4=q41*(q43+q42/(D_d1*(n*pi/L)^2+eta1));
q=q1+q2+q3+q4;
end
function q=jifen_p3i_air(t52,t53,s,n,bi1,bi2)
global D_d1
global L
global inter1
global period1
global eta1
d=inter1;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(2*bi1/d^3*(t52-t53)^2+2*bi1/d^2*(t52-t53));
q131=2*bi1/(3*d^3)*(s-t53)^3+bi1/d^2*(s-t53)^2+4*bi1/d^3*(s^3/3-(t52+t53)*s^2/2+t52*t53*s);
q132=-bi2/(3*d^3)*(s-t52)^3+bi2*(s-t52)^2-4*bi2/d^3*(s^3/3-(t53+t52)*s^2/2+t53*t52*s);
q1=q11*(q131+q132+q12/(D_d1*(n*pi/L)^2+eta1));
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(8*bi1/d^3*(t52-t53)+2*bi1/d^2+2*bi2/d^2*(1-2/d*(t52-t53)));
q23=4*bi1/d^3*(s-t53)^2+2*bi1/d^2*(s+(s-t52)^2/d)-4*bi2/d^3*(s-t52)^2+2*bi2/d^2*(s-(s-t53)^2/d);
q2=q21*(q23+q22/(D_d1*(n*pi/L)^2+eta1));
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(12*bi1/d^3-12*bi2/d^3);
q33=(12*bi1/d^3-12*bi2/d^3)*s;
q3=q31*(q33+q32/(D_d1*(n*pi/L)^2+eta1));
q=q1+q2+q3;
end
function q=jifen_q1i_air(t51,t52,s,n,bi)
global D_d1
global L
global inter1
global eta1
d=inter1;
q11=bi/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t51))-exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52));
q13=1/(D_d1*(n*pi/L)^2+eta1);
q=q11*q12*q13;
end

function q=jifen_q2i_air(t52,t53,s,n,bi1,bi2)
global D_d1
global L
global inter1
global period1
global eta1
d=inter1;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*bi1/d^2*(t52-t53)^2;
q13=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(bi2/d^2*(t53-t52)^2);
q14=1/(D_d1*(n*pi/L)^2+eta1);
q1=q11*(q12-q13)*q14;
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(2*bi1/d^3*(t52-t53)^2+2*bi1/d^2*(t52-t53));
q23=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(-2*bi2/d^3*(t53-t52)^2+2*bi2/d^2*(t53-t52));
q24=1/(D_d1*(n*pi/L)^2+eta1);
q2=q21*(q22-q23)*q24;
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(8*bi1/d^3*(t52-t53)+2*bi1/d^2+2*bi2/d^2*(1-2/d*(t52-t53)));
q33=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(2*bi1/d^2*(1+2/d*(t53-t52))-8*bi2/d^3*(t53-t52)+2*bi2/d^2);
q34=1/(D_d1*(n*pi/L)^2+eta1);
q3=q31*(q32-q33)*q34;
q41=-1/(D_d1*(n*pi/L)^2+eta1)^4;
q42=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(12*bi1/d^3-12*bi2/d^3);
q43=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(12*bi1/d^3-12*bi2/d^3);
q44=1/(D_d1*(n*pi/L)^2+eta1);
q4=q41*(q42-q43)*q44;
q=q1+q2+q3+q4;
end
function q=jifen_q3i_air(t52,t53,s,n,bi1,bi2)
global D_d1
global L
global inter1
global period1
global eta1
d=inter1;
q11=1/(D_d1*(n*pi/L)^2+eta1);
q12=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(2*bi1/d^3*(t52-t53)^2+2*bi1/d^2*(t52-t53));
q13=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(-2*bi2/d^3*(t53-t52)^2+2*bi2/d^2*(t53-t52));
q14=1/(D_d1*(n*pi/L)^2+eta1);
q1=q11*(q12-q13)*q14;
q21=-1/(D_d1*(n*pi/L)^2+eta1)^2;
q22=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(8*bi1/d^3*(t52-t53)+2*bi1/d^2+2*bi2/d^2*(1-2/d*(t52-t53)));
q23=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(2*bi1/d^2*(1+2/d*(t53-t52))-8*bi2/d^3*(t53-t52)+2*bi2/d^2);
q24=1/(D_d1*(n*pi/L)^2+eta1);
q2=q21*(q22-q23)*q24;
q31=1/(D_d1*(n*pi/L)^2+eta1)^3;
q32=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t52))*(12*bi1/d^3-12*bi2/d^3);
q33=exp(-(D_d1*(n*pi/L)^2+eta1)*(s-t53))*(12*bi1/d^3-12*bi2/d^3);
q34=1/(D_d1*(n*pi/L)^2+eta1);
q3=q31*(q32-q33)*q34;
q=q1+q2+q3;
end

function q=jifen_d1i_air(s,bi)
global D_d1
global L
global inter1
d=inter1;
q=bi*s;
end

function q=jifen_d2i_air(t52,t53,s,bi1,bi2)
global D_d1
global L
global inter1
d=inter1;
q1=bi1/d^2*((s-t53)^3/3+2/d*(s^4/4-2*t53*s^3/3+t53^2*s^2/2-t52*s^3/3+t52*t53*s^2-t52*t53^2*s));
q2=bi2/d^2*((s-t52)^3/3-2/d*(s^4/4-2*t52*s^3/3+t52^2*s^2/2-t53*s^3/3+t53*t52*s^2-t53*t52^2*s));
q=q1+q2;
end

function q = TnA_5(t,n)
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
q111=0;q112=0;q121=0;q122=0;q131=0;q132=0;

if aa3==0
        q1=eta1*(jifen_p1i_air(t5(1),t,n,b(1))-jifen_p1i_air(t5(1),t5(1),n,b(1)));
elseif aa3==1
    if aa4==0
    q111=eta1*(jifen_p1i_air(t5(1),t5(2),n,b(1))-jifen_p1i_air(t5(1),t5(1),n,b(1)));
    q112=eta1*(jifen_q1i_air(t5(1),t5(2),t,n,b(1))-jifen_q1i_air(t5(1),t5(2),t5(2),n,b(1)));
    q12=eta1*(jifen_p2i_air(t5(2),t5(3),t,n,b(1),b(2))-jifen_p2i_air(t5(2),t5(3),t5(2),n,b(1),b(2)));
    q13=(jifen_p3i_air(t5(2),t5(3),t,n,b(1),b(2))-jifen_p3i_air(t5(2),t5(3),t5(2),n,b(1),b(2)));
    q1=q111+q112+q12+q13;
    else
    q111=eta1*(jifen_p1i_air(t5(1),t5(2),n,b(1))-jifen_p1i_air(t5(1),t5(1),n,b(1)));
    q112=eta1*(jifen_q1i_air(t5(1),t5(2),t,n,b(1))-jifen_q1i_air(t5(1),t5(2),t5(2),n,b(1)));
    q113=eta1*(jifen_p1i_air(t5(3),t,n,b(2))-jifen_p1i_air(t5(3),t5(3),n,b(2)));
    q121=eta1*(jifen_p2i_air(t5(2),t5(3),t5(3),n,b(1),b(2))-jifen_p2i_air(t5(2),t5(3),t5(2),n,b(1),b(2)));
    q122=eta1*(jifen_q2i_air(t5(2),t5(3),t,n,b(1),b(2))-jifen_q2i_air(t5(2),t5(3),t5(3),n,b(1),b(2)));
    q131=(jifen_p3i_air(t5(2),t5(3),t5(3),n,b(1),b(2))-jifen_p3i_air(t5(2),t5(3),t5(2),n,b(1),b(2)));
    q132=(jifen_q3i_air(t5(2),t5(3),t,n,b(1),b(2))-jifen_q3i_air(t5(2),t5(3),t5(3),n,b(1),b(2)));
    q1=q111+q112+q113+q121+q122+q131+q132;
    end
else
    if aa4==0
        q111=q111+eta1*(jifen_p1i_air(t5(1),t5(2),n,b(1))-jifen_p1i_air(t5(1),t5(1),n,b(1)));
        q112=q112+eta1*(jifen_q1i_air(t5(1),t5(2),t,n,b(1))-jifen_q1i_air(t5(1),t5(2),t5(2),n,b(1)));
        for i=1:1:aa3-1
            q121=q121+eta1*(jifen_p2i_air(t5(i*2),t5(i*2+1),t5(i*2+1),n,b(i),b(i+1))-jifen_p2i_air(t5(i*2),t5(i*2+1),t5(i*2),n,b(i),b(i+1)));
            q122=q122+eta1*(jifen_q2i_air(t5(i*2),t5(i*2+1),t,n,b(i),b(i+1))-jifen_q2i_air(t5(i*2),t5(i*2+1),t5(i*2+1),n,b(i),b(i+1)));
            q131=q131+(jifen_p3i_air(t5(i*2),t5(i*2+1),t5(i*2+1),n,b(i),b(i+1))-jifen_p3i_air(t5(i*2),t5(i*2+1),t5(i*2),n,b(i),b(i+1)));
            q132=q132+(jifen_q3i_air(t5(i*2),t5(i*2+1),t,n,b(i),b(i+1))-jifen_q3i_air(t5(i*2),t5(i*2+1),t5(i*2+1),n,b(i),b(i+1)));
            q111=q111+eta1*(jifen_p1i_air(t5(i*2+1),t5(i*2+2),n,b(i+1))-jifen_p1i_air(t5(i*2+1),t5(i*2+1),n,b(i+1)));
            q112=q112+eta1*(jifen_q1i_air(t5(i*2+1),t5(i*2+2),t,n,b(i+1))-jifen_q1i_air(t5(i*2+1),t5(i*2+2),t5(i*2+2),n,b(i+1)));
        end
            q123=eta1*(jifen_p2i_air(t5(aa2),t5(aa2+1),t,n,b(aa3),b(aa3+1))-jifen_p2i_air(t5(aa2),t5(aa2+1),t5(aa2),n,b(aa3),b(aa3+1)));
            q133=(jifen_p3i_air(t5(aa2),t5(aa2+1),t,n,b(aa3),b(aa3+1))-jifen_p3i_air(t5(aa2),t5(aa2+1),t5(aa2),n,b(aa3),b(aa3+1)));
            q1=q111+q112+q121+q122+q123+q131+q132+q133;
    else
        q111=q111+eta1*(jifen_p1i_air(t5(1),t5(2),n,b(1))-jifen_p1i_air(t5(1),t5(1),n,b(1)));
        q112=q112+eta1*(jifen_q1i_air(t5(1),t5(2),t,n,b(1))-jifen_q1i_air(t5(1),t5(2),t5(2),n,b(1)));
        for i=1:1:aa3-1         
            q121=q121+eta1*(jifen_p2i_air(t5(i*2),t5(i*2+1),t5(i*2+1),n,b(i),b(i+1))-jifen_p2i_air(t5(i*2),t5(i*2+1),t5(i*2),n,b(i),b(i+1)));
            q122=q122+eta1*(jifen_q2i_air(t5(i*2),t5(i*2+1),t,n,b(i),b(i+1))-jifen_q2i_air(t5(i*2),t5(i*2+1),t5(i*2+1),n,b(i),b(i+1)));
            q131=q131+(jifen_p3i_air(t5(i*2),t5(i*2+1),t5(i*2+1),n,b(i),b(i+1))-jifen_p3i_air(t5(i*2),t5(i*2+1),t5(i*2),n,b(i),b(i+1)));
            q132=q132+(jifen_q3i_air(t5(i*2),t5(i*2+1),t,n,b(i),b(i+1))-jifen_q3i_air(t5(i*2),t5(i*2+1),t5(i*2+1),n,b(i),b(i+1)));
            q111=q111+eta1*(jifen_p1i_air(t5(i*2+1),t5(i*2+2),n,b(i+1))-jifen_p1i_air(t5(i*2+1),t5(i*2+1),n,b(i+1)));
            q112=q112+eta1*(jifen_q1i_air(t5(i*2+1),t5(i*2+2),t,n,b(i+1))-jifen_q1i_air(t5(i*2+1),t5(i*2+2),t5(i*2+2),n,b(i+1)));
        end
            q121=q121+eta1*(jifen_p2i_air(t5(aa3*2),t5(aa3*2+1),t5(aa3*2+1),n,b(aa3),b(aa3+1))-jifen_p2i_air(t5(aa3*2),t5(aa3*2+1),t5(aa3*2),n,b(aa3),b(aa3+1)));
            q122=q122+eta1*(jifen_q2i_air(t5(aa3*2),t5(aa3*2+1),t,n,b(aa3),b(aa3+1))-jifen_q2i_air(t5(aa3*2),t5(aa3*2+1),t5(aa3*2+1),n,b(aa3),b(aa3+1)));
            q131=q131+(jifen_p3i_air(t5(aa3*2),t5(aa3*2+1),t5(aa3*2+1),n,b(aa3),b(aa3+1))-jifen_p3i_air(t5(aa3*2),t5(aa3*2+1),t5(aa3*2),n,b(aa3),b(aa3+1)));
            q132=q132+(jifen_q3i_air(t5(aa3*2),t5(aa3*2+1),t,n,b(aa3),b(aa3+1))-jifen_q3i_air(t5(aa3*2),t5(aa3*2+1),t5(aa3*2+1),n,b(aa3),b(aa3+1)));
            q113=eta1*(jifen_p1i_air(t5(aa2),t,n,b(aa3+1))-jifen_p1i_air(t5(aa2),t5(aa2),n,b(aa3+1)));
            q1=q111+q112+q113+q121+q122+q131+q132;
    end
end

fi_n=(2*d_initial*L^2/n/pi*(-1)^(n+1)+4*d_initial*L^2/(n*pi)^3*((-1)^n-1)+d_0*(2*L)/n/pi*(-1)^n);
q4=2*L/(n*pi)*(-1)^n*q1;
q=fi_n*(1-exp(-(D_d1*(n*pi/L)^2+eta1)*t))/(D_d1*(n*pi/L)^2+eta1)+q4;
end

%Below is the adaptive medication function, situation 5 of drug administration.
function q = A_5(x,t) 
global L
global t51 t52 t53 t54 t55 t56 t57 t58 
global b1 b2 b3 b4
global A_0
global B_3 
global b t5

aa1=find(t>=t5);
aa2=aa1(end);
aa3=floor(aa2/2);
aa4=aa2-aa3*2;

q11=0;q12=0;
if aa3==0
        q1=jifen_d1i_air(t,b(1))-jifen_d1i_air(t5(1),b(1));
elseif aa3==1
    if aa4==0
        q11=jifen_d1i_air(t5(2),b(1))-jifen_d1i_air(t5(1),b(1));
        q12=jifen_d2i_air(t5(2),t5(3),t,b(1),b(2))-jifen_d2i_air(t5(2),t5(3),t5(2),b(1),b(2));
        q1=q11+q12;
    else
        q11=jifen_d1i_air(t5(2),b(1))-jifen_d1i_air(t5(1),b(1));
        q12=jifen_d2i_air(t5(2),t5(3),t5(3),b(1),b(2))-jifen_d2i_air(t5(2),t5(3),t5(2),b(1),b(2));
        q13=jifen_d1i_air(t,b(2))-jifen_d1i_air(t5(3),b(2));
        q1=q11+q12+q13;
    end
else
    if aa4==0
        q11=q11+jifen_d1i_air(t5(2),b(1))-jifen_d1i_air(t5(1),b(1));
        for i=1:1:aa3-1
            q12=q12+jifen_d2i_air(t5(i*2),t5(i*2+1),t5(i*2+1),b(i),b(i+1))-jifen_d2i_air(t5(i*2),t5(i*2+1),t5(i*2),b(i),b(i+1));
            q11=q11+jifen_d1i_air(t5(i*2+2),b(i+1))-jifen_d1i_air(t5(i*2+1),b(i+1));
        end
        q12=q12+jifen_d2i_air(t5(aa2),t5(aa2+1),t,b(aa3),b(aa3+1))-jifen_d2i_air(t5(aa2),t5(aa2+1),t5(aa2),b(aa3),b(aa3+1));
        q1=q11+q12;
    else
        q11=q11+jifen_d1i_air(t5(2),b(1))-jifen_d1i_air(t5(1),b(1));
        for i=1:1:aa3-1
            q12=q12+jifen_d2i_air(t5(i*2),t5(i*2+1),t5(i*2+1),b(i),b(i+1))-jifen_d2i_air(t5(i*2),t5(i*2+1),t5(i*2),b(i),b(i+1));
            q11=q11+jifen_d1i_air(t5(i*2+2),b(i+1))-jifen_d1i_air(t5(i*2+1),b(i+1));
        end
        q12=q12+jifen_d2i_air(t5(aa3*2),t5(aa3*2+1),t5(aa3*2+1),b(aa3),b(aa3+1))-jifen_d2i_air(t5(aa3*2),t5(aa3*2+1),t5(aa3*2),b(aa3),b(aa3+1));
        q11=q11+jifen_d1i_air(t,b(aa3+1))-jifen_d1i_air(t5(aa2),b(aa3+1));
        q1=q11+q12;
    end
end
dnt=q1;

q0=0;
    for n=1:1:20
    q0=q0+TnA_5(t,n)*sin(n*pi*x/L);  
    end
q=A_0+B_3/x*q0+B_3*dnt;
end