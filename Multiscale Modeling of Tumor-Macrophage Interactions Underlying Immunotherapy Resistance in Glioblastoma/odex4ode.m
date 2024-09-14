function dudt_ode = odex4ode(x,t,u)
global V7 a V81 V82 V9 V10 d7 d8 d9 d10 K71 K72 K81 K82 K83 K91 K92 K101 K102 K73 K93 n
global select_EGFR_I select_IGF1R_I dmax2 dmax3 delt_t TimeLength
    nx=length(x);
    EGFR_I=zeros(1,nx);
    IGF1R_I=zeros(1,nx);
    for i=1:1:nx
    EGFR_I(i)=Drug2(x(i),t,select_EGFR_I,dmax2);
    IGF1R_I(i)=Drug2(x(i),t,select_IGF1R_I,dmax3);
    end

    up=zeros(4,nx);
    for jj=1:1:nx
    up(1,jj) = V7*u(5,jj)/5/(K71+u(5,jj)/5)/(1+u(8,jj)/K72)/(1+EGFR_I(jj)/K73)*(1-u(7,jj))-d7*u(7,jj);
    up(2,jj) = a*(1+V81*u(7,jj)^n/(K81^n+u(7,jj)^n))*(1+V82*u(9,jj)/(K82+u(9,jj)))/(1+u(10,jj)/K83)*(1-u(8,jj))-d8*u(8,jj); 
    up(3,jj) = V9*u(6,jj)/3000/(K91+u(6,jj)/3000)/(1+u(8,jj)/K92)/(1+IGF1R_I(jj)/K93)*(1-u(9,jj))-d9*u(9,jj);
    up(4,jj) = V10*u(7,jj)/(K101+u(7,jj))*u(9,jj)/(K102+u(9,jj))*(1-u(10,jj))-d10*u(10,jj); 
    end
    dudt_ode=up;
end