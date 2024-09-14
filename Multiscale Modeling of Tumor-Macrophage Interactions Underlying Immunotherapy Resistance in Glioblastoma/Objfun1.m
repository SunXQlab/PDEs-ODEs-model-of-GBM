function ObjV=Objfun1(rvChrom)
[m,n]=size(rvChrom);
global time
global Ptest1     
ObjV=zeros(m,1);
for iii=1:m

Ptest1=rvChrom(iii,:);
Conditions_ODES_2;
%%%%%
time=time+1;

ObjV(iii)=sum(sum((Y2(:,1)-pY1068_EGFR).^2))+sum(sum((Y2(:,2)-p_ERK).^2))+sum(sum((Y2(:,3)-p_IGF1Rbeta).^2))+sum(sum((Y2(:,4)-pS473_AKT).^2));   % y is simulated data; y0 is experimental data

end
% plot(ObjV)