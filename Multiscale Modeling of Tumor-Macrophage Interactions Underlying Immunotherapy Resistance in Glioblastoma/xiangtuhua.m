   load('CSF1RI_1_a_ERK_1_a_A_0.05_grow.mat','u1','u2','u3','u4','u5','u6','u7','u8','u9','u10','C_T_mean','M1_mean','M2_mean');
   C_T_mean1=C_T_mean;M1_mean1=M1_mean;M2_mean1=M2_mean;
   [aaa1,aaa2]=min(C_T_mean1(25*24+1:200*24+1,1));
   aaa3=25*24+aaa2;
   aaa4=find(C_T_mean1(25*24+1:200*24+1,1)<=0.01);
   aaa5=25*24+aaa4(1,1);
   load('CSF1RI_1_a_ERK_1_a_A_2.5_grow.mat','u1','u2','u3','u4','u5','u6','u7','u8','u9','u10','C_T_mean','M1_mean','M2_mean');
   C_T_mean2=C_T_mean;M1_mean2=M1_mean;M2_mean2=M2_mean;
   [bbb1,bbb2]=min(C_T_mean2(25*24+1:200*24+1,1));
   bbb3=25*24+bbb2;
    


  
   %%%Elephant Diagram
color_r=zeros(4801,1);
color_g=zeros(4801,1);
color_b=zeros(4801,1);
% c1=[250 227 113];c2=[137 157 192];c3=[217 104 49];c4=[56 13 48];
c1=[0.1 0.7 0.1]*255;c2=[137 157 192];c3=[0.201 0.5 1]*255;c4=[0.2 0.2 0.7]*255;
% c1=[1 1 0.1]*255;c2=[0.8 0.8 0.2]*255;c3=[0 0.7 0.7]*255;c4=[0.2 0.2 0.71]*255;
color_r(1:601,1)=c1(1)/255:(c2(1)-c1(1))/255/600:c2(1)/255;
color_g(1:601,1)=c1(2)/255:(c2(2)-c1(2))/255/600:c2(2)/255;
color_b(1:601,1)=c1(3)/255:(c2(3)-c1(3))/255/600:c2(3)/255;
color_r(601:bbb3,1)=c2(1)/255:(c3(1)-c2(1))/255/(bbb3-601):c3(1)/255;
color_g(601:bbb3,1)=c2(2)/255:(c3(2)-c2(2))/255/(bbb3-601):c3(2)/255;
color_b(601:bbb3,1)=c2(3)/255:(c3(3)-c2(3))/255/(bbb3-601):c3(3)/255;
color_r(bbb3:4801,1)=c3(1)/255:(c4(1)-c3(1))/255/(4801-bbb3):c4(1)/255;
color_g(bbb3:4801,1)=c3(2)/255:(c4(2)-c3(2))/255/(4801-bbb3):c4(2)/255;
color_b(bbb3:4801,1)=c3(3)/255:(c4(3)-c3(3))/255/(4801-bbb3):c4(3)/255;
c=zeros(4801,3);
c(:,1)=color_r;
c(:,2)=color_g;
c(:,3)=color_b;
%Responder
figure
sz=25;
scatter(M1_mean1,C_T_mean1,sz,c,'filled');
map = colormap(c);
colorbar;
   caxis([0 200]);
   h=colorbar('xtick',[0,50,100,150,200]);
   xlim([0 1.2]),ylim([0 1]);
   set(gca,'xtick',0:0.2:1.2); set(gca,'xticklabel',0:0.2:1.2);
   xlabel('Average M1 macrophage density','FontWeight','Bold','FontSize',18,'FontName','Arial');
   ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18,'FontName','Arial');
%    title('Positive patient','FontWeight','Bold','FontSize',18);
   title('Responder','FontWeight','Bold','FontSize',18,'FontName','Arial');
   set(gca,'FontWeight','Bold','FontSize',18,'FontName','Arial');  
  
figure
sz=25;
% map = colormap(c);
scatter(M2_mean1,C_T_mean1,sz,c,'filled');
colorbar
map = colormap(c);
   caxis([0 200])
   h=colorbar('xtick',[0,50,100,150,200])
   xlim([0 1.2]),ylim([0 1]);
   set(gca,'xtick',0:0.2:1.2); set(gca,'xticklabel',0:0.2:1.2);
   xlabel('Average M2 macrophage density','FontWeight','Bold','FontSize',18);ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18);
%    title('Positive patient','FontWeight','Bold','FontSize',18);
   title('Responder','FontWeight','Bold','FontSize',18);
   set(gca,'FontWeight','Bold','FontSize',18);    
 %Nonresponder
   figure
sz=25;
scatter(M1_mean2,C_T_mean2,sz,c,'filled');
colorbar
   caxis([0 200])
   h=colorbar('xtick',[0,50,100,150,200])
   map = colormap(c);
   xlim([0 1.2]),ylim([0 1]);
   set(gca,'xtick',0:0.2:1.2); set(gca,'xticklabel',0:0.2:1.2);
   xlabel('Average M1 macrophage density','FontWeight','Bold','FontSize',18,'FontName','Arial');
   ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18,'FontName','Arial');
%    title('Positive patient','FontWeight','Bold','FontSize',18);
   title('Nonresponder','FontWeight','Bold','FontSize',18,'FontName','Arial');
   set(gca,'FontWeight','Bold','FontSize',18,'FontName','Arial');  
  
figure
sz=25;
scatter(M2_mean2,C_T_mean2,sz,c,'filled');
colorbar
map = colormap(c);
   caxis([0 200])
   h=colorbar('xtick',[0,50,100,150,200])
   xlim([0 1.2]),ylim([0 1]);
   set(gca,'xtick',0:0.2:1.2); set(gca,'xticklabel',0:0.2:1.2);
   xlabel('Average M2 macrophage density','FontWeight','Bold','FontSize',18,'FontName','Arial');
   ylabel('Average cancer cell density','FontWeight','Bold','FontSize',18,'FontName','Arial');
%    title('Positive patient','FontWeight','Bold','FontSize',18);
   title('Nonresponder','FontWeight','Bold','FontSize',18,'FontName','Arial');
   set(gca,'FontWeight','Bold','FontSize',18,'FontName','Arial');