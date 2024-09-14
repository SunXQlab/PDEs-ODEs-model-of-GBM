% === This test program is design for estimating ODE system parameters in Signaling Pathway====
% === based on data.              ====
%
% Written by Xiaoqiang Sun, Haofeng Lin   Apr 05, 2021
% School of mathematics, Sun Yat-sen University, Guangzhou 510275, China
tic

clc
clear all
close all
figure,
global time
time=0;
global Ptest1 

group=10;%Using a genetic algorithm to estimate 10 sets of parameters.
NVAR =22;% Number of variables, according to the original setting of Para0
p=ones(group,NVAR+1);     %Let's use matrix P to store the 22 parameters for each of the 10 sets, along with their evaluation values.
for q=1:1:group
NIND = 20;           % Number of individuals per subpopulations    
MAXGEN = 400;        % maximum Number of generations     
GGAP = 0.90;         % Generation gap, how many new individuals are created
PRECI = 8;           % Precision of binary representation  
Para0= ones(1,NVAR);  % initial values in parameter space
Ptest1=Para0;
shu_min=0.000000001;           %Minimum search range    
shu_max=50;          %Maximum search range        
% Build field descriptor
FieldD = [rep(PRECI,[1, NVAR]); [shu_min*ones(1,21) 9;shu_max*ones(1,10) 1*ones(1,9) 10*ones(1,2) 10.2]
              rep([1; 0; 1 ;1], [1, NVAR])];%Decoding function.
% Initialise population
Chrom = crtbp(NIND, NVAR*PRECI);   %8-bit binary.
% Reset counters
Best = NaN*ones(MAXGEN,1);	% best in current population
gen = 1;			        % generational counter
rvChrom=bs2rv(Chrom,FieldD);  %Translation
ObjV = Objfun1(rvChrom);  %  %  %%% Objective function % Evaluate initial population
for i=1:1:MAXGEN          %
    gen=i  %To view the current iteration count in the command prompt.
     % Assign fitness-value to entire population
       FitnV = ranking(ObjV);
    % Select individuals for breeding   
       SelCh = select('sus', Chrom, FitnV, GGAP);
    % Recombine selected individuals (crossover)  
       SelCh = recombin('xovsp',SelCh,0.7);
    % Perform mutation on offspring  
       SelCh = mut(SelCh);
    % Evaluate offspring, call objective function   
       ObjVSel = Objfun1(bs2rv(SelCh,FieldD));
    % Reinsert offspring into current population  
       [Chrom ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVSel);


    % Update display and record current best individual 
       [Best(i),ind] = min(ObjV);
      %%画图
       plot(log10(Best),'ro'); xlabel('generation'); ylabel('log10(f(x))');
       text(0.5,0.95,['Best = ', num2str(Best(gen))],'Units','normalized');
       drawnow;     
        % Increment generational counter
end
chrom_value=bs2rv(Chrom,FieldD)   % Show the last set of subpopulation 
%The following is the reset of the parameter search range, with maximum and minimum values.
[value ind]=min(ObjV);
Ptest2=chrom_value(ind,:)         
for j=1:1:10
    FieldD(2,j)=max(Ptest2(j)-1,0.000000001);
    FieldD(3,j)=min(Ptest2(j)+1,100);
end
for j=11:1:19
    FieldD(2,j)=max(Ptest2(j)-0.1,0.000000001);
    FieldD(3,j)=min(Ptest2(j)+0.1,1);     
end
for j=20:1:21
    FieldD(2,j)=max(Ptest2(j)-1,0.000000001);
    FieldD(3,j)=min(Ptest2(j)+1,10);       
end
for j=22:1:22
    FieldD(2,j)=max(Ptest2(j)-2,1);
    FieldD(3,j)=min(Ptest2(j)+2,10.2);
end
b_min=FieldD(2,:);
b_max=FieldD(3,:);
%min_old+(max_old-min_old)*rand=min_new+(max_new-min_new)*rand_new
chrom_new=(Ptest2-b_min)./(b_max-b_min);
chrom_new2=chrom_new*2^8;
chrom_new2=floor(chrom_new2);
for i=0:1:7
    b=floor(chrom_new2/2^(7-i));
    for j=0:1:NVAR-1   
    Chrom(ind,((8-(7-i))+8*j))=b(j+1);
    end
    chrom_new2=rem(chrom_new2,2^(7-i));
end
rvChrom=bs2rv(Chrom,FieldD);
ObjV = Objfun1(rvChrom);  %  %  %%% Objective function % Evaluate initial population
%Below is the initiation of minor iteration after updating the search range, followed by another update of the search range.
for i=1:1:30
    gen=1;
    MAXGEN=30;
    Best = NaN*ones(MAXGEN,1);	% best in current population
    rvChrom=bs2rv(Chrom,FieldD);
    ObjV = Objfun1(rvChrom);  %  %  %%% Objective function % Evaluate initial population
while gen < MAXGEN+1
    gen%It's convenient to see where it's running to
    % Assign fitness-value to entire population
       FitnV = ranking(ObjV);
    % Select individuals for breeding
       SelCh = select('sus', Chrom, FitnV, GGAP);
    % Recombine selected individuals (crossover)
       SelCh = recombin('xovsp',SelCh,0.7);
    % Perform mutation on offspring
       SelCh = mut(SelCh);
    % Evaluate offspring, call objective function
       ObjVSel = Objfun1(bs2rv(SelCh,FieldD));

    % Reinsert offspring into current population
       [Chrom ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVSel);


    % Update display and record current best individual
       [Best(gen),ind] = min(ObjV);
      
       
       plot(log10(Best),'ro'); xlabel('generation'); ylabel('log10(f(x))');
       text(0.5,0.95,['Best = ', num2str(Best(gen))],'Units','normalized');
       drawnow;       
        % Increment generational counter
       gen = gen+1;
end 
chrom_value=bs2rv(Chrom,FieldD)   % Show the last set of subpopulation
[value,ind]=min(ObjV);
Ptest2=chrom_value(ind,:)         %Find the optimal individual.
Ptest1=Ptest2
for j=1:1:10
    FieldD(2,j)=max(Ptest2(j)-1,0.000000001);
    FieldD(3,j)=min(Ptest2(j)+1,100);
end
for j=11:1:19
    FieldD(2,j)=max(Ptest2(j)-0.1,0.000000001);
    FieldD(3,j)=min(Ptest2(j)+0.1,1);     
end
for j=20:1:21
    FieldD(2,j)=max(Ptest2(j)-1,0.000000001);
    FieldD(3,j)=min(Ptest2(j)+1,10);       
end
for j=22:1:22
    FieldD(2,j)=max(Ptest2(j)-2,1);
    FieldD(3,j)=min(Ptest2(j)+2,10.2);
end
b_min=FieldD(2,:);
b_max=FieldD(3,:);
%min_old+(max_old-min_old)*rand=min_new+(max_new-min_new)*rand_new
chrom_new=(Ptest2-b_min)./(b_max-b_min);
chrom_new2=chrom_new*2^8;
chrom_new2=floor(chrom_new2);
for k=0:1:7
    b=floor(chrom_new2/2^(7-k));
    for j=0:1:NVAR-1   %NVAR represents the number of parameters.
    Chrom(ind,((8-(7-k))+8*j))=b(j+1);
    end
    chrom_new2=rem(chrom_new2,2^(7-k));
end
 
end

 Conditions_ODES_2;
 tick={'Control group 1','Gefitinib','IGF1','IGF1+Gefitinib','EGF','EGF+Gefitinib','Control group 2','Gefitinib','5 times Gefitinib','OSI-906','5 times OSI-906'};
 b1=[147/255,224/255,255/255];b2=[153/255,77/255,82/255];
 figure,plot(1:c2,Y2(:,1),'-o','color',b1,'MarkerFaceColor',b1);hold on;plot(1:c2,pY1068_EGFR,'-o','color',b2,'MarkerFaceColor',b2);
  %set(gca,'ytick',[0.1:2]);
%  xlabel('Condition','FontWeight','Bold','FontSize',18);
 ylabel('Relative value','FontWeight','Bold','FontSize',18);
 set(gca,'XTick',1:1:11);
 set(gca,'XTickLabelRotation',45);
 legend('Model','Experiments');
 set(gca,'FontWeight','Bold','FontSize',18);
 title('pY1068-EGFR','FontWeight','Bold','FontSize',18);
 set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',14);
 figure,plot(1:c2,Y2(:,2),'-o','color',b1,'MarkerFaceColor',b1);hold on;plot(1:c2,p_ERK,'-o','color',b2,'MarkerFaceColor',b2);
 set(gca,'XTick',1:1:11);
 set(gca,'XTickLabelRotation',45);
 set(gca,'YTick',0.4:0.1:1); 
 set(gca,'Yticklabel',0.4:0.1:1,'FontWeight','Bold','FontSize',18);
 ylabel('Relative value','FontWeight','Bold','FontSize',18);
 legend('Model','Experiments');
 set(gca,'FontWeight','Bold','FontSize',18);
 title('p-ERK','FontWeight','Bold','FontSize',18);
 set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',14);
 figure,plot(1:c2,Y2(:,3),'-o','color',b1,'MarkerFaceColor',b1);hold on;plot(1:c2,p_IGF1Rbeta,'-o','color',b2,'MarkerFaceColor',b2);
 set(gca,'XTick',1:1:11);
 set(gca,'XTickLabelRotation',45);
%  xlabel('Condition','FontWeight','Bold','FontSize',18);
 ylabel('Relative value','FontWeight','Bold','FontSize',18);
 legend('Model','Experiments');
 set(gca,'FontWeight','Bold','FontSize',18);
 title('p-IGF1Rbeta','FontWeight','Bold','FontSize',18);
 set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',14);
 figure,plot(1:c2,Y2(:,4),'-o','color',b1,'MarkerFaceColor',b1);hold on;plot(1:c2,pS473_AKT,'-o','color',b2,'MarkerFaceColor',b2); 
  set(gca,'XTick',1:1:11);
 set(gca,'XTickLabelRotation',45);
%  xlabel('Condition','FontWeight','Bold','FontSize',18);
 ylabel('Relative value','FontWeight','Bold','FontSize',18);
 legend('Model','Experiments');
 set(gca,'FontWeight','Bold','FontSize',18);
 title('pS473-AKT','FontWeight','Bold','FontSize',18);
  set(gca,'Xticklabel',tick,'FontWeight','Bold','FontSize',14);
 figure
 p(q,1:NVAR)=Ptest2;
 p(q,NVAR+1)=value;
end
save('canshu1-24hours-ode15s_wode4.mat','p');


 

