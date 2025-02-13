clc
clear
load data

K1=mean(squeeze(Obj_all(2,:,:)));
K2=mean(squeeze(Obj_all(3,:,:)));
K3=mean(squeeze(Obj_all(4,:,:)));
index=0:8;
index2=1:8;
load data_sensing.mat
K0=stream_test;

load data_SDR_sensing.mat
K0_SDR=-CRB_all;

load data_SDR_Ns.mat

obj_all_new=zeros(3,100,9);
for k_par=1:300
channel=mod((k_par-1),I_out)+1;
weight=floor((k_par-1)/I_out)+1;
obj_all_new(weight,channel,:)=0.25*SR_all(:,k_par)-CRB_all(:,k_par);
end

K1_SDR=mean(squeeze(obj_all_new(1,:,:)));
K2_SDR=mean(squeeze(obj_all_new(2,:,:)));
K3_SDR=mean(squeeze(obj_all_new(3,:,1)))*ones(1,9);



figure
xlabel('Number of sensing streams')
ylabel('Objective value')
grid on
hold on

plot(index,K1,'Color',color(1),'LineWidth',1.5,'LineStyle','-','Marker','s')
plot(index,K2,'Color',color(2),'LineWidth',1.5,'LineStyle','-','Marker','o')
plot(index,K3,'Color',color(3),'LineWidth',1.5,'LineStyle','-','Marker','v')
plot(index2,K0,'Color',color(4),'LineWidth',1.5,'LineStyle','-','Marker','d')

legend({'K=1','K=2','K=3','K=0'},'Location','southeast')
% set(gcf,'Units','centimeters','Position',[10 5 7 5])
set(gca,'FontName','Times New Rome','FontSize',10)


plot(index,K1_SDR,'Color',color(11),'LineWidth',1.5,'LineStyle','--','Marker','s')
plot(index,K2_SDR,'Color',color(12),'LineWidth',1.5,'LineStyle','--','Marker','o')
plot(index,K3_SDR,'Color',color(13),'LineWidth',1.5,'LineStyle','--','Marker','v')
plot(index2,K0_SDR,'Color',color(14),'LineWidth',1.5,'LineStyle','--','Marker','d')

legend({'Algorithm 1, K=1','Algorithm 1, K=2','Algorithm 1, K=3','Algorithm 1, K=0','WMMSE-SDR, K=1','WMMSE-SDR, K=2','WMMSE-SDR, K=3','WMMSE-SDR, K=0'},'Location','southeast')
xlim([0, 8]);
ylim([-12, 0]);


% figure(2)
% xlabel('Number of Iterations')
% ylabel('Objective value')
% grid on
% hold on
% plot(Z2,'Color',color(2),'LineWidth',1.5,'LineStyle','-')

% legend({'SCA-SGPI'},'Location','southeast')
% set(gcf,'Units','centimeters','Position',[10 5 7 5])
% set(gca,'FontName','Times New Rome','FontSize',10)