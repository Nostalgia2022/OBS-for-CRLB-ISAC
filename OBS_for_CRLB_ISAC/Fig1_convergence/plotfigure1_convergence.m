clc
clear
load data_convergence



Z1=Con{1}(:,3);
Z2=Con{2}(:,3);
Z3=Con{3}(:,3);


figure(1)
xlabel('Number of Iterations')
ylabel('Objective value')
grid on
hold on
plot(Z1,'Color',color(1),'LineWidth',1.5,'LineStyle','-')
plot(Z2,'Color',color(2),'LineWidth',1.5,'LineStyle','-')
plot(Z3,'Color',color(3),'LineWidth',1.5,'LineStyle','-')

legend({'\delta_s=1, \delta_c=0.05','\delta_s=1, \delta_c=0.1','\delta_s=1, \delta_c=0.15'},'Location','southeast')
% set(gcf,'Units','centimeters','Position',[10 5 7 5])
set(gca,'FontName','Times New Rome','FontSize',10)

% figure(2)
% xlabel('Number of Iterations')
% ylabel('Objective value')
% grid on
% hold on
% plot(Z2,'Color',color(2),'LineWidth',1.5,'LineStyle','-')

% legend({'SCA-SGPI'},'Location','southeast')
% set(gcf,'Units','centimeters','Position',[10 5 7 5])
% set(gca,'FontName','Times New Rome','FontSize',10)