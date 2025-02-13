


% title('average time vs. SNR')

load("data_WMMSESDR_Nt.mat")
T1=Time_all;
Y1=Obj_all;
load("data_FP_SGDA.mat")
T2=Time_all;
Y2=Obj_all;
load("data_SCA_Ns=3M.mat")
T3=Time_all;
Y3=Obj_all;
load("data_SCA_Ns=0.mat")
T4=Time_all;
Y4=Obj_all;
load("data_SCA_LD_Ns=3M.mat")
T5=Time_all;
Y5=Obj_all;

x=[8,16,32];
x2=[8,16,32,64,128];

y1=mean(T1);
y2=mean(T2);
y3=mean(T3);
y4=mean(T4);
y5=mean(T5);
% 
figure
slg=loglog(x,y1,'-s',x2,y2,'-*',x2,y3,'-v',x2,y4,'-^',x2,y5,'--d');
slg(1).LineWidth=1.5;
slg(2).LineWidth=1.5;
slg(3).LineWidth=1.5;
slg(4).LineWidth=1.5;
slg(5).LineWidth=1.5;

slg(1).Color=color(1);
slg(2).Color=color(2);
slg(3).Color=color(3);
slg(4).Color=color(4);
slg(5).Color=color(5);

xlabel('Number of antennas')
ylabel('Average CPU Time (sec)')

xlim([8,128]);
xticks(x2);
xticklabels(x2)
% ylim([1e-3,1e3]);
% yticks([1e-3,1e-2,1e-1,1,1e1,1e2,1e3]);
% yticklabels({'10^{-3}','10^{-2}','10^{-1}','10^{0}','10^1','10^2','10^3'})
grid on

legend({'WMMSE-SDR','FP-SGDA','Algorithm 1, Ns=3M','Algorithm 1, Ns=0','LD Algorithm1, Ns=3M'},...
    'Location','northwest')
% axesNew=axes('position',get(gca,'position'),'visible','off');
% 
% legend(axesNew,[slg(5),slg(6),slg(7),slg(8)],{'RFP-MFPI-s','FP','WMMSE','SCA'},'Location','northeast')
