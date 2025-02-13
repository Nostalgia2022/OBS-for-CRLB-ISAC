clc
clear

load data_WMMSESDR_to4

CRB_WMMSE_SDR_10dB=mean(CRB_all);
SR_WMMSE_SDR_10dB=mean(SR_all);
load data_proposed_SCA

CRB_PSLA_10dB=mean(CRB_all);
SR_PSLA_10dB=mean(SR_all);

load data_FP_SGDA.mat
CRB_FP_10dB=mean(CRB_all);
SR_FP_10dB=mean(SR_all);

figure
hold on
grid on
% slg1=plot(CRB_SCA_SGPI_20dB,SR_SCA_SGPI_20dB,'-v',CRB_SCA_SDR_20dB,SR_SCA_SDR_20dB,'--s');
% slg1(1).LineWidth=1.5;
%  slg1(1).MarkerIndices=[1,4:4:188,191];
% slg1(1).Color=color(1);
% slg1(2).LineWidth=1.5;
% slg1(2).MarkerIndices=[1,8:8:80,80:4:160,161];
% slg1(2).Color=color(2);

xlabel('Trace of the Inverse of the FIM')
ylabel('Sum Rate (nats/Hz)')


slg1=plot(CRB_WMMSE_SDR_10dB,SR_WMMSE_SDR_10dB,'-v',CRB_FP_10dB,SR_FP_10dB,'-d',CRB_PSLA_10dB,SR_PSLA_10dB,'--s');
slg1(1).LineWidth=1.5;
slg1(1).MarkerIndices=[1:2:60];
slg1(1).Color=color(1);
slg1(2).LineWidth=1.5;
slg1(2).MarkerIndices=[1:2:60];
slg1(2).Color=color(2);
slg1(3).LineWidth=1.5;
slg1(3).MarkerIndices=[1:2:60];
slg1(3).Color=color(3);


legend({'WMMSE-SDR','FP-SGDA','Algorithm 1'},'Location','southeast')