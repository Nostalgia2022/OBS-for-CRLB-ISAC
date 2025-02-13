
clc
clear

% Load data
data1 = load('data_WMMSESDR_pt.mat');
data2 = load('data_FP_SGDA.mat');
data3 = load('data_SCA_Ns=3M.mat');

% Extract and average SR_all and CRB_all
SR1 = mean(data1.SR_all, 1);  
SR2 = mean(data2.SR_all, 1);
SR3 = mean(data3.SR_all, 1);

CRB1 = mean(data1.CRB_all, 1);
CRB2 = mean(data2.CRB_all, 1);
CRB3 = mean(data3.CRB_all, 1);

% Compute the weighted sum
WS1 = 0.25 * SR1 - CRB1;
WS2 = 0.25 * SR2 - CRB2;
WS3 = 0.25 * SR3 - CRB3;

% X-axis labels (assuming indices from 1 to 6 correspond to methods or parameters)
x_labels = 1:6;  

% Create figure
figure; 
hold on;
grid on;

% Plot SR
line(1)=plot(x_labels, SR1, '-.o', 'LineWidth', 1.5,'Color', color(1), 'DisplayName', 'WMMSE-SDR, SR');
line(2)=plot(x_labels, SR2, '-.s', 'LineWidth', 1.5, 'Color', color(2),'DisplayName', 'FP-SGDA, SR');
line(3)=plot(x_labels, SR3, '-.^', 'LineWidth', 1.5, 'Color', color(3),'DisplayName', 'Algorithm 1, SR');

% Plot CRB
line_CRB(1)=plot(x_labels, CRB1, '--<', 'LineWidth', 1.5,'Color', color(1), 'DisplayName', 'WMMSE-SDR, CRLB');
line_CRB(2)=plot(x_labels, CRB2, '-->', 'LineWidth', 1.5,'Color', color(2), 'DisplayName', 'FP-SGDA, CRLB');
line_CRB(3)=plot(x_labels, CRB3, '--v', 'LineWidth', 1.5,'Color', color(3), 'DisplayName', 'Algorithm 1, CRLB');

% Plot Weighted Sum (WS)
line_OV(1)=plot(x_labels, WS1, '-d', 'LineWidth', 1.5, 'Color', color(1), 'DisplayName', 'WMMSE-SDR, OV');
line_OV(2)=plot(x_labels, WS2, '-x', 'LineWidth', 1.5, 'Color', color(2), 'DisplayName', 'FP-SGDA, OV');
line_OV(3)=plot(x_labels, WS3, '-p', 'LineWidth', 1.5, 'Color', color(3), 'DisplayName', 'Algorithm 1, OV');

% Labeling
xlabel('Transmit Power (dBm)');
ylabel('Metric Values');

lgd1 = legend(line, 'Location', 'NorthWest');
title(lgd1, 'Sum Rate');
axesNew1=axes('position',get(gca,'position'),'visible','off');
lgd2 = legend(axesNew1,line_CRB, 'Location', 'NorthEast');
title(lgd2, 'CRLB');
axesNew=axes('position',get(gca,'position'),'visible','off');
lgd3 = legend(axesNew,line_OV, 'Location', 'SouthEast');
title(lgd3, 'Objective Value');

