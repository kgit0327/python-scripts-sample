iData = 1;
iPeriod = 2;
% [JT, JAV, JTP] = GetJATP(ExperimentDatas(iData));
JT.sh = ExperimentDatas(iData).JT(4:6, :);
JAV.sh = ExperimentDatas(iData).JAV(4:6, :);
JTP.sh = ExperimentDatas(iData).JTP(4:6, :);
[GC, MinHB, MaxER, RM] = GetTiming(ExperimentDatas, GC2RMFrame, iData, iPeriod);
l = RM-GC;


figure;
subplot(311)
plot((-l-20)/200:1/200:20/200, JT.sh(1, GC-20:RM+20), '-', 'LineWidth', 3, 'Color', 'r')
hold on;
plot((-l-20)/200:1/200:20/200, JT.sh(2, GC-20:RM+20), '-', 'LineWidth', 3, 'Color', 'g')
plot((-l-20)/200:1/200:20/200, JT.sh(3, GC-20:RM+20), '-', 'LineWidth', 3, 'Color', 'b')
% plot((-l-20)/200:1/200:20/200, segdat(14).E(1, GC-20:RM+20), '-', 'LineWidth', 3, 'Color', 'k')

xline(-l/200)
xline((GC-MinHB)/200)
xline((GC-MaxER)/200)
xline((GC-GC)/200)

xlim([(-l-20)/200 20/200])
set(gca, 'FontSize', 20)
ylabel('Joint Torque [N m]', 'FontSize', 20)
legend('水平屈曲', '内転', '内旋')

subplot(312)
plot((-l-20)/200:1/200:20/200, JAV.sh(1, GC-20:RM+20), '-', 'LineWidth', 3, 'Color', 'r')
hold on;
plot((-l-20)/200:1/200:20/200, -JAV.sh(2, GC-20:RM+20), '-', 'LineWidth', 3, 'Color', 'g')
plot((-l-20)/200:1/200:20/200, JAV.sh(3, GC-20:RM+20), '-', 'LineWidth', 3, 'Color', 'b')
% plot((-l-20)/200:1/200:20/200, segdat(14).E(1, GC-20:RM+20), '-', 'LineWidth', 3, 'Color', 'k')

xline(-l/200)
xline((GC-MinHB)/200)
xline((GC-MaxER)/200)
xline((GC-GC)/200)

xlim([(-l-20)/200 20/200])
set(gca, 'FontSize', 20)
ylabel('Joint Angular Velocity [deg/s]', 'FontSize', 20)
legend('水平屈曲', '内転', '内旋')

JTP.sh(2, :) = JT.sh(2, :).*-JAV.sh(2, :);

subplot(313)
plot((-l-20)/200:1/200:20/200, JTP.sh(1, GC-20:RM+20), '-', 'LineWidth', 3, 'Color', 'r')
hold on;
plot((-l-20)/200:1/200:20/200, JTP.sh(2, GC-20:RM+20), '-', 'LineWidth', 3, 'Color', 'g')
plot((-l-20)/200:1/200:20/200, JTP.sh(3, GC-20:RM+20), '-', 'LineWidth', 3, 'Color', 'b')
plot((-l-20)/200:1/200:20/200, JTP.sh(1, GC-20:RM+20)+JTP.sh(2, GC-20:RM+20)+JTP.sh(3, GC-20:RM+20), '-', 'LineWidth', 3, 'Color', 'k')

xline(-l/200)
xline((GC-MinHB)/200)
xline((GC-MaxER)/200)
xline((GC-GC)/200)

xlim([(-l-20)/200 20/200])
% legend('hand', 'forearm', 'upper arm')
set(gca, 'FontSize', 20)
xlabel('Time [s]', 'FontSize', 20)
ylabel('Joint Torque Power [W]', 'FontSize', 20)
legend('水平屈曲', '内転', '内旋', '合計')

samexaxis('XMinorTick','on','Join','YTickAntiClash','YLabelDistance',1.0)

function [JT, JAV, JTP] = GetJATP(LoadedData)
    % 関節トルク
    JT.sh = LoadedData.FN(1,1).N; % 肩関節周りの上腕トルク uaARsh
    JT.el = LoadedData.FN(1,2).N; % faARel
    JT.wr = LoadedData.FN(1,3).N; % hdARwr
    
    JT.uaARel = -JT.el;
    JT.faARwr = -JT.wr;


    
    % 関節角速度
    JAV.sh = LoadedData.segdat(3).seganV_GCS - LoadedData.segdat(14).seganV_GCS;
    JAV.el = LoadedData.segdat(2).seganV_GCS - LoadedData.segdat(3).seganV_GCS;
    JAV.wr = LoadedData.segdat(1).seganV_GCS - LoadedData.segdat(2).seganV_GCS;

    
    
    % 関節トルクパワー
    JTP.sh = JT.sh.*JAV.sh; % sh2ua
    JTP.el = JT.el.*JAV.el; % el2fa
    JTP.wr = JT.wr.*JAV.wr; % wr2hd
    
    JTP.el2ua = JT.el.*JAV.el;
    JTP.wr2fa = JT.wr.*JAV.wr;
    
  

end

function [GC, MinHB, MaxER, RM] = GetTiming(ExperimentDatas, GC2RMFrame, iData, iPeriod)
    GC = GC2RMFrame(iData).GC2RM(1, iPeriod);

    top = ExperimentDatas(iData).n.Top;
    top_hv = dif3(top(1:2, :), size(top, 2), 1/200);
    top_hv_r = sqrt(sum(top_hv.*top_hv));

    [~, MinHB] = min(top_hv_r(1, GC2RMFrame(iData).GC2RM(1, iPeriod):GC2RMFrame(iData).GC2RM(2, iPeriod)), [], 2);
    MinHB = MinHB+GC2RMFrame(iData).GC2RM(1, iPeriod)-1;

    sh_ex = ExperimentDatas(iData).Unit(1).JA(3, :);

    [~, MaxER] = max(sh_ex(1, GC2RMFrame(iData).GC2RM(1, iPeriod):GC2RMFrame(iData).GC2RM(2, iPeriod)));
    MaxER = MaxER+GC2RMFrame(iData).GC2RM(1, iPeriod)-1;

    RM = GC2RMFrame(iData).GC2RM(2, iPeriod);
end