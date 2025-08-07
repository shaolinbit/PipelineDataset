clear all;
close all;
clc
ts=0.005;
dir=fileparts(mfilename('fullpath'));
filePath = fullfile(dir, 'Sim_result3_2', 'my_motion_def_tube3_2', 'ref_pva_part.txt');
ref_pva = dlmread(filePath, ' ', 1, 0);

ref_enu=lla2enu(ref_pva(:,2:4),ref_pva(1,2:4));

filePathmy_pva = fullfile(dir, 'Sim_result3_2', 'opt', 'opt_result.txt');
my_pva=dlmread(filePathmy_pva);
my_enu=lla2enu(my_pva(:,2:4),ref_pva(1,2:4));


filePathEKF_pva = fullfile(dir, 'Sim_result3_2', 'EKF', 'EKF_results_smooth.txt');
EKF_pva=dlmread(filePathEKF_pva);
EKF_enu=lla2enu(EKF_pva(:,8:10),ref_pva(1,2:4));


row=0;
for i=1:size(my_pva,1)
    [index,~]=find(abs(ref_pva(row+1:end,1)-my_pva(i,1))<(0.5*ts));
    if isempty(index)
        index=0;
    end
    row=index+row;
    time_diff=ref_pva(row,1)-my_pva(i,1);
    
    my_enu_error(i,:)=my_enu(i,:)-ref_enu(row,:);
end

row=0;
for i=1:size(EKF_pva,1)
    [index,~]=find(abs(ref_pva(row+1:end,1)-EKF_pva(i,1))<(0.5*ts));
    if isempty(index)
        index=0;
    end
    row=index+row;
    time_diff=ref_pva(row,1)-EKF_pva(i,1);
    
    EKF_enu_error(i,:)=EKF_enu(i,:)-ref_enu(row,:);
end

ylabels={'Longitude error(m)','Latitude error(m)','Height error(m)'};
labelVarwin1 = 'Extended Kalman Filter';
labelVarwin2 = 'Factor Graph';

for i=1:3
    figure;
    plot(my_pva(:,1),my_enu_error(:,i),'LineWidth',2);hold on;
    xlim([300,801]);
    plot(EKF_pva(:,1),EKF_enu_error(:,i),'LineWidth',2);hold on;
    xlim([300,801]);
    set(gca, 'FontSize', 24);
    xlabel('Time(s)');ylabel(ylabels{i});
    box off;
    %grid on;
    legend(labelVarwin1,labelVarwin2','FontSize', 20,'Location', 'north');%
end


function enu=lla2enu(lla,lla0)
ecef = lla2ecef(lla, 'WGS84');
wgs84 = wgs84Ellipsoid('meter');
[enu(:,1),enu(:,2),enu(:,3)] = ecef2enu(ecef(:,1),ecef(:,2),ecef(:,3),lla0(1),lla0(2),lla0(3),wgs84 );
end