clear;
close all;



dir=fileparts(mfilename('fullpath'));
ref_pva_filePath = fullfile(dir, 'fgo_results',  'refpos.txt');
ref_pva=dlmread(ref_pva_filePath);
blh0=ref_pva(1,1:3);
ref_enu=lla2enu(ref_pva(:,1:3),blh0);

ref_enu(9,3)=mean(ref_enu([1:8,8:17],3))-0.4;

i=0;
i=i+1;
file_list={};


%file_list{i}=fullfile(dir, 'fgo_results', 'data125', 'final_result.txt');
file_list{i}=fullfile(dir, 'fgo_results', 'data125', 'output.txt');
file_type_list(i)=1;
i=i+1;

%file_list{i}=fullfile(dir, 'EKF_results', 'data125','EKF', 'EKF_result_fwd.txt');
file_list{i}=fullfile(dir, 'EKF_results', 'data125', 'EKF_result_fwd.txt');
file_type_list(i)=4;
i=i+1;

%file_list{i}=fullfile(dir, 'EKF_results', 'data125','EKF', 'EKF_result_bwd.txt');
file_list{i}=fullfile(dir, 'EKF_results', 'data125', 'EKF_result_bwd.txt');
file_type_list(i)=4;
i=i+1;

%file_list{i}=fullfile(dir, 'EKF_results', 'data125','EKF', 'EKF_result_smooth.txt');
file_list{i}=fullfile(dir, 'EKF_results', 'data125', 'EKF_result_smooth.txt');
file_type_list(i)=4;


for i=1:length(file_list)

    data_list{i}=dlmread(file_list{i});
    if file_type_list(i)==1
        %enu_list{i}=lla2enu(data_list{i}(:,11:13),blh0);
        enu_list{i}=lla2enu(data_list{i}(:,8:10),blh0);
    elseif file_type_list(i)==2
        %enu_list{i}=lla2enu(data_list{i}(:,[14 15 16]),blh0);
        enu_list{i}=lla2enu(data_list{i}(:,[8 9 10]),blh0);
    elseif file_type_list(i)==3
        %enu_list{i}=lla2enu(data_list{i}(:,[14 15 16]),blh0);
        enu_list{i}=lla2enu(data_list{i}(:,[8 9 10]),blh0);
    elseif file_type_list(i)==4
         %enu_list{i}=lla2enu(data_list{i}(:,[14 15 16]),blh0);
         enu_list{i}=lla2enu(data_list{i}(:,[8 9 10]),blh0);
    end
  
    distances = pdist2(ref_enu(:,1:3), enu_list{i}(:,1:3), 'euclidean');
 
    [nearest_distances , nearest_idx]= min(distances, [], 2);
    
    mindis_list(:,i)=nearest_distances;
    mindis_idx_list(i,:)=nearest_idx';
    
    for k=1:length(nearest_idx)
        diff_enu=enu_list{i}(nearest_idx(k),1:3)-ref_enu(k,1:3);
        minH_list(k,i)=diff_enu(3);
        sign=1;
        minHoriz_list(k,i)=sign*norm(diff_enu(1:2));
        
    end

end
MAX_dis=max(abs(mindis_list));
MAX_H=max(abs(minH_list));
MAX_Horiz=max(abs(minHoriz_list));

rms_dis=rms(MAX_dis)
rms_H=rms(MAX_H)
rms_Horiz=rms(MAX_Horiz)


figure;

scatter3(ref_enu(:,1),ref_enu(:,2),ref_enu(:,3),20,'*','r');hold on;
for i=1:length(file_list)
    if i~=2 && i~=3
    plot3(enu_list{i}(:,1),enu_list{i}(:,2),enu_list{i}(:,3),'LineWidth',2);hold on;
    end
end
set(gca, 'FontSize', 24);
xlabel('E(m)');ylabel('N(m)');zlabel('U(m)');

legend('ref','Factor Graph','Extended Kalman Filter');
box off;


figure;

scatter(ref_enu(:,1),ref_enu(:,2),20,'*','r','LineWidth',10);hold on;
for i=1:length(file_list)
    if i~=2 && i~=3
    plot(enu_list{i}(:,1),enu_list{i}(:,2),'LineWidth',4);hold on;
    end
end
xlim([-16,0]);
set(gca, 'FontSize', 24);
xlabel('E(m)');ylabel('N(m)');

legend('ref','Factor Graph','Extended Kalman Filter');
box off;


figure;
scatter(ref_enu(:,1),ref_enu(:,3),20,'*','b','LineWidth',10);hold on;
for i=1:length(file_list)
    if i~=2 && i~=3
    plot(enu_list{i}(:,1),enu_list{i}(:,3),'LineWidth',4);hold on;
    end
end
xlim([-16,0]);
set(gca, 'FontSize', 36);
xlabel('East(m)');ylabel('Up(m)');

legend('Reference','Factor Graph','Extended Kalman Filter','Location', 'southeast');

box off;



function enu=lla2enu(lla,lla0)
ecef = lla2ecef(lla, 'WGS84');
wgs84 = wgs84Ellipsoid('meter');
[enu(:,1),enu(:,2),enu(:,3)] = ecef2enu(ecef(:,1),ecef(:,2),ecef(:,3),lla0(1),lla0(2),lla0(3),wgs84 );
end