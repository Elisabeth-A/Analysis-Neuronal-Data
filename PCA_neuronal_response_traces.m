%% PCA of dendritic data
clear
clc
% 1. make array with average traces of all ROIs (neurons), combining the
% response traces for all time points (hab and ret) and stimuli (CS- and CS+)
%   in the format time bins x neurons
% (2. mean substract --> not needed for matlab PCA function, uses singular
%     value decomposition to centre the data)
% 3. do PCA
% 4. project averages response traces for the different stimuli and time
% points into the space of the first 3 principal components


%% load data to get all df-f traces for all subjects (mice)
mouse_ID={}

% Matrix combining all responses across sessions and stimuli for all ROIs
% across all subjects
M=[];

% mean traces of all ROIs per stimulus and session 
minus_hab=[];
minus_ret=[];

plus_hab=[];
plus_ret=[];


% load files containing response traces and combine data from all ROIs and subjects 
for i=1:length(mouse_ID)
   
    filename=char(strcat(mouse_ID(i), '_df_f_bin.mat'))    
    load(filename, 'df_f_bin', 'time_bins')
    
    for j=1:length(df_f_bin{1,1})
        M=[M nanmean([df_f_bin{1,2}{1,j}, (df_f_bin{1,6}{1,j}), df_f_bin{1,1}{1,j}, (df_f_bin{1,5}{1,j})],2)];
        minus_hab=[minus_hab nanmean(df_f_bin{1,2}{1,j},2)];
        minus_ret=[minus_ret nanmean(df_f_bin{1,6}{1,j},2)];
        plus_hab=[plus_hab nanmean(df_f_bin{1,1}{1,j},2)];
        plus_ret=[plus_ret nanmean(df_f_bin{1,5}{1,j},2)];
        
    end
      
end

time=time_bins;
%% PCA
% MATLAB function includes centering of data

[coeff, score, latent, tsquared, explained, mu]=pca(M);



%% centre test data by substracting column means of full matrix 
minus_hab_c=bsxfun(@minus,minus_hab, nanmean(M));
minus_ret_c=bsxfun(@minus,minus_ret, nanmean(M));
plus_hab_c=bsxfun(@minus,plus_hab, nanmean(M));
plus_ret_c=bsxfun(@minus,plus_ret, nanmean(M));


%% project data onto first 3 principal component
minus_hab_project= [minus_hab_c* coeff(:,1),minus_hab_c* coeff(:,2),minus_hab_c* coeff(:,3)]
minus_ret_project= [minus_ret_c* coeff(:,1),minus_ret_c* coeff(:,2),minus_ret_c* coeff(:,3)]
plus_hab_project= [plus_hab_c* coeff(:,1),plus_hab_c* coeff(:,2),plus_hab_c* coeff(:,3)]
plus_ret_project= [plus_ret_c* coeff(:,1),plus_ret_c* coeff(:,2),plus_ret_c* coeff(:,3)]


%% plot data in 3D plot comparing sessions (habituation & 24h recall)

subplot(1,2,1)
p1=plot3(minus_hab_project(:,1),minus_hab_project(:,2),minus_hab_project(:,3), 'LineWidth', 2)
hold on
p2=plot3(minus_ret_project(:,1),minus_ret_project(:,2),minus_ret_project(:,3),'LineWidth', 2)
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
title('CS-')
hold off
%ylim([-0.3 0.4])
%xlim([-0.4 0.8])
%zlim([-0.2 0.3])
grid on

subplot(1,2,2)
p3=plot3(plus_hab_project(:,1),plus_hab_project(:,2),plus_hab_project(:,3),'LineWidth', 2)
hold on
p4=plot3(plus_ret_project(:,1),plus_ret_project(:,2),plus_ret_project(:,3),'LineWidth', 2)
legend('habituation', '24h recall', 'Location','northwestoutside')
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
title('CS+')
%ylim([-0.3 0.4])
%xlim([-0.4 0.8])
%zlim([-0.2 0.3])
grid on
propertyeditor on


%% plot data in 3D plot comparing stimuli (CS- vs CS+)
figure
subplot(1,2,1)
p1=plot3(minus_hab_project(:,1),minus_hab_project(:,2),minus_hab_project(:,3), 'LineWidth', 2)
hold on
p2=plot3(plus_hab_project(:,1),plus_hab_project(:,2),plus_hab_project(:,3),'LineWidth', 2)
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
title('habituation')
% ylim([-1.5 2.5])
% xlim([-2 4])
% zlim([-0.8 0.5])
hold off
grid on

subplot(1,2,2)
p3=plot3(minus_ret_project(:,1),minus_ret_project(:,2),minus_ret_project(:,3),'LineWidth', 2)
hold on
p4=plot3(plus_ret_project(:,1),plus_ret_project(:,2),plus_ret_project(:,3),'LineWidth', 2)
legend('CS-', 'CS+', 'Location','northwestoutside')
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
title('24h recall')
% ylim([-1.4 2.5])
% xlim([-2 4])
% zlim([-0.8 0.5])
grid on
propertyeditor on


%% calculate euclidean distance between response trajectories and plot over time

% comparing stimuli distance for each session
d_CS_hab= sqrt(sum((minus_hab_project-plus_hab_project).^2,2))
d_CS_ret= sqrt(sum((minus_ret_project-plus_ret_project).^2,2))

figure
plot(time,d_CS_hab,'g', 'LineWidth',2)
hold on
plot(time, d_CS_ret,'r', 'LineWidth',2)
legend('habituation', '24h recall')
vline([20 25]) % mark stimulus interval
xlabel('time')
ylabel('distance')
title('distance traj. stimuli')

% comparing sessions distance for each stimulus
d_minus= sqrt(sum((minus_hab_project-minus_ret_project).^2,2))
d_plus= sqrt(sum((plus_hab_project-plus_ret_project).^2,2))

figure
plot(time,d_minus, 'LineWidth',2)
hold on
plot(time, d_plus, 'LineWidth',2)
legend('CS-', 'CS+')
vline([20 25]) % mark stimulus interval
xlabel('time')
ylabel('distance')
title('distance traj. session')
