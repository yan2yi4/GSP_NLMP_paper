clc;clear;
% close all;
load('Ext_Real_graph_and_sampling.mat')

runs_total = 100;

rng(12345);%for reproducibility 

L = diag(sum(A)) - A;
[n_nodes,numberSamples] = size(completeDatasetMatrix);


gs = completeDatasetMatrix;

alpha = 1.2;
gam = 0.1;

p = alpha-0.05;

names = {'GLMS','GLMP','GNLMS','GNLMP'};
u_lms = 0.8;
u_nlms = 0.8;
u_lmp = 1.6;
u_nlmp = 0.55;

alg_selection_vec = [1 2 3 4]; 
alg_param_vec = [u_lms u_lmp u_nlms u_nlmp]; 

[~,num_algs] = size(alg_param_vec);
max_it = numberSamples;

sigma_diag = zeros(n_nodes,1);
sigma_diag(used_indices) = 1;

mean_MSD_mat_comp = [];

mean_elapsedTime_mat_comp = [];

mean_MSD_mat_comp_tv = [];
threshold = FLOM( p-1, alpha,gam)*M;
for i = 1:num_algs 
    [mean_MSD_mat_comp(:,i), ...
        mean_x_vec_matrix(:,:,i), mean_elapsedTime_mat_comp(:,i),mean_MSD_mat_comp_tv(:,i)] = ...
      Ext_run_loop_NLMP( gs, ... 
            D_s, U_f, [alpha gam], ...
            alg_param_vec(i), alg_selection_vec(i), ...
            max_it, runs_total, p,0,sigma_diag,U,threshold) ;
end

    % ------------------------------------------------------------------- %
%%
 
sampled_position = find(diag(D_s)==1); 

figure
x_values = 1:numberSamples;
pos = 80;
hold on
for i = 1:num_algs 
    y_values = mean_x_vec_matrix(sampled_position(pos),:,i); 
    plot( x_values, y_values , 'LineWidth', 2 )
end
plot( x_values, gs(sampled_position(pos),:) ,'black.', 'LineWidth', 2)
leg = legend(names,'Location','southeast','NumColumns',2);
hold off
grid on
xlabel('Iteration, k')
ylabel('Amplitude')

figure
% left = 1.7*k_scaling/2;	% normalized left margin
% set(0,'defaultAxesPosition',[left/width bottom/hight (width-left-right)/width  (hight-bottom-top)/hight]);
hold on
for i = 1:num_algs
    y_values = mean_MSD_mat_comp_tv(x_values,i);
    plot( x_values, 10*log10( y_values ),'LineWidth', 2);
end
leg = legend(names,'NumColumns',2);
set(gca,'gridlinestyle','-');
grid on
hold off
xlabel('Iteration, k')
ylabel('NMSD_t(dB)')

% figure
% left = 1.5*k_scaling/2;	% normalized left margin
% set(0,'defaultAxesPosition',[left/width bottom/hight (width-left-right)/width  (hight-bottom-top)/hight]);
% x_values = 1:numberSamples;
% pos = 1;
% hold on
% for i = 1:num_algs 
%     sampleRate = 1;
%     y_values = zeros(1,numberSamples);
%     if i == 2 || i == 4
%         for j = 1:numberSamples
%             if j ==1
%                 y_values(j) =mean_elapsedTime_mat_comp(j,i); 
%             else
%                 y_values(j) = y_values(j-1)+mean_elapsedTime_mat_comp(j,i);
%             end
%         end
%         plot( x_values, y_values , 'LineWidth', 2 )
%     end
% end
% leg = legend('LMP','NLMP','Location','southeast','NumColumns',2);
% hold off
% xlabel('Iteration, k')
% ylabel('Time(s)')
% grid on

%%
sum(mean_elapsedTime_mat_comp)