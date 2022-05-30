clc;clear;
% close all;
load('Ext_Real_graph_and_sampling.mat')

ensemble = 100;

rng(12345);%for reproducibility 

L = diag(sum(A)) - A;
[n_nodes,numberSamples] = size(completeDatasetMatrix);
gs = cat(3, completeDatasetMatrix, completeDatasetMatrix2);

alpha = 1.2;
gam = 0.1;

p = alpha-0.05;

names = {'LMS','LMP','NLMS','NLMP'};
u_nlmp_1 = 0.55;
u_nlmp_2 = 0.475;

u_mat = diag([u_nlmp_1,u_nlmp_2]);
[~,features] = size(u_mat);

max_it = numberSamples;

sigma_diag = zeros(n_nodes,1);
sigma_diag(used_indices) = 1;

mean_MSD_mat_comp = [];

mean_elapsedTime_mat_comp = [];

mean_MSD_mat_comp_tv = [];
threshold = FLOM( p-1, alpha,gam)*M;
%%
num_algs = 1;
for i = 1:num_algs 
    [~, ...
        mean_x_vec_matrix(:,:,:,i), mean_elapsedTime_mat_comp(:,i),~] = ...
      Ext_run_loop_NLMP_mult( gs, ... 
            D_s, U_f, [alpha gam], ...
            u_mat, ...
            max_it, ensemble, p,0,sigma_diag,U,threshold) ;
end

%%
    
sampled_position = find(diag(D_s)==1); 

figure
x_values = 1:numberSamples;
pos = 80;
hold on
for i = 1:features 
    y_values = mean_x_vec_matrix(sampled_position(pos),:,i); 
    plot( x_values, y_values , 'LineWidth', 2 )
end
plot( x_values, gs(sampled_position(pos),:,1),'black.', 'LineWidth', 2)
plot( x_values, gs(sampled_position(pos),:,2),'black--', 'LineWidth', 2)
legend('Temp.', 'Wind Speed','Ground Truth Temp.','Ground Wind Speed','NumColumns', 2);
hold off
xlabel('Iteration, k')
ylabel('Amplitude')
set(gca,'gridlinestyle','-');
grid on

%%
sum(mean_elapsedTime_mat_comp)