clc;clear;
% close all;
load('sensor_graph_50_nodes.mat')
clear D_s

runs_total = 100;

rng(846213);%for reproducibility 

nun_unsampled = 20;
L = diag(sum(A)) - A;
[U,D] = eig(L); 
n_nodes = 50;
numberSamples = 1200;

F = n_nodes-nun_unsampled-10;        
M = n_nodes-nun_unsampled;

gs = node_values';
freq = U'*gs;

alpha = 1.5;
gam = 0.1;

p = alpha-0.05;

[sorted_values, sorted_indices] = sort(abs(freq),'descend'); 


used_indices = sort( sorted_indices(1:F), 'ascend' );
sigma_diag = zeros(n_nodes,1);
sigma_diag(used_indices) = 1;

U_f = [];
s_f = [];
for counter = 1:F
    U_f = [U_f U( : , used_indices(counter) ) ];
    s_f = [s_f; freq( used_indices(counter) )];
end

%BANDLIMITED GS~~~~~~~~~~~~~~~~
gs_BL = U_f*s_f;
gs = gs_BL;


%%%
% a = U_f * s_f; %code representation
% b = [U_f zeros(n_nodes,F)] * [s_f;zeros(F,1)]; %sparse representation
% %sparsity worksbecause a-b is 0

names = {'Equation (12)', 'Algorithm 1' 'Equation (13)'};
u_nlmp = 0.01;
alg_selection_vec = 4; 


num_algs = 3;
max_it = numberSamples;

amp_mat = ones(1,max_it) ;
gs = gs*amp_mat;

D_s = eig_sampling_strategy( M, U_f);

mean_MSD_mat_comp = [];

mean_elapsedTime_mat_comp = [];
D_s = double(D_s);

threshold = [0, FLOM( p-1, alpha,gam)*M, inf];
for i = 1:num_algs 
    [mean_MSD_mat_comp(:,i), ...
        mean_x_vec_matrix(:,:,i), mean_elapsedTime_mat_comp(:,i),~] = ...
      Ext_run_loop_NLMP( gs, ... 
            D_s, U_f, [alpha gam], ...
            u_nlmp, alg_selection_vec, ...
            max_it, runs_total, p,1,sigma_diag,U,threshold(i)) ;
end

    % ------------------------------------------------------------------- %
%%
sampled_position = find(diag(D_s)==1); 
hold on
x_values = 1:numberSamples;

for i = 1:num_algs
    y_values = mean_MSD_mat_comp(x_values,i);
    plot( x_values, 10*log10( y_values ));
end

leg = legend(names);
set(gca,'gridlinestyle','-');
grid on
hold off
xlabel('Iteration, k')
ylabel('MSD(dB)')
sum(mean_elapsedTime_mat_comp)
rng('default')