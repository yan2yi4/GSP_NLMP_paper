clc;clear;
% close all;
load('sensor_graph_50_nodes.mat')
clear D_s

runs_total = 100;

rng(123456);%for reproducibility 

nun_unsampled = 20;
L = diag(sum(A)) - A;
[U,D] = eig(L); 
n_nodes = 50;
numberSamples = 6000;

F = n_nodes-nun_unsampled-10;        
M = n_nodes-nun_unsampled;

gs = node_values';
freq = U'*gs;
alpha = [2,1.8,1.5,1.2];
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

%%%
% a = U_f * s_f; %code representation
% b = [U_f zeros(n_nodes,F)] * [s_f;zeros(F,1)]; %sparse representation
% %sparsity worksbecause a-b is 0

%BANDLIMITED Graph signal~~~~~~~~~~~~~~~~
gs_BL = U_f*s_f;
gs = gs_BL;

names = {'\alpha = 2','\alpha = 1.8','\alpha = 1.5','\alpha = 1.2'};
u_lmp = 0.005;
u_nlmp = 0.00375;
alg_selection_vec = [4 4 4 4 2 2 2 2 ];  
alg_param_vec = [u_nlmp u_nlmp*0.875 u_nlmp*0.675 u_nlmp*0.675 u_lmp u_lmp u_lmp u_lmp ]; 

[~,num_algs] = size(alg_param_vec);
max_it = numberSamples;

amp_mat = ones(1,max_it) ;
gs = gs*amp_mat;

D_s = eig_sampling_strategy( M, U_f);

mean_MSD_mat_comp = [];

mean_elapsedTime_mat_comp = [];
D_s = double(D_s);
NLMP_counter = 1;
LMP_counter = 1;

for i = 1:num_algs 
    if alg_selection_vec(i) == 4
        threshold = FLOM(p(NLMP_counter)-1, alpha(NLMP_counter),gam)*M;
        [mean_MSD_mat_comp(:,i), mean_x_vec_matrix(:,:,i), mean_elapsedTime_mat_comp(:,i),~] = ...
        Ext_run_loop_NLMP( gs, D_s, U_f, [alpha(NLMP_counter) gam],alg_param_vec(i), alg_selection_vec(i), ...
            max_it, runs_total, p(NLMP_counter),1,sigma_diag,U,threshold) ;
        NLMP_counter = NLMP_counter+1;
    end
    if alg_selection_vec(i) == 2
        threshold = FLOM(p(LMP_counter)-1, alpha(LMP_counter),gam)*M;
        alpha_stable = makedist('Stable','alpha',alpha(LMP_counter),'beta',0,'gam',gam,'delta',0);
        [mean_MSD_mat_comp(:,i), mean_x_vec_matrix(:,:,i), mean_elapsedTime_mat_comp(:,i),~] = ...
        Ext_run_loop_NLMP( gs, D_s, U_f, [alpha(LMP_counter) gam],alg_param_vec(i), alg_selection_vec(i), ...
            max_it, runs_total, p(LMP_counter),1,sigma_diag,U,threshold) ;
        LMP_counter = LMP_counter+1;
    end
end

    % ------------------------------------------------------------------- %
%%
figure
hold on
ColOrd = get(gca,'ColorOrder');
sampled_position = find(diag(D_s)==1); 
pos = 7;
x_values = 1:numberSamples;
for i = 1:num_algs
    y_values = mean_x_vec_matrix(sampled_position(pos),:,i); 
    if i>4
        plot(x_values, y_values,'Color',ColOrd(i-4,:))
        hline = findobj(gcf, 'type', 'line');
        set(hline(1),'LineStyle',':');    
    else
        plot( x_values, y_values)
    end
end
plot( x_values, gs(sampled_position(pos),:) ,'Color','black')
leg = legend([names, 'Target Node'],'Location','southeast','NumColumns',2);
hold off
xlabel('Iteration, k')
ylabel('Amplitude')
grid on

figure
hold on
ColOrd = get(gca,'ColorOrder');
for i = 1:num_algs
    y_values = mean_MSD_mat_comp(x_values,i);   
    if i>4
        plot(x_values,10*log10( y_values),'Color',ColOrd(i-4,:))
        hline = findobj(gcf, 'type', 'line');
        set(hline(1),'LineStyle',':');    
    else
        plot( x_values, 10*log10( y_values ))
    end
end
leg = legend(names,'NumColumns',2);
hold off
xlabel('Iteration, k')
ylabel('MSD(dB)')
grid on
%%
sum(mean_elapsedTime_mat_comp)