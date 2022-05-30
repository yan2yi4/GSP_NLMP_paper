clc;clear;
% close all;
load('sensor_graph_50_nodes.mat')
clear D_s

runs_total = 100;

rng(12345);%for reproducibility 

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
%alpha_stable = makedist('Stable','alpha',lpha,'beta',0,'gam',gam,'delta',0);

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

%BANDLIMITED NOW~~~~~~~~~~~~~~~~
gs_BL = U_f*s_f;
gs = gs_BL;


%%%
% a = U_f * s_f; %code representation
% b = [U_f zeros(n_nodes,F)] * [s_f;zeros(F,1)]; %sparse representation
% %sparsity worksbecause a-b is 0

names = {'\mu = 0.05','\mu = 0.01','\mu = 0.005'};

alg_selection_vec = [4 4 4 2 2 2];  
alg_param_vec = [0.05 0.01 0.005 0.05 0.01 0.005]; 



[~,num_algs] = size(alg_param_vec);
max_it = numberSamples;

amp_mat = ones(1,max_it) ;
gs = gs*amp_mat;


% D_s = sampling_strategy( M, U_f,R);
D_s = eig_sampling_strategy( M, U_f);

mean_MSD_mat_comp = [];

mean_elapsedTime_mat_comp = [];
D_s = double(D_s);
threshold = FLOM( p-1, alpha,gam)*M;
for i = 1:num_algs 
    [mean_MSD_mat_comp(:,i), ...
        mean_x_vec_matrix(:,:,i), mean_elapsedTime_mat_comp(:,i),~] = ...
      Ext_run_loop_NLMP( gs, ... 
            D_s, U_f, [alpha gam], ...
            alg_param_vec(i), alg_selection_vec(i), ...
            max_it, runs_total, p,1,sigma_diag,U,threshold) ;

end

    % ------------------------------------------------------------------- %
%%
sampled_position = find(diag(D_s)==1); 
x_values = 1:numberSamples;

% B_l = U_f*inv(U_f'*D_s*(FLOM(p, alpha,  gam)^(p-2))*U_f)*U_f';
% % B_l = U_f*inv(U_f'*D_s*(FLOM(p-2, alpha,  gam))*U_f)*U_f';
% vec_I = eye(n_nodes);
% vec_I = vec_I(:);
% % C = (0.95^2)*(FLOM( 0.95*2-2, 1,gam))*eye(n_nodes);
% % C = ones(n_nodes);
% C = FLOM( 2*p-2, alpha,gam);
% G = B_l'*D_s*C*D_s*B_l;
% % R = (p_sign-1)*FLOM( p_sign-2, 1,gam);
% R = FLOM( p-2, alpha,gam);

% for i = 1:num_algs
%     Q = kron(eye(n_nodes)-alg_param_vec(i)*B_l*D_s*R,eye(n_nodes)-alg_param_vec(i)*B_l*D_s*R);
%     [Q_dim,~] = size(Q); 
%     MSD(i) = alg_param_vec(i)^2*(G(:)'*(pinv(eye(Q_dim)-Q)*vec_I));
% end

figure
ColOrd = get(gca,'ColorOrder');

hold on
for i = 1:3
    plot( x_values, 10*log10( mean_MSD_mat_comp(x_values,i) ),'Color',ColOrd(i,:));
end
% for i = 1:num_algs
%     plot(x_values, 10*log10(MSD(i))*ones(numberSamples,1),'Color',ColOrd(i,:),'LineStyle',':');
% end
for i = 4:6
    plot( x_values, 10*log10( mean_MSD_mat_comp(x_values,i) ),'Color',ColOrd(i-3,:),'LineStyle',':');
end
leg = legend([names,'','','']);

grid on
hold off
xlabel('Iteration, k')
ylabel('MSD(dB)')


%%

% G = graph(A);
% G.Nodes.Value = node_values';
% figure()
% hight = width;
% set(0,'defaultFigurePosition',[0 0 width hight]);
% top = 0.5;  % normalized top margin
% bottom = 0.5;	% normalized bottom margin
% left = 0.5;	% normalized left margin
% right = 0.5;  % normalized right margin
% set(0,'defaultAxesPosition',[left/width bottom/hight (width-left-right)/width  (hight-bottom-top)/hight]);
% plt = plot(G,'NodeLabelMode','manual');
% plt.NodeCData = G.Nodes.Value;
% plt.NodeLabel = {};
% box off
% set(gca,'Visible','off')
% colorbar
% colormap('jet')
%%
sum(mean_elapsedTime_mat_comp)