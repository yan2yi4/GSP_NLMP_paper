function [ mean_MSD_vector, mean_x_vec_matrix, ...
           mean_elapsedTime_vector,mean_MSD_t ] = ...
            Ext_run_loop_NLMP_mult( x_ref , ...
            D_s, U_f, stable_dist_para, alg_factor, max_it, ...
            ensemble, p,static,sigma,U,th)
[classes, steps,dims] = size(x_ref); 
F = size(U_f, 2);
alpha = stable_dist_para(1);
gam = stable_dist_para(2);
sigma_diag = diag(sigma);
num_sampled = sum(sum(D_s));
MSD_vector_matrix = []; 
NMSD_vector_matrix = [];
MSD_t_ensemble_matrix = []; 
elapsedTime_vector_matrix = [];

x_vec_tensor = [];
for j = 1:ensemble
    stable_dist = makedist('Stable','alpha',alpha,'beta',0,'gam',gam,'delta',0);
    MSD_t = []; MSD_vector = []; elapsedTime_vector = [];
    x_vec = zeros(classes,dims);
    x_vec_matrix = [];
    for i = 1:max_it
        x_current_obs = zeros(classes,dims);
        for k = 1:dims
            x_current_obs(:,k) = D_s*squeeze(x_ref(:,i,k))+random(stable_dist,classes,1);
        end
                threshold = th;
                B_l = U_f*inv(U_f'*D_s*(FLOM(p, alpha,  gam)^(p-2))*U_f)*U_f';
                tic;                            
                signed_error_vec = (x_current_obs - x_vec);
                abs_error_vec = abs(signed_error_vec);
                abs_error_vec_2 = abs_error_vec.^(p-1);
                error_vec = D_s*abs_error_vec_2.*sign(signed_error_vec);
                x_vec = x_vec +  B_l * error_vec * alg_factor;

                elapsedTime = toc;
        x_vec_matrix = cat(3,x_vec_matrix,x_vec);
        MSD = sum((vecnorm(squeeze(x_ref(:,i,:)) - x_vec,2,2).^2));
        MSD_vector = [MSD_vector;                      MSD];
        if static ==0  
        MSD_t = [MSD_t mean(MSD_vector)];

        end
        elapsedTime_vector = [ elapsedTime_vector ;
                               elapsedTime ];
    ed
    if static ==0  %time-varying
        MSD_t_ensemble_matrix = [MSD_t_ensemble_matrix;MSD_t];
    end
    MSD_vector_matrix = [MSD_vector_matrix MSD_vector];
    x_vec_tensor(:,:,:,j) = x_vec_matrix;
    elapsedTime_vector_matrix = [ elapsedTime_vector_matrix elapsedTime_vector ];   
end
mean_MSD_t = [];
mean_MSD_vector = zeros(max_it,1);
mean_elapsedTime_vector = zeros(max_it,1);
mean_x_vec_matrix = [];
mean_x_vec_matrix = mean(x_vec_tensor,4);
mean_x_vec_matrix = permute(mean_x_vec_matrix,[1 3 2]);

for counter = 1:max_it
    if static ==0
        mean_MSD_t(counter) = mean( MSD_t_ensemble_matrix(:,counter) );
    end
    mean_MSD_vector(counter) = mean( MSD_vector_matrix(counter,:) )';
    mean_elapsedTime_vector(counter) = mean( elapsedTime_vector_matrix(counter,:) )';

end

end
