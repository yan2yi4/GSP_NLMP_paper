function [ mean_MSD_vector, mean_x_vec_matrix, ...
           mean_elapsedTime_vector,mean_MSD_t ] = ...
            Ext_run_loop_NLMP( x_ref , ...
            D_s, U_f, stable_dist_para, alg_factor, alg_selection, max_it, ...
            runs_total, p,static,sigma,U,th)
[classes, steps] = size(x_ref); 
F = size(U_f, 2);
alpha = stable_dist_para(1);
gam = stable_dist_para(2);
sigma_diag = diag(sigma);
num_sampled = sum(sum(D_s));
MSD_vector_matrix = []; 
NMSD_vector_matrix = [];
MSD_t_runs_total_matrix = []; 
elapsedTime_vector_matrix = [];

x_vec_tensor = [];

if alg_selection ==4
    B_1 = alg_factor * U_f *pinv(U_f'*D_s*U_f)*U_f';
    B_2 = U_f*U_f'; 
    FLOM_noise_weight = FLOM( p, alpha,gam);
    aux_sum = ones(1,classes);
end

    
for j = 1:runs_total
    stable_dist = makedist('Stable','alpha',alpha,'beta',0,'gam',gam,'delta',0);
    a_noise = random(stable_dist,classes,max_it);
    x_noisy = x_ref + a_noise;
    MSD_t = []; MSD_vector = []; elapsedTime_vector = []; NMSD_vector = []; NMSD_t = [];
    x_vec = zeros(classes,1);
    x_vec_matrix = [];
    for i = 1:max_it

       
        switch(alg_selection)

            case 1 %LMS
                
                B_l = U*sigma_diag*U';
                y_k = D_s*x_noisy(:,i);
                tic;

                
                error_vec = y_k - x_vec;
                x_vec = x_vec + alg_factor * B_l * D_s* error_vec;
                elapsedTime = toc;
            case 2 %LMP
                
                B_l = U*sigma_diag*U';
                y_k = D_s*x_noisy(:,i);
                tic;
                signed_error_vec = (y_k - x_vec);
                error_vec = abs( signed_error_vec).^(p-1).*sign(signed_error_vec);
                x_vec = x_vec + alg_factor * B_l * D_s * error_vec;
                elapsedTime = toc;
            case 3 %NLMS
                tic;
                mu_Bn_matrix = alg_factor * U_f * inv( U_f' * D_s * U_f ) * U_f';
                error_vec = D_s*(x_noisy(:,i) - x_vec);
                x_vec = x_vec + mu_Bn_matrix * error_vec ;
                elapsedTime = toc;
            case 4 %NLMP
                threshold = th;
                y_k = x_noisy(:,i);
                B_l = U_f*inv(U_f'*D_s*(FLOM(p, alpha,  gam)^(p-2))*U_f)*U_f';
                tic;                            
                signed_error_vec = (y_k - x_vec);
                abs_error_vec = abs(signed_error_vec);
                abs_error_vec_2 = abs_error_vec.^(p-1);
                error_vec = abs_error_vec_2.*sign(signed_error_vec);
                if sum(D_s*abs_error_vec_2)<threshold

                    x_vec = x_vec + alg_factor * B_l * D_s * error_vec;
                else
                    M_inv = U_f'*D_s*diag((abs_error_vec).^(p-2))*U_f;
                    x_vec = x_vec + alg_factor * U_f*(M_inv\(U_f' * D_s *error_vec));
                end

                elapsedTime = toc;
        end          

        x_vec_matrix = [x_vec_matrix x_vec];
        MSD = sum((vecnorm(x_ref(:,i) - x_vec,2,2).^2));
        
        MSD_vector = [MSD_vector;                      MSD];

        
        if static ==0 
            NMSD = MSD/sum((vecnorm(x_ref(:,i),2,2).^2));
            NMSD_vector = [NMSD_vector;                      NMSD];
%             MSD_t = [MSD_t immse(x_vec_matrix(:,1:i),x_ref(:,1:i))];
%         MSD_t = [MSD_t sum((vecnorm(x_ref(:,1:i) - x_vec_matrix,2,1).^2))/(i*classes)];
        NMSD_t = [NMSD_t mean(NMSD_vector)];

        end
        elapsedTime_vector = [ elapsedTime_vector ;
                               elapsedTime ];
    end
    if static ==0  
        MSD_t_runs_total_matrix = [MSD_t_runs_total_matrix;NMSD_t];
    end
    MSD_vector_matrix = [MSD_vector_matrix MSD_vector];
    x_vec_tensor(:,:,j) = x_vec_matrix;
    elapsedTime_vector_matrix = [ elapsedTime_vector_matrix elapsedTime_vector ];   
end
mean_MSD_t = [];
mean_MSD_vector = zeros(max_it,1);
mean_elapsedTime_vector = zeros(max_it,1);
mean_x_vec_matrix = [];

for counter = 1:max_it
    if static ==0
        mean_MSD_t(counter) = mean( MSD_t_runs_total_matrix(:,counter) );
    end
    mean_MSD_vector(counter) = mean( MSD_vector_matrix(counter,:) )';
    mean_elapsedTime_vector(counter) = mean( elapsedTime_vector_matrix(counter,:) )';

end
for i = 1:classes
    for j = 1:steps
        mean_x_vec_matrix(i,j) = mean(x_vec_tensor(i,j,:));    
    end
end

end
