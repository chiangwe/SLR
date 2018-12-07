%%Version 2018 12 06
%clear all; clc;
addpath(genpath('./'))

%%%%%%%% meta parameters%%%%%%%%

%InputPath
display(InputPath)
load(InputPath);
size(train_data);

Input_data = train_data;

tic;

%%%%%%%%get training data and texting data%%%%
%train_index = setdiff((1:2000),cell2mat(TEST(fold)));
%TEST_all = Input_data(cell2mat(TEST(fold)),:);

TRAIN_all = Input_data;
Label = [ones(size(TRAIN_all,1)/2,1);-ones(size(TRAIN_all,1)/2,1)];
Label_p = [ ones(size(TRAIN_all,1)/2,1)];
Label_n = [-ones(size(TRAIN_all,1)/2,1)];

A = TRAIN_all;
inst_num = size(A,1);

A_p = A(              1 : inst_num/2 , :);
A_n = A( inst_num/2 + 1 : inst_num   , :);

Mask = sparse(A);
Mask_p = sparse(A_p);
Mask_n = sparse(A_n);

%Mask = sparse(ones(size(A,1),size(A,2)));
phi = PHI;
omega_p = OMEGA_p; omega_n = OMEGA_n;
lambda_l1_p = LAMBDA_l1_p; lambda_l1_n = LAMBDA_l1_n;
beta_l2_p = BETA_l2_p; beta_l2_n = BETA_l2_n;
gamma_l1 = GAMMA_l1; alpha_l2 = ALPHA_l2;

%%change initial point
[ x_weight,c_intercept ] = LogR( A, Label ,alpha_l2,gamma_l1);
x_weight = sparse(x_weight);

[ W_rec_p ] = SLIM(A_p,beta_l2_p,lambda_l1_p);
W_rec_p = sparse(W_rec_p);

[ W_rec_n ] = SLIM(A_n,beta_l2_n,lambda_l1_n);
W_rec_n = sparse(W_rec_n);



RHO_admm_p = 100; rho_initial_p = RHO_admm_p;
RHO_admm_n = 100; rho_initial_n = RHO_admm_n;
Alter_itr= 1; admm_itr = 700; admm_itr_min = 100; check_converge_itr = 5; max_iteration = 100;
Mu_for_ADMM = 10; tao_incr = 2; tao_decr = 2;
abs_tolerance = 1*10^(-6); rel_tolerance = 1*10^(-6);
initial_learn_rate = 0.01;   converge_percentage	= 0.0001; 
        

cost_total =  (1/2) * norm( A_p - A_p * W_rec_p,'fro')^2 + (1/2) * phi * norm( A_n - A_n * W_rec_n,'fro')^2 ...
            + omega_p * sum (log( 1 + exp( -( Label_p .* ( x_weight.' * ( ( ( A_p * W_rec_p ) .' ) ) + c_intercept ).' ) ) ) ) ...
            + omega_n * sum (log( 1 + exp( -( Label_n .* ( x_weight.' * ( ( ( A_n * W_rec_n ) .' ) ) + c_intercept ).' ) ) ) ) ...
            + beta_l2_p / 2 * norm( W_rec_p , 'fro' )^2 + beta_l2_n / 2 * norm( W_rec_n , 'fro' )^2 ...
            + lambda_l1_p * norm( W_rec_p ,1) + lambda_l1_n * norm( W_rec_n ,1) ...
            + alpha_l2 / 2 * norm(x_weight,'fro')^2 + gamma_l1 * norm(x_weight,1);
        
SLIM_total = (1/2) * norm( A_p - A_p * W_rec_p,'fro')^2 + (1/2) * phi * norm( A_n - A_n * W_rec_n,'fro')^2 ...
             + beta_l2_p / 2 * norm( W_rec_p , 'fro' )^2 + beta_l2_n / 2 * norm( W_rec_n , 'fro' )^2 ...
             + lambda_l1_p * norm( W_rec_p ,1) + lambda_l1_n * norm( W_rec_n ,1);
         
LogR_total =  omega_p * sum (log( 1 + exp( -( Label_p .* ( x_weight.' * ( ( ( A_p * W_rec_p ) .' )  ) + c_intercept ).' ) ) ) ) ...
            + omega_n * sum (log( 1 + exp( -( Label_n .* ( x_weight.' * ( ( ( A_n * W_rec_n ) .' )  ) + c_intercept ).' ) ) ) ) ...
            + alpha_l2 / 2 * norm(x_weight,'fro')^2 + gamma_l1 * norm(x_weight,1);

cost_all = zeros(2*Alter_itr+1,1);
cost_record_all = zeros(Alter_itr+1,1);
SLIM_all = zeros(2*Alter_itr+1,1);
LogR_all = zeros(2*Alter_itr+1,1);

cost_all(1) =   cost_total; 
SLIM_all(1) =   SLIM_total;
LogR_all(1) =   LogR_total; 

string_out = ['Start cost with : ' num2str(cost_total) '\n'];
display(string_out)

itr = 1;
cost_all_admm = zeros(admm_itr+1,1);
cost_all_admm(1) = cost_total;

%clc;
%A_in = kron(A,x_weight.');
W_rec_old_p = W_rec_p; W_rec_old_n = W_rec_n; 
x_weight_old = x_weight; c_intercept_old = c_intercept;

%%%%%%%%%%%%% FIX x.c solve W ADMM %%%%%%%%%%%%%%%%%%%%
 
vec_Wt_p = reshape(W_rec_p.',size(W_rec_p,1)*size(W_rec_p,2),1);
vec_Wt_n = reshape(W_rec_n.',size(W_rec_n,1)*size(W_rec_n,2),1);

%%%%%%% ADMM start point %%%%%%%%%
w_k_p = W_rec_p;
z_k_p = vec_Wt_p;
u_k_p = zeros(size(vec_Wt_p,1),size(vec_Wt_p,2));

w_k_n = W_rec_n;
z_k_n = vec_Wt_n;
u_k_n = zeros(size(vec_Wt_n,1),size(vec_Wt_n,2));

do = true; 

itr_AMDD = 1;
z_k_old_p = z_k_p;
w_k_old_p = w_k_p;
u_k_old_p = u_k_p; 

z_k_old_n = z_k_n;
w_k_old_n = w_k_n;
u_k_old_n = u_k_n; 


rho_p = rho_initial_p;
rho_n = rho_initial_n;

%%
while(do && itr_AMDD <= admm_itr )
    
    %%% calculate error p

    error_pri_p  = sqrt( size( z_k_p , 1 ) ) * abs_tolerance + rel_tolerance *  max( norm( w_k_p ,'fro') , norm( z_k_p , 'fro' ) );
    error_dual_p = sqrt( size( u_k_p , 1 ) ) * abs_tolerance + rel_tolerance * norm( rho_p * u_k_p ,'fro');

    %%% calculate error n
  
    error_nri_n  = sqrt( size( z_k_n , 1 ) ) * abs_tolerance + rel_tolerance *  max( norm( w_k_n ,'fro') , norm( z_k_n , 'fro' ) );
    error_dual_n = sqrt( size( u_k_n , 1 ) ) * abs_tolerance + rel_tolerance * norm( rho_n * u_k_n ,'fro');

    %%%% update z^k+1
    
    s_p = reshape( w_k_p.' , size( w_k_p ,1) * size(w_k_p,2),1) + u_k_p;
    k_p = -s_p;
    
    s_n = reshape( w_k_n.' , size( w_k_n ,1) * size(w_k_n,2),1) + u_k_n;
    k_n = -s_n;

    %%%%%%%%%%%%%%%%

    %%%%% z function input %%%%%

    initial_learn_rate = 0.01;  converge_percentage	= 0.0001; 
    %%%%%%%%%%%%%%%%%%%%%%;
    
    %%%%%%p%%%%%%
    [ z_k_p cost_out_all_p] = LogR_for_ADMM_concise( z_k_p, c_intercept, k_p , Label_p, rho_p , omega_p , max_iteration , initial_learn_rate , A_p , x_weight , converge_percentage  );
	z_k_p = sparse(z_k_p);

    [ w_k_p ] = SLIM_aug( A_p , z_k_p , u_k_p , w_k_p , beta_l2_p , lambda_l1_p , rho_p );
	w_k_p = sparse(w_k_p);

    %%%%%%n%%%%%%
 
    [ z_k_n cost_out_all_n] = LogR_for_ADMM_concise( z_k_n, c_intercept, k_n , Label_n, rho_n , omega_n , max_iteration , initial_learn_rate , A_n , x_weight , converge_percentage  );
	z_k_n = sparse(z_k_n);

    [ w_k_n ] = SLIM_aug( A_n , z_k_n , u_k_n , w_k_n , beta_l2_n , lambda_l1_n , rho_n );
	w_k_n = sparse(w_k_n);
    
    %%%%%%%%%%%%% FIX W solve  x.c logR %%%%%%%%%%%%%%%%%%%% 
    opts.rsL2 = alpha_l2/omega_p;
    opts.rFlag = 0;
    opts.c0 = c_intercept;
    opts.x0 = x_weight;
    opts.init=1;
    opts.mFlag = 0;
    opts.lFlag = 0;
    opts.sWeight = [1;1*omega_n/omega_p];
    
    %A_w_k_Mask = [ (A_p * w_k_p) .* Mask_p ;  ( A_n * w_k_n ) .* Mask_n ];
    A_w_k = [ (A_p * w_k_p);  ( A_n * w_k_n )];
    
    [ x_weight, c_intercept ] = LogisticR( A_w_k , Label ,gamma_l1/omega_p, opts);

    u_k_p = u_k_p + reshape(w_k_p.',size(w_k_p,1)*size(w_k_p,2),1) - z_k_p;
    u_k_p = sparse(u_k_p);

    u_k_n = u_k_n + reshape(w_k_n.',size(w_k_n,1)*size(w_k_n,2),1) - z_k_n;
    u_k_n = sparse(u_k_n);
    
    R_k_p = norm(reshape(w_k_p.',size(w_k_p,1)*size(w_k_p,2),1)-z_k_p,'fro');
    S_k_p = norm(z_k_p-z_k_old_p,'fro');
    
    R_k_n = norm(reshape(w_k_n.',size(w_k_n,1)*size(w_k_n,2),1)-z_k_n,'fro');
    S_k_n = norm(z_k_n-z_k_old_n,'fro');
	
    %%%%rho update%%%%
    if(R_k_p > Mu_for_ADMM * S_k_p )
		rho_p = tao_incr * rho_p; 
		u_k_p = u_k_p/tao_incr;
    elseif(S_k_p > Mu_for_ADMM * R_k_p)
		rho_p = rho_p / tao_decr; 
		u_k_p = u_k_p * tao_decr;
    else
		rho_p = rho_p;
    end
    %%%%%%%%%%%%%%%%%%
    %%%%rho update%%%%
    if(R_k_n > Mu_for_ADMM * S_k_n )
		rho_n = tao_incr * rho_n; 
		u_k_n = u_k_n/tao_incr;
    elseif(S_k_n > Mu_for_ADMM * R_k_n)
		rho_n = rho_n / tao_decr; 
		u_k_n = u_k_n * tao_decr;
    else
		rho_n = rho_n;
    end
    %%%%%%%%%%%%%%%%%%
    
    
    W_rec_p = w_k_p;
    W_rec_n = w_k_n;
    
    cost_total =  (1/2) * norm( A_p - A_p * W_rec_p,'fro')^2 + (1/2) * phi * norm( A_n - A_n * W_rec_n,'fro')^2 ...
            + omega_p * sum (log( 1 + exp( -( Label_p .* ( x_weight.' * ( ( ( A_p * W_rec_p ) .' )  ) + c_intercept ).' ) ) ) ) ...
            + omega_n * sum (log( 1 + exp( -( Label_n .* ( x_weight.' * ( ( ( A_n * W_rec_n ) .' )  ) + c_intercept ).' ) ) ) ) ...
            + beta_l2_p / 2 * norm( W_rec_p , 'fro' )^2 + beta_l2_n / 2 * norm( W_rec_n , 'fro' )^2 ...
            + lambda_l1_p * norm( W_rec_p ,1) + lambda_l1_n * norm( W_rec_n ,1) ...
            + alpha_l2 / 2 * norm(x_weight,'fro')^2 + gamma_l1 * norm(x_weight,1);

    cost_all_admm(itr_AMDD+1) = cost_total;

	
   if(itr_AMDD > 10)
		check_diff = cost_all_admm(itr_AMDD+1 - check_converge_itr +1  : itr_AMDD+1) - cost_all_admm(itr_AMDD+1 - check_converge_itr   : itr_AMDD) ;
		percentage = abs( check_diff ) ./  cost_all_admm(itr_AMDD+1 - check_converge_itr   : itr_AMDD);
		string_out = ['  cost_ADMM: ' num2str(cost_total) ' rho_p:' num2str(rho_p) '  itr_AMDD: ' num2str(itr_AMDD) '  R_k_p: ' num2str(R_k_p) '  error_pri: ' num2str(error_pri_p) '  S_k_p: ' num2str(S_k_p) '  error_dual_p: ' num2str(error_dual_p) '  convergence percentage ' num2str(percentage.'*100,'%1.4f\t') '\n'];
		display(string_out);
		fprintf(fid, string_out);
	
		if( sum(percentage<0.0005) == check_converge_itr && itr_AMDD >=check_converge_itr+2 && itr_AMDD >=admm_itr_min)
			break;
		end
	else
		string_out = ['  cost_ADMM: ' num2str(cost_total) ' rho_p:' num2str(rho_p) '  itr_AMDD: ' num2str(itr_AMDD) '  R_k: ' num2str(R_k_p) '  error_pri_p: ' num2str(error_pri_p) '  S_k_p: ' num2str(S_k_p) '  error_dual_p: ' num2str(error_dual_p) '\n'];
		display(string_out);
	%	fprintf(fid, string_out);
   end		

	z_k_old_p = z_k_p ;u_k_old_p  = u_k_p ; w_k_old_p  = w_k_p ; 
    z_k_old_n = z_k_n ;u_k_old_n  = u_k_n ; w_k_old_n  = w_k_n ;
    
    itr_AMDD = itr_AMDD+1;
	 
end %%%ADMM



time_out =toc;
string_out = ['Compelte time_out: ' num2str(time_out) '\n'];
display(string_out)
fclose(fid);

%exit()
save( OutputPath,'W_rec_p','W_rec_n','x_weight','c_intercept','cost_all_admm','itr_AMDD');

exit();
