function [ W_rec ] = SLIM_aug( A,z_k,u_k,W_ini,beta_l2,lambda_l1,rho)
	W_rec = zeros(size(A,2),size(A,2));
    z_k_matrix = reshape(z_k,size(A,2),size(A,2)); z_k_matrix = z_k_matrix.';
    u_k_matrix = reshape(u_k,size(A,2),size(A,2)); u_k_matrix = u_k_matrix.';
	for i = 1: size(A,2)
        z_k_j = [z_k_matrix(1:(i-1),i) ; z_k_matrix(i+1:end,i)];
        u_k_j = [u_k_matrix(1:(i-1),i) ; u_k_matrix(i+1:end,i)];
        w_init = [W_ini(1:(i-1),i) ; W_ini(i+1:end,i)];
        a_j =  A(:,i);
        A_remain = [A(:,1:(i-1)) A(:,i+1 :end )];
        opts.rsL2 = beta_l2;
        opts.rFlag = 0;
        opts.init = 1;
        opts.x0 = w_init;
        a_j_in = [a_j ; sqrt(rho)*(z_k_j-u_k_j)];
        A_remain_in = [A_remain ; sqrt(rho)*eye(size(z_k_j,1),size(u_k_j,1) )] ;
        [w_j, funVal]=nnLeastR( A_remain_in , a_j_in,lambda_l1, opts);
        w_out = [w_j(1:(i-1));0;w_j(i :end)];  
        W_rec(:,i) = w_out;
    end
end

