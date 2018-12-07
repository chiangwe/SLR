function [ W_rec ] = SLIM( A,beta_l2,lambda_l1)
	W_rec = zeros(size(A,2),size(A,2));
	for i = 1: size(A,2)
        a_j =  A(:,i);
        A_remain = [A(:,1:(i-1)) A(:,i+1 :end )];
        opts.rsL2 = beta_l2;
        opts.rFlag = 0;
        [w_j, funVal]=nnLeastR( A_remain , a_j,lambda_l1, opts);
        w_out = [w_j(1:(i-1));0;w_j(i :end)];  
        W_rec(:,i) = w_out;
    end
end

