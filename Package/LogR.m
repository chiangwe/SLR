function [ x_weight c_intercept ] = LogR( A, Label ,alpha_l2,gamma_l1)
	        opts.rsL2 = alpha_l2;
            opts.rFlag = 0;
            [x_weight, c_intercept, funVal]=LogisticR(A, Label ,gamma_l1, opts);
end

