function [ z_k cost_out_all] = LogR_for_ADMM_concise( z_k_in,const_intercept, k_const , Label, rho, omega,max_iteration, initial_learn_rate, A_orig ,x_weight,converge_percentage,Mask)

	x_weight = sparse(x_weight);
	A = kron(A_orig,x_weight.');
    learning_rate = initial_learn_rate;   
   
    z_k = z_k_in;
   
    total_gra = floor(max_iteration * 5);
    gra = 1;total_itr = 1;
	step_ahead = 1;
	
	x_weight_rep = repmat(x_weight.',size(A_orig,1),1); 
	g = sparse(1,size(A_orig,2)^2);
	%

    while(gra <= max_iteration && gra <= total_gra)
    
     if(gra == 1) 
	 
		aa = A * z_k + const_intercept;
		aa = -( Label .* ( aa));
		bb = (log( 1 + exp( aa ) ) );
			
        cost_out = omega * sum(bb) + rho/2 * (z_k + k_const)' * (z_k + k_const);
        cost_out_last = cost_out;
        z_k_last = z_k;
		 
     else
		  
        if(step_ahead ==1)
            prob = 1./( 1 + exp(aa) );
            b = -Label .* (1-prob);
            %g = omega * A' * b + (rho/2) * (z_k + k_const) ;

			x_weight_b = x_weight_rep .* repmat(b,1,size(x_weight_rep,2));
			A_x_b = x_weight_b.' * A_orig;
		
			%A_x_b = A_orig.' * x_weight_b;
			%A_x_b = A_x_b.';
			
			g = omega * A_x_b(:) + (rho/2) * (z_k + k_const);
				
            g_last = g;
		
        else
		
			g = g_last;   
		  
        end
		

        z_k = z_k - learning_rate * g;
        aa = A * z_k + const_intercept;
        aa = -( Label .* ( aa));
        bb = (log( 1 + exp( aa ) ) ); 
        cost_out = omega * sum(bb) + rho/2 * (z_k + k_const)' * (z_k + k_const);
 
	end
       
	 %step ahead
     if((cost_out_last - cost_out) > 0)
         learning_rate = learning_rate * 1.5;
         cost_out_all(gra) = cost_out;
         cost_out_last = cost_out;
         z_k_last = z_k;
         gra = gra+1;
         step_ahead = 1;
		 
	 %go back
     elseif((cost_out_last - cost_out) < 0)
         learning_rate = learning_rate * 0.5;
         z_k = z_k_last;
         cost_out = cost_out_last;
         step_ahead = 0;
         
	 %steat ahead learning rate no change
	 else
         cost_out_all(gra) = cost_out;
         cost_out_last = cost_out;
         z_k_last = z_k;
         gra = gra+1; 
         step_ahead = 1;
     end
        
     check_converge_itr = 10;
    if( gra > (check_converge_itr+2) )
        check_diff = cost_out_all( gra - 1 - check_converge_itr + 1  :  gra-1) - cost_out_all(gra-2 - check_converge_itr + 1  : gra - 2) ;
        percentage = abs( check_diff ) ./  cost_out_all(gra-2 - check_converge_itr +1  : gra-2);

		if( sum(percentage<converge_percentage) == check_converge_itr)
            gra;
            break;
        end
    end
      total_itr = total_itr+ 1;
  end
  

end


