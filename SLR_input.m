function [  ] = SLR_input( InputPath, OutputPath, OMEGA, ALPHA, LAMBDA, BETA, GAMMA)
PHI = 1;
OMEGA_p = OMEGA;
OMEGA_n = OMEGA;
LAMBDA_l1_p = LAMBDA;
LAMBDA_l1_n = LAMBDA;
BETA_l2_p = BETA;
BETA_l2_n = BETA;
GAMMA_l1 = GAMMA;
ALPHA_l2 = ALPHA;
maxNumCompThreads(1);
run('SLR.m');
end

