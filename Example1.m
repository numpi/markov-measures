function Example1
%EXAMPLE1 Test the computation of the inst. average number of clients
%
% This m-file contains the setup required to replicate Experiment 1, which
% compute the average number of clients at time 1 for a
% quasi-Birth-and-Death queue with rates rho_1 and rho_2 to move left and
% right, respectively. 
%
% The computation is done at time T = 1, and is repeated for the same model
% with an increasing number of states, ranging from 2^10 to 2^20. 

if ~exist('funm_quad', 'file')
    error('funm_quad toolbox not found. Please add it to your path.');
end

exps = 10 : 20;
times = zeros(1, length(exps));

for i = 1 : length(exps)
    n = 2^exps(i); 
    
    pi0 = zeros(1, n); 
    pi0(1) = 1;

    v = 0 : (n-1);
    v = v';

    rho1 = 0.4; 
    rho2 = 0.6;
    
    Q = spdiags(ones(n,1) * [ rho1, -rho1-rho2, rho2 ], -1:1, n, n);
    Q(1,1) = -rho2;
    Q(end,end) = -rho1;
    
    f = funm_markov(pi0, Q, v, 'exp', 1);
    tt = timeit(@() funm_markov(pi0, Q, v, 'exp', 1));
    
    fprintf('N = %d, time = %e (f = %e)\n', n, tt, f);
    times(i) = tt;
end

dlmwrite('clients.dat', [ 2.^(exps)', times' ], '\t');