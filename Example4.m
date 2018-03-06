function Example4
%EXAMPLE4 Estimate sensitivity for the model of Example1
%
% Here we consider the same example of Example1 -- but our focus is compute
% the derivative with respect to one of the parameters. To this aim, we
% rely on the function funm_markov_sensitivity. 

% We choose a small n for testing purposes
n = 128; 

pi0 = zeros(1, n); 
pi0(1) = 1;

v = 0 : (n-1);
v = v';

rho1 = 0.6; 
rho2 = 0.4;
Q = spdiags(ones(n,1) * [ rho1, -rho1-rho2, rho2 ], -1:1, n, n);
Q(1,1) = -rho2;
Q(end,end) = -rho1;

% Construct the derivative with respect to rho2
dQ = spdiags(ones(n,1) * [ -1, 1 ], 0:1, n, n);
dQ(end,end) = 0;

% We evaluate the function at Q and Q with rho2 slightly perturbed, in
% order to get an approximation of the change in the measure. 
h = 1e-7;
d0 = funm_markov(pi0, Q, v, 'exp', 1);
Q1 = Q + h * dQ;
d1 = funm_markov(pi0, Q1, v, 'exp', 1);

de = (d1 - d0) / h; 

% Compute the derivatives using matrix functions
ds = funm_markov_sensitivity(pi0, Q, v, 'exp', 1, dQ);

fprintf(' - Derivative computed by finite differences: %e\n', de);
fprintf(' - Derivative computed by our approach: %e\n', ds);

% Test the scalability of this approach when n increases
Ns = 2.^(8 : 18);
times = zeros(1, length(Ns));

for i = 1 : length(Ns)
    n = Ns(i);
    pi0 = zeros(1, n); 
    pi0(1) = 1;

    v = 0 : (n-1);
    v = v';

    rho1 = 0.6; 
    rho2 = 0.4;
    Q = spdiags(ones(n,1) * [ rho1, -rho1-rho2, rho2 ], -1:1, n, n);
    Q(1,1) = -rho2;
    Q(end,end) = -rho1;

    % Construct the derivative with respect to rho2
    dQ = spdiags(ones(n,1) * [ -1, 1 ], 0:1, n, n);
    dQ(end,end) = 0;
    
    times(i) = timeit(@() funm_markov_sensitivity(pi0, Q, v, 'exp', 1, dQ));
    fprintf(' - N = %d, time required: %e seconds\n', n, times(i));
end

dlmwrite('sensitivities.dat', [ Ns', times' ], '\t');


end

