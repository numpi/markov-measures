function Example3
%EXAMPLE3 Compute ???
%   Detailed explanation goes here

% Use the following to load Q
Q = loadArm2('data/example3/Experiment_80.arm');

pi0 = [ 1 , zeros(1, size(Q,2)-1) ];
r = double(full(sum(abs(Q')) == 0)');

% [pi0, Q, r] = deflate_zero_eigs(pi0, Q, r);

T = 40;

n = size(Q, 2);

tic; 
ee = funm_markov(pi0, Q, r, 'phi', T, 'alg', 'higham');
toc

% NG = 95, T = 30, 3.090795e-4
% NG = 95, T = 25, 5.257271e-7

end

