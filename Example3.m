function Example3(Q)
%EXAMPLE3 Compute ???
%   Detailed explanation goes here

% Use the following to load Q
% Q = loadArm('example3/Example3/Transformer/SS/Experiment_10.arm');
pi0 = [ 1 , zeros(1, size(Q,2)-1) ];
r = double(full(sum(abs(Q')) == 0)');

T = 10;

tic; ee = funm_markov(pi0, Q, r, 'exp', T); toc
ee


end

