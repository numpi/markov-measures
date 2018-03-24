function Example3
%EXAMPLE3 Compute ???
%   Detailed explanation goes here

% Maximum index of the experiment to load.
N = 10;

times = zeros(1, N);
sizes = zeros(1, N);

T = 40;

% Use the following to load Q
for i = 1 : N
    Q = loadArm2(sprintf('data/example3/Experiment_%d.arm', i));
    n = size(Q, 2);

    % Construct the initial probability and the reward vector
    pi0 = [ 1 , zeros(1, n-1) ];

    % You might want to use this vector to mark all the absorbing states
    % r = double(full(sum(abs(Q')) == 0)');

    r = zeros(n, 1);
    r(dlmread(sprintf('data/example3/r_%d.dat', i))) = 1;

    tic;
    ee = funm_markov(pi0, Q, r, 'phi', T, ...
        'alg', 'quad', 'restarts', 50);
    times(i) = toc;

    sizes(i) = n;
    fprintf(' - N = %d, time = %d (f = %e)\n', n, times(i), ee);
end

dlmwrite('example3.dat', [ sizes', times' ], '\t');

end
