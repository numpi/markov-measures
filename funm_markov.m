function r = funm_markov(pi0, Q, v, f, t, varargin)
%FUNM_MARKOV Estimate the bilinear form pi0' * f(Q) * v. 
%
% R = FUNM_MARKOV(PI0, Q, V, F, T) estimates the bilinear form defined by the
%     matrix function f(Q). The function f needs to be expressed as a
%     string, and the valid options are: 
%     - F == 'exp': Matrix exponential exp(tQ)
%     - F == 'phi': The function (exp(tQ) - 1) / Q, properly extended by
%          continuity at 0. Evaluation is performed by computing the matrix
%          exponential of a larger matrix. 
% 
% R = FUNM_MARKOV(PI0, Q, V, F, T, 'tol', TOL) sets the tolerance for the
%     computation to the specified value. 
% 
% Author: Leonardo Robol <leonardo.robol@isti.cnr.it>

% Lambda parameter uesd for shifting
l = 0.0; % 0.01 * norm(Q, 1);

% Option parsing
p = inputParser;
addParameter(p, 'tol', 1e-8);

parse(p, varargin{:});
opts = p.Results;

if strcmp(f, 'exp') || strcmp(f, 'phi')
    param.function = 'exp';
end

if strcmp(f, 'phi_pade')
    error('This scheme is currently unsupported');
    param.function = 'expz'; % @(t, z) 1 ./ (t - z) .* exp(t) ./ (2 * pi);
end

param.stopping_accuracy = opts.tol;
param.waitbar = false;
param.verbose = false;
param.tol = opts.tol;
param.thick = [];
param.inner_product = @(x,y) x' * y;
param.V_full = false;
param.H_full = true;
param.restart_length = 10;
param.max_restarts = 10;
param.hermitian = false;
param.reorth_number = 0;
param.exact = [];
param.min_decay = 0.95;
param.truncation_length = inf;
param.bound = false;

if strcmp(f, 'exp') || strcmp(f, 'phi_pade')
    A = t * Q - l * speye(size(Q, 1));
    r = v;
end

if strcmp(f, 'phi')
    A = [ t*Q - l * speye(size(Q, 1)), v ; zeros(1, size(Q,2) + 1) ];
    r = [ zeros(size(Q,2),1) ; 1 ];
end

param = param_init_quad(param);
ff = funm_quad(A, r, param);

ff = ff(1 : size(Q, 2));

% Ensure positivity: no matter what, the vector ff is the product of a
% positive matrix (expm(tQ)) with a positive vector, the ones containing
% the rewards -- therefore we can enforce positivity of ff here. 
if min(v) >= 0
    ff = max(ff, 0);
end

if strcmp(f, 'phi') || strcmp(f, 'phi_pade')
    ff = t * ff;
end

r = pi0 * ff * exp(l);


end

