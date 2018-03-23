function [r, ff] = funm_markov(pi0, Q, v, f, t, varargin)
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
addParameter(p, 'alg', 'quad');
addParameter(p, 'restarts', 15);

parse(p, varargin{:});
opts = p.Results;

alg = opts.alg; 

if strcmp(f, 'exp') || strcmp(f, 'phi') || strcmp(f, 'phi2')
    param.function = 'exp';
end

if strcmp(alg, 'higham')
	% f = @(x) expmv(1, Q', pi0', []);
else
	param.stopping_accuracy = opts.tol;
	param.waitbar = false;
	param.verbose = false;
	param.tol = opts.tol;
	param.thick = [];
	param.inner_product = @(x,y) x' * y;
	param.V_full = false;
	param.H_full = true;
	param.restart_length = opts.restarts;
	param.max_restarts = 25;
	param.hermitian = false;
	param.reorth_number = 0;
	param.exact = [];
	param.min_decay = 0.95;
	param.truncation_length = inf;
	param.bound = false;
	param = param_init_quad(param);
end

if strcmp(f, 'exp') || strcmp(f, 'phi_pade')
    A = t * Q';
    r = pi0';
end

if strcmp(f, 'phi')
    A = [ t*Q', pi0' ; zeros(1, size(Q,2) + 1) ];
    r = [ zeros(size(Q,2),1) ; 1 ];
end

if strcmp(f, 'phi2')
    A = [ t*Q', pi0' zeros(size(Q, 2), 1)  ; ...
        zeros(1, size(Q,2) + 1) , 1 ; ...
        zeros(1, size(Q, 2) + 2) ];
    
    r = [ zeros(size(Q,2),1) ; 0 ; 1 ];
end

if strcmp(alg, 'higham')
	ff = expmv(1.0, A, r, [], 'double');
else
	ff = funm_quad(A, r, param);
end

ff = ff(1 : size(Q, 2));

% Ensure positivity: no matter what, the vector ff is the product of a
% positive matrix (expm(tQ)) with a positive vector, the ones containing
% the rewards -- therefore we can enforce positivity of ff here. 
% ff = max(ff, 0);

if strcmp(f, 'phi') || strcmp(f, 'phi_pade') || strcmp(f, 'phi2')
    ff = t * ff;
end

if strcmp(f, 'phi2')
    ff = t * ff;
end

r = v' * ff * exp(l);


end

