function d = funm_markov_sensitivity(pi0, Q, v, f, T, dQ)
%FUNM_MARKOV_SENSITIVITY Compute the sensitivity of a measure. 
%
% D = FUNM_MARKOV_SENSITIVITY(PI0, Q, V, F, T, DQ) is the derivative of the
%     specified measure with respect to a step in the direction DQ. In
%     practice, 
%
%       F(Q + h * DQ) = F(Q) + h * D + O(h^2), 
%
%     where F = @(Q) FUNM_MARKOV(PI0, Q, V, F, T) and D is the output of
%     this function. The meaning of the other parameters is the same of the
%     function FUNM_MARKOV. 
%
% Author: Leonardo Robol <leonardo.robol@isti.cnr.it>

% Build the agumented matrix
QQ = [ Q , dQ ; sparse(size(Q,1), size(Q,2)) , Q ]; 

d = funm_markov([ pi0, zeros(1,size(Q,2)) ], QQ, ...
        [ zeros(size(Q,2), 1) ; v ], f, T);

end

