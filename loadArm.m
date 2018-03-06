function Q = loadArm(armfile)
%LOADARM Loads a Q matrix from the ARM format of Moebius. 
%
% Q = LOADARM(ARMFILE) loads the infinitesimal generator Q of a Markov
% chains stored in the .arm format used by Moebius. Therefore, you can use
% this function to interface FUNM_MARKOV and FUNM_MARKOV_SENSITIVITY with
% models obtained working in that environment. 
%
% This function is _not_ tuned for performance. In particular, it might be
% quite slow for larger files. 
%
% Author: Leonardo Robol <leonardo.robol@isti.cnr.it>

f = fopen(armfile);

n = fscanf(f, '%d', 1);
n = fscanf(f, '%d', 1);

Q = sparse(n, n);

current_row = 0;

while ~feof(f)
    if current_row == 0
        current_row = fscanf(f, '%d', 1);
    end
    
    j = fscanf(f, '%d', 1);
    
    if j == 0
        current_row = 0;
    else
        ee = fscanf(f, '%e', 1);
        Q(current_row, j) = ee;
    end

end

e = ones(n, 1);
d = Q * e;
Q = Q + spdiags(-d, 0, n, n);

end

