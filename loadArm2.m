function Q = loadArm2(armfile)
%LOADARM2 Loads a Q matrix from the ARM format of Moebius. 
%
% Q = LOADARM2(ARMFILE) loads the infinitesimal generator Q of a Markov
% chains stored in the .arm format used by Moebius. Therefore, you can use
% this function to interface FUNM_MARKOV and FUNM_MARKOV_SENSITIVITY with
% models obtained working in that environment. 
%
% This function is weakly tuned for performances, compared to LOADARM. 
%
% Author: Leonardo Robol <leonardo.robol@isti.cnr.it>

f = dlmread(armfile); 
n = f(2);
current_row = 0;

k = 3;
D = zeros(length(f) - 2, 3);
kk = 1; 

while k <= length(f)
    if current_row == 0
        current_row = f(k);
        k = k + 1;
    end
    
    j = f(k); 
    k = k + 1;
    
    if j == 0
        current_row = 0;
    else
        D(kk, 1) = current_row; 
        D(kk, 2) = j; 
        D(kk, 3) = f(k);
        kk = kk + 1;
        k = k + 1; 
    end
end

Q = spconvert(D(1:kk-1,:));

if min(size(Q)) < n
    Q(n,n) = 0;
end

e = ones(n, 1);
d = Q * e;
Q = Q + spdiags(-d, 0, n, n);

end

