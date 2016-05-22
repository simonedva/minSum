function vHat = minSumOptimized(gamma, P, d, iteration) 
%init
z = P*gamma;

for k=1:iteration
    
    y = minSumOperation(z,d);
    u = gamma + P'*y;
    
    %decisione
    vHat = u<0;
    
    z = P*u - y;
end