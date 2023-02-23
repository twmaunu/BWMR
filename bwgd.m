function [barys, dists, bary, times] = bwgd(X, Y, r, iter, U0, Uz)
%BWGD

[n, d] = size(X);
CX = X' * X / n;
CXsqrt = sqrtm(CX);
CXsqrtinv = inv(CXsqrt);
X2 = X * CXsqrtinv;

bary = U0;
dists = zeros(iter, 1);
YX = repmat(Y, 1, d) .* X2;

times = zeros(iter,1);
tic
for i=1:iter
    
    step = 1;
    
    denoms = vecnorm(YX * bary, 2, 2);
%     keyboard
    bary = (step) * YX' * ((YX * bary) ./ repmat(denoms, 1, r)) * (1/n) + (1-step) * bary;
    
    Ub = CXsqrtinv * bary;
    %keyboard
    if size(Ub, 2) ~= size(Uz, 2)
         dists(i) = norm(sqrtm(Ub * Ub') - sqrtm(Uz * Uz'), 'fro');
    else
        
        [u, ~, v] = svd(Uz' * Ub);
         dists(i) = norm(u * v' * Ub'  - Uz', 'fro');
    end
    
    
    times(i) = toc;
%     i
end


barys = CXsqrtinv * bary * bary' * CXsqrtinv;
bary = CXsqrtinv * bary;

end

