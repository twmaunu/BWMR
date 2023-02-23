function [barys, dists] = egd_sqrt_fr(X, Y, step, iter, U0, SigmaZ)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[n, d] = size(X);
CX = X' * X / n;
CXsqrt = sqrtm(CX);
CXsqrtinv = inv(CXsqrt);
X2 = X * CXsqrtinv;

YX = repmat(Y, 1, d) .* X2;


bary = U0*U0';

dists = zeros(iter, 1);

stepi = step;
for i=1:iter
    
    
    denoms = sum((YX * bary).*YX, 2).^.5;
    T = YX' * ((YX) ./ repmat(denoms, 1, d))/n ;
    
%     while cost(X2, Y,  bary - (stepi) * (eye(d) - T)) > cost(X2, Y, bary)
%         stepi = stepi / 2;
% %         stepi
%     end

    if mod(i, 200)==0
        stepi = stepi/2;
    end
    
    bary = bary - (stepi) * (eye(d) - (T));
    [U, D] = eig(bary);
    D = diag(D);
    D(D < 0) = 0.0;
    bary = U * diag(D) * U';
    
    
    
    
    
    

    [ub, sb, ~] = svd(CXsqrtinv * bary * CXsqrtinv);
    [uz, sz, ~] = svd(SigmaZ);
    
    dists(i) = norm(ub * (sb.^.5) * ub' - uz * (sz.^.5) * uz', 'fro');

end

barys = bary;

end


function [c] = cost(X, Y, A)

Y2 = sum((X * A).*X, 2).^.5;

c = norm(Y - Y2);

end
