function [barys, dists, bary] = bwsgd(X, Y, r, s0, bary0, SigmaZ)
%BWSGD 

[n, d] = size(X);
YX = repmat(Y, 1, d) .* X;
bary = bary0;
dists = zeros(n, 1);

for i=1:n
    [ub, sb, ~] = svd(bary * bary');
    [uz, sz, ~] = svd(SigmaZ);
    dists(i) = norm(ub * (sb.^.5) * ub' - uz * (sz.^.5) * uz', 'fro');

    T = T_samp_LR(bary, YX(i,:)');
    step = s0/(s0 + i );
    bary = ((1-step) * eye(d) + step * T ) * bary;
    
    
end

barys = bary * bary';

end

