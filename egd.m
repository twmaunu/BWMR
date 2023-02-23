function [A, dists, Ai, times] = egd(X, Y, r, step, iter, U0, Uz, ls)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[n, d] = size(X);

% YX = repmat(Y.^.5, 1, d) .* X; 
% 
% 
% Z = YX' * YX /(2*n) ;
% [uz, sz, vz] = svd(Z);
% 
% sz2 = diag(sz(1:r, 1:r).^.5);
% l = sum(Y)/(2*n);
% sz2 = diag(sz2 - l);
% 
% Ai = uz(:, 1:r) * sz2;


% Ai = randn(d, r)/sqrt(d);
Ai = U0;

dists = zeros(iter, 1);

stepi = step;
times = zeros(iter, 1);
tic;
for i=1:iter
    %keyboard
    
    norms2 = vecnorm(X * Ai, 2, 2).^2;
    
    Xn = X .* repmat(norms2-Y, 1, d);
    
%     cost(X, Y, Ai)
if ls ==1
    while cost(X, Y, Ai - (stepi/n) * Xn' * (X * Ai)) > cost(X, Y, Ai)
        stepi = stepi / 2;
        %stepi
    end
        Ai = Ai - (stepi/n) * Xn' * (X * Ai);

else
    
    Ai = Ai - (step/n) * Xn' * (X * Ai);
    
end
    %Ai
%     [ub, sb, ~] = svd(Ai * Ai');
    
    [u, ~, v] = svd(Uz' * Ai);
    dists(i) = norm(u * v' * Ai'  - Uz', 'fro');
    
%         dists(i) = norm(Ai * Ai' - SigmaZ, 'fro')/norm(SigmaZ, 'fro');

%     if mod(i, 100)==0
%         i
%     end
% i
times(i) = toc;
end


A = Ai * Ai';

end


function [c] = cost(X, Y, A)

Y2 = sum(((X * A) * A').*X, 2);

c = norm(Y - Y2);

end
