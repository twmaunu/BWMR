%% Load data

load('../Abalone/abalone.mat');
addpath('../')

T2 = table2array(T(:, 2:end));

Y = T2(:, end);
X = T2(:, 1:end-1);
X = [X, zeros(size(X, 1), 1)];
C = table2array(T(:, 1));
for i=1:length(C)
   if strcmp(C{i}, 'M')
       X(i, end) = 1;
   elseif strcmp(C{i}, 'F')
       X(i, end) = 2;
   else
       X(i, end) = 3;
   end
end

%%


X = X(:, 2:end-1); 

X = X - repmat(mean(X), size(X, 1), 1);
% X = X * sqrtm(inv(cov(X)));
X = [X, ones(size(X, 1), 1)];

Ya = abs(Y - mean(Y));

%% Run methods


[n, d] = size(X);

r = 2;
iter = 4000;

SigmaZ = eye(d);


U0 = 1 * randn(d, r);
U0FR = 1 * randn(d, d);


[B_bwgd,~,~] = bwgd(X, Ya, r, iter, U0, SigmaZ);

[B_bwgdfr,~,~] = bwgd(X, Ya, d, iter, U0FR, SigmaZ);

[B_bwsgd,~,~] = bwsgd(X, Ya, r, 1, U0, SigmaZ);

[B_egdsqrt,~] = egd_sqrt_fr(X, Ya, 1, iter, U0FR, SigmaZ);

[B_egd,~, ~] = egd(X, Ya.^2, r, 0.01, iter, U0, SigmaZ);

[B_egdfr,~, ~] = egd(X, Ya.^2, d, 0.01, iter, U0FR, SigmaZ);


% [M] = sdp_solver(X, Ya);

C1 = 0.5 * ((X .* repmat(Ya.^2, 1, d))' * X / n - sum(Ya.^2) * eye(d)/n);
[u, s, v] = svd(C1);
s = diag(s);
Cspec = u(:, 1:r) * diag(s(1:r)) * u(:, 1:r)';

%% compare barycenter to spectral

figure
imagesc(real(((B_bwgdfr))))
colormap(gray)
colorbar
axis off
title("Barycenter")
set(gca,'FontSize', 18);
set(gca,'FontName', 'Avenir');

figure
plot(svd(B_bwgdfr), 'linewidth', 4)
title("Barycenter")
set(gca,'FontSize', 18);
set(gca,'FontName', 'Avenir');


figure
imagesc(real(((C1))))
colormap(gray)
colorbar
axis off
title("Spectral")

set(gca,'FontSize', 18);
set(gca,'FontName', 'Avenir');

figure
plot(svd(C1), 'linewidth', 4)
title("Spectral")
set(gca,'FontSize', 18);
set(gca,'FontName', 'Avenir');

%% compare barycenter to spectral

figure
imagesc(real(((B_bwgd))))
colormap(gray)
colorbar
axis off
title("Barycenter")
set(gca,'FontSize', 18);
set(gca,'FontName', 'Avenir');

figure
plot(svd(B_bwgd), 'linewidth', 4)
title("Barycenter")
set(gca,'FontSize', 18);
set(gca,'FontName', 'Avenir');


figure
imagesc(real(((Cspec))))
colormap(gray)
colorbar
axis off
title("Spectral")

set(gca,'FontSize', 18);
set(gca,'FontName', 'Avenir');

figure
plot(svd(Cspec), 'linewidth', 4)
title("Spectral")
set(gca,'FontSize', 18);
set(gca,'FontName', 'Avenir');



%% plot Y versus x projected into this space


figure
set(gcf, 'Position',  [100, 100, 900, 600])

tiledlayout(2,3);

[u, s, ~] = svd(B_bwgd);
% s = diag(s);
% Pu = u(:, 1:2) * diag(s(1:2).^-.5);
Pu = -u(:, 1:2);
XP = X * Pu;
nexttile
plot3(XP(:, 1), -XP(:, 2), Y, '.')
xlabel('$P X_1$','interpreter', 'Latex')
ylabel('$P X_2$','interpreter', 'Latex')
% ylabel('xp2')
zlabel('$y$','interpreter', 'latex')
% title('BWGD')
set(gca,'FontSize', 18);
set(gca,'FontName', 'Times');


[u, s, ~] = svd(B_egd);
% s = diag(s);
% Pu = u(:, 1:2) * diag(s(1:2).^-.5);
Pu = u(:, 1:2);
XP = X * Pu;
nexttile
plot3(XP(:, 1), -XP(:, 2), Y, '.')
xlabel('$P X_1$','interpreter', 'Latex')
ylabel('$P X_2$','interpreter', 'Latex')
% ylabel('xp2')
zlabel('$y$','interpreter', 'latex')
% title('GD')
set(gca,'FontSize', 18);
set(gca,'FontName', 'Times');


[u, s, ~] = svd(X' * X / n);
% s = diag(s);
Pu = -u(:, 2:3) ;
XP = X * Pu;
nexttile
plot3(XP(:, 1), -XP(:, 2), Y, '.')
xlabel('$P X_1$','interpreter', 'Latex')
ylabel('$P X_2$','interpreter', 'Latex')
% ylabel('xp2')
zlabel('$y$','interpreter', 'latex')
set(gca,'FontSize', 18);
set(gca,'FontName', 'Times');



[u, s, ~] = svd(B_bwgdfr);
% s = diag(s);
% Pu = u(:, 1:2) * diag(s(1:2).^-.5);
Pu = u(:, 1:2);
XP = X * Pu;
nexttile
plot3(XP(:, 1), -XP(:, 2), Y, '.')
xlabel('$P X_1$','interpreter', 'Latex')
ylabel('$P X_2$','interpreter', 'Latex')
% ylabel('xp2')
zlabel('$y$','interpreter', 'latex')
set(gca,'FontSize', 18);
set(gca,'FontName', 'Times');

[u, s, ~] = svd(B_egdfr);
% s = diag(s);
% Pu = u(:, 1:2) * diag(s(1:2).^-.5);
Pu = u(:, 1:2);
XP = X * Pu;
nexttile
plot3(XP(:, 1), -XP(:, 2), Y, '.')
xlabel('$P X_1$','interpreter', 'Latex')
ylabel('$P X_2$','interpreter', 'Latex')
zlabel('$y$','interpreter', 'latex')
% title('Full-rank BWGD')
set(gca,'FontSize', 18);
set(gca,'FontName', 'Times');

[u, s, ~] = svd(C1);
% s = diag(s);
% Pu = u(:, 1:2) * diag(s(1:2).^-.5);
Pu = u(:, 1:2)*diag([1,1]);
XP = X * Pu;
nexttile
plot3(XP(:, 1), -XP(:, 2), Y, '.')
xlabel('$P X_1$','interpreter', 'Latex')
ylabel('$P X_2$','interpreter', 'Latex')
% ylabel('xp2')
zlabel('$y$','interpreter', 'latex')
% title('Spectral')
set(gca,'FontSize', 18);
set(gca,'FontName', 'Times');

f = gcf;
exportgraphics(f,'aba_proj.png','Resolution',300)



