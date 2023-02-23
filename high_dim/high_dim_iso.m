
addpath('../')

ds = [512];
r = 8;
f = 3;
iter = 2000;
itere = 600;
reps = 1;

conv_b = zeros(iter, length(ds), reps);
conv_e1 = zeros(itere, length(ds), reps);
conv_e2 = zeros(itere, length(ds), reps);
conv_e3 = zeros(itere, length(ds), reps);
conv_e4 = zeros(itere, length(ds), reps);

time_b = zeros(iter, length(ds), reps);
time_e1 = zeros(itere, length(ds), reps);
time_e2 = zeros(itere, length(ds), reps);
time_e3 = zeros(itere, length(ds), reps);
time_e4 = zeros(itere, length(ds), reps);

for k=1:reps
for i=1:length(ds)
    d = ds(i);
    beta = 1;
    
    % rs = [1];
  
   
    n = f * d * r;
    
    
    UZ = randn(d, r);
    scales = eye(r);
%     scales = diag(r:-1:1)/r;
    SigmaZ = UZ  * scales * UZ';
    [u, s, v] = svd(SigmaZ);
    UZ = UZ * scales.^.5;
    
    SigmaX = eye(d);
    X = randn(n, d) * SigmaX.^.5;
    Xn = vecnorm(X, 2, 2);
    Y = (sum((X * SigmaZ).*X, 2)).^.5;
    
    U0 = 1 * randn(d, r);
    U0FR = 1 * randn(d, d);
    
    
    [B,distsB,B2, timesb] = bwgd(X, Y, r, iter, U0, UZ);
    
    [E,distsE1, E2, timese1] = egd(X, Y.^2, r, .01/d, itere, U0, UZ,0);
    [E,distsE2, E2, timese2] = egd(X, Y.^2, r, .05/d, itere, U0, UZ,0);
    [E,distsE3, E2, timese3] = egd(X, Y.^2, r, .1/d, itere, U0, UZ,0);
    [E,distsE4, E2, timese4] = egd(X, Y.^2, r, 10/d, itere, U0, UZ,1);
    
    conv_b(:,i,k) = distsB;
    conv_e1(:,i,k) = distsE1;
    conv_e2(:,i,k) = distsE2;
    conv_e3(:,i,k) = distsE3;
    conv_e4(:,i,k) = distsE4;
    
    time_b(:,i,k) = timesb;
    time_e1(:,i,k) = timese1;
    time_e2(:,i,k) = timese2;
    time_e3(:,i,k) = timese3;
    time_e4(:,i,k) = timese4;

    
    [k, i]
end
    
end

save('high_dim_iso')

%%


figure, hold on
cmap = [0 0.1470 0.7410; 0.0 0.4470 0.7410; 0.0 0.7470 0.7410];

cmap2 = [0.8500 0.0250 0.0980; 0.8500 0.3250 0.0980; 0.8500 0.6250 0.0980];


dcord = [0 0.1470 0.7410;...
    0.8500 0.0250 0.0980;...
    0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;...
    0.4660 0.6740 0.1880];

linestyles = ['-', ':'];
for i=1:length(ds)
    plot(mean(time_b(:, i,:), 3), log10(mean(conv_b(:, i, :), 3)), '-', 'color', dcord(1,:), 'linewidth', 2)
    plot(mean(time_e1(:, i,:), 3), log10(mean(conv_e1(:, i, :), 3)), '-', 'color', dcord(2,:), 'linewidth', 2)
    plot(mean(time_e2(:, i,:), 3), log10(mean(conv_e2(:, i, :), 3)), '--', 'color', dcord(2,:), 'linewidth', 2)
%     plot(mean(time_e3(:, i,:), 3), log10(mean(conv_e3(:, i, :), 3)), ':', 'color', dcord(2,:), 'linewidth', 2)
    plot(mean(time_e4(:, i,:), 3), log10(mean(conv_e4(:, i, :), 3)), '-.', 'color', dcord(2,:), 'linewidth', 2)
%         plot(log10(mean(conv_b(:, i, :), 3)), '-', 'color', dcord(1,:), 'linewidth', 2)
%     plot(log10(mean(conv_e1(:, i, :), 3)), '-', 'color', dcord(2,:), 'linewidth', 2)
%     plot(log10(mean(conv_e2(:, i, :), 3)), '--', 'color', dcord(2,:), 'linewidth', 2)
% %     plot(mean(time_e3(:, i,:), 3), log10(mean(conv_e3(:, i, :), 3)), ':', 'color', dcord(2,:), 'linewidth', 2)
%     plot(log10(mean(conv_e4(:, i, :), 3)), '-.', 'color', dcord(2,:), 'linewidth', 2)
end
xlim([0,20])
grid on
xlabel('Time (s)', 'interpreter', 'latex')
ylabel('$\log_{10} \|\Sigma_k^{1/2} - S^{1/2} \|$','Interpreter','latex')
set(gca,'FontSize', 18);
legend('BWGD', 'GD, $\eta = 0.01$', 'GD, $\eta =0.1$', 'GD with linesearch', 'interpreter', 'latex', 'location', 'southwest')

exportgraphics(gcf,'exp_highdim.png','Resolution',300)