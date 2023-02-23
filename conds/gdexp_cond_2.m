%% run experiment
clc; clear

addpath("../")

d = 32;
beta = 1;

% rs = [1];
rs = [2, 4, 16];

%facts = [5];
% facts = [3, 5, 10];

scales = [0, 1, 2];

nr = length(rs);
nf = length(scales);

% n = 3 * d * d * r;
%n = 100;

reps = 20;



BGD = cell(nr, nf);
BFRGD = cell(nr, nf);
EGD = cell(nr, nf);
ESQRT = cell(nr, nf);
BSGD = cell(nr, nf);
SDP = cell(nr, nf);
SPEC = cell(nr, nf);


TBGD = cell(nr, nf);
TFRBGD = cell(nr, nf);
TEGD = cell(nr, nf);
TESQRT = cell(nr, nf);
TBSGD = cell(nr, nf);
TSDP = cell(nr, nf);
TSPEC = cell(nr, nf);

for i=1:length(rs)
    
    for j=1:length(scales)
        
        
        
        r = rs(i);
        f = scales(j);
        n = 5 * d * r;
        
        
        iter = 500;
        
        distBGDreps = zeros(reps, iter);
        distBFRGDreps = zeros(reps, iter);
        distBSGDreps = zeros(reps, n);
        distEGDreps = zeros(reps, iter);
        distESQRTreps = zeros(reps, iter);
        distSDPreps = zeros(reps, 1);

        timesBGD = zeros(reps, 1);
        timesBFRGD = zeros(reps, 1);
        timesEGD = zeros(reps, 1);
        timesESQRT = zeros(reps, 1);
        timesSDP = zeros(reps, 1);
        timesBSGD = zeros(reps, 1);



        for k=1:reps
            
            UZ = randn(d, r);
            s1 = diag((r:-1:1).^f);
            SigmaZ = UZ  * s1 * UZ';
            
            norm(UZ);
            
            SigmaX = eye(d);
            X = randn(n, d) * SigmaX.^.5;
            Xn = vecnorm(X, 2, 2);
            Y = (sum((X * SigmaZ).*X, 2)).^.5;
            
            U0 = 1 * randn(d, r);
            U0FR = 1 * randn(d, d);
            
            
            tic;
            [B,distsB,B2] = bwgd(X, Y, r, iter, U0, UZ * s1.^.5);
            timesBGD(k) = toc;
            distsB;
            
            tic;
            [B,distsBFR,B2] = bwgd(X, Y, d, iter, U0FR, UZ * s1.^.5);
            timesBFRGD(k) = toc;
            distsBFR;
            
            
            
            tic;
            [B,distsBS,B2] = bwsgd(X, Y, r, 1, U0, SigmaZ);
            timesBSGD(k) = toc;
            distsBS;
            
            tic;
            [B,distsESQRT] = egd_sqrt_fr(X, Y, 1, iter, U0FR, SigmaZ);
            timesESQRT(k) = toc;
            distsBFR;
            
            tic;
%             [E,distsE, E2] = egd(X, Y.^2, r, .1*.1^f/(d*r), iter, U0, UZ * s1.^.5, 0);
            [E,distsE, E2] = egd(X, Y.^2, r, 1, iter, U0, UZ * s1.^.5, 1);
            timesEGD(k) = toc;
            
%             tic;
%             [M] = sdp_solver(X, Y);
%             [ub, sb, ~] = svd(M);
%             [uz, sz, ~] = svd(SigmaZ);
%             distsSDP = norm(ub * (sb.^.5) * ub' - uz * (sz.^.5) * uz', 'fro');
%             timesSDP(k) = toc;
            distsSDP = 1;
            timesSDP = 0;
            
            tic;
            C1 = 0.5 * ((X .* repmat(Y.^2, 1, d))' * X / n - sum(Y.^2) * eye(d)/n);
            [u, s, v] = svd(C1);
            s = diag(s);
            Cspec = u(:, 1:r) * diag(s(1:r)) * u(:, 1:r)';
            
            [ub, sb, ~] = svd(Cspec);
            [uz, sz, ~] = svd(SigmaZ);
            distsSPEC = norm(ub * (sb.^.5) * ub' - uz * (sz.^.5) * uz', 'fro');
            timesSPEC(k) = toc;
            
            
            distBGDreps(k, :) = distsB;
            distBFRGDreps(k, :) = distsBFR;
            distBSGDreps(k, :) = distsBS;
            distEGDreps(k, :) = distsE;
            distESQRTreps(k, :) = distsESQRT;
            
            
            distSDPreps(k) = distsSDP;
            distSPECreps(k) = distsSPEC;
            
            k
        end
        
        BGD{i, j} = distBGDreps;
        BFRGD{i, j} = distBFRGDreps;
        BSGD{i, j} = distBSGDreps;
        EGD{i, j} = distEGDreps;
        ESQRT{i, j} = distESQRTreps;
        SDP{i, j} = distSDPreps;
        SPEC{i, j} = distSPECreps;
        
        TBGD{i, j} = timesBGD;
        TFRBGD{i, j} = timesBFRGD;
        TEGD{i, j} = timesEGD;
        TESQRT{i, j} = timesESQRT;
        TBSGD{i, j} = timesEGD;
        TSDP{i, j} = timesSDP;
        TSPEC{i, j} = timesSPEC;
        
        [i,j]
    end
    
end


% sum(distBGDreps(:, end) <1e-1)/reps
% sum(distBSGDreps(:, end) <1e-1)/reps
% sum(distEGDreps(:, end) <1e-1)/reps
% sum(distSDPreps(:, end) <1e-1)/reps


% figure
% plot(costs)

save('exp_cond_2')

%%
load('exp_cond_2')

figure
set(gcf, 'Position',  [100, 100, 900, 600])

tiledlayout(nr,nf);

for i=1:(nr)
    
    for j=1:nf
        nexttile
        
        if i==1 && j==1
%             b1 = BGD{i, j};
%             e1 = EGD{i, j};
%             bf1 = BFRGD{i, j};
%             es1 = ESQRT{i, j};
%             spec = SPEC{i, j};
%             
            pb = plot(log10(mean(BGD{i, j})), '-', 'color', [0, 0.4470, 0.7410], 'linewidth', 4);
            hold on
            pe = plot(log10(mean(EGD{i, j})), '-', 'color', [0.8500, 0.6250, 0.1980], 'linewidth', 4);
            pfr = plot(log10(mean(BFRGD{i, j})), '-', 'color', [0.4940, 0.1840, 0.5560], 'linewidth', 4);
            psp = plot(log10(mean(spec)), '-', 'color', [.1, .6, .1], 'linewidth', 2);
%             esq = plot(log10(es1(1,:)), '-', 'color', [0.4660 0.6740 0.1880, 2/3], 'linewidth', 2);
            
            set(gca,'FontSize', 18);
            set(gca,'FontName', 'Times');
            grid on
        end
        
        plot(log10(mean(BGD{i, j})'), '-', 'color', [0, 0.4470, 0.7410], 'linewidth', 4)
        hold on
        plot(log10(mean(EGD{i, j})'), '-', 'color', [0.8500, 0.6250, 0.1980], 'linewidth', 4)
        plot(log10(mean(BFRGD{i, j})'), '-', 'color', [0.4940, 0.1840, 0.5560], 'linewidth', 4)
        
        plot(repmat([1;iter], 1, reps), log10([mean(SPEC{i,j}); mean(SPEC{i,j})]),'color', [.1, .6, .1], 'linewidth', 4)
%         plot(log10(ESQRT{i, j}'), '-', 'color', [0.4660 0.6740 0.1880, 1/3], 'linewidth', 1.3)
        set(gca,'FontSize', 18);
        set(gca,'FontName', 'Times');
        grid on
        
        if i==(nr) && j==2
            xlabel('Iteration $k$','Interpreter','latex')
        end
        if j==1 && i==2
            ylabel('$\log \|\Sigma_k^{1/2} - S^{1/2} \|$','Interpreter','latex')
        end
        if i==1 
            xlim([0, 500])
            ylim([-12, 2])
        end
        if i==2 
             xlim([0, 500])
            ylim([-12, 2])
        end
        if i==3 
             xlim([0, 500])
            ylim([-2, 2])
        end
        
    end
end

% cb = legend([pb, pe, pfr, esq], 'BWGD', 'EGD', 'Full-rank BWGD', 'EGD SQRT');
cb = legend([pb, pfr, pe, psp], 'BWGD', 'Full-rank BWGD', 'GD', 'Spectral');
cb.Layout.Tile = 'east';
f = gcf;
exportgraphics(f,'exp_scales_e2.png','Resolution',300)


%%

% figure
% set(gcf, 'Position',  [100, 100, 900, 600])
% 
% t = tiledlayout(nr,nf,'TileSpacing','Compact');
% 
% for i=1:(nr)
%     
%     for j=1:nf
%         nexttile
%         
%         if i==1 && j==1
%             b1 = BGD{i, j};
%             e1 = EGD{i, j};
%             bf1 = BFRGD{i, j};
%             es1 = ESQRT{i, j};
%             spec = SPEC{i, j};
%             
%             pb = plot(log10(b1(1,:)), '-', 'color', [0, 0.4470, 0.7410, 2/3], 'linewidth', 2);
%             hold on
%             pe = plot(log10(e1(1,:)), '-', 'color', [0.8500, 0.3250, 0.0980, 2/3], 'linewidth', 2);
%             pfr = plot(log10(bf1(1,:)), '-', 'color', [0.4940, 0.1840, 0.5560, 2/3], 'linewidth', 2);
%             psp = plot(log10(spec), '-', 'color', [1, 0, 0, 2/3], 'linewidth', 2);
% %             esq = plot(log10(es1(1,:)), '-', 'color', [0.4660 0.6740 0.1880, 2/3], 'linewidth', 2);
%             
%             set(gca,'FontSize', 18);
%             set(gca,'FontName', 'Times');
%             grid on
%         end
%         
%         plot(log10(BGD{i, j}'), '-', 'color', [0, 0.4470, 0.7410, 1/3], 'linewidth', 1.3)
%         hold on
%         plot(log10(EGD{i, j}'), '-', 'color', [0.8500, 0.3250, 0.0980, 1/3], 'linewidth', 1.3)
%         plot(log10(BFRGD{i, j}'), '-', 'color', [0.4940, 0.1840, 0.5560, 1/3], 'linewidth', 1.3)
%         
%         plot(repmat([1;iter], 1, reps), log10([SPEC{i,j}; SPEC{i,j}]),'color', [1, 0, 0, 1/6], 'linewidth', 1.3)
% %         plot(log10(ESQRT{i, j}'), '-', 'color', [0.4660 0.6740 0.1880, 1/3], 'linewidth', 1.3)
%         set(gca,'FontSize', 18);
%         set(gca,'FontName', 'Times');
%         grid on
%         
%         if i==(nr) && j==2
% %             xlabel('Iteration $k$','Interpreter','latex')
%         end
%         if j==1 && i==2
%         end
%         ylim([-8, 1.5])
%         
%     end
% end
% 
% xlabel(t,'Iteration $k$','Interpreter','latex','FontSize', 18)
% ylabel(t,'$\log \|\Sigma_k^{1/2} - S^{1/2} \|$','Interpreter','latex','FontSize', 18)
% set(gca,'FontSize', 18);
%         set(gca,'FontName', 'Times');
% % cb = legend([pb, pe, pfr, esq], 'BWGD', 'EGD', 'Full-rank BWGD', 'EGD SQRT');
% cb = legend([pb, pfr, pe, psp], 'BWGD', 'Full-rank BWGD', 'GD', 'Spectral');
% cb.Layout.Tile = 'east';
% f = gcf;
% exportgraphics(f,'exp_scales_e.png','Resolution',300)

%%

% errs_b = zeros(nr, nf);
% errs_e = zeros(nr, nf);
% errs_bsgd = zeros(nr, nf);
% errs_sdp = zeros(nr, nf);
% 
% for i=1:nr
%     
%     for j=1:nf
%         tmp = BGD{i, j};
%         errs_b(i, j) = log10(mean(tmp(:, end)));
%         
%         tmp = EGD{i, j};
%         errs_e(i, j) = log10(mean(tmp(:, end)));
%         
% %         tmp = BSGD{i, j};
% %         errs_bsgd(i, j) = log10(mean(tmp(:, end)));
%         
%         tmp = SPEC{i, j};
%         errs_spec(i, j) = log10(mean(tmp(:, end)));
%         
%         
%         tmp = SDP{i, j};        
%         errs_sdp(i, j) = log10(mean(tmp(:, end)));
%     end
%     
% end
% 
% m1 = min(min(errs_b));
% m2 = min(min(errs_e));
% m3 = min(min(errs_spec));
% m4 = min(min(errs_sdp));
% 
% M1 = max(max(errs_b));
% M2 = max(max(errs_e));
% M3 = max(max(errs_spec));
% M4 = max(max(errs_sdp));
% 
% c1 = min([m1, m2, m3, m4]);
% C1 = max([M1, M2, M3, M4]);
% 
% figure
% 
% tiledlayout(2,2);
% nexttile
% imagesc(errs_b)
% caxis([c1, C1]);
% set(gca,'FontSize', 18);
% set(gca,'FontName', 'Times');
% xticklabels(facts);
% yticklabels(rs);
% title('(a)')
% 
% nexttile
% imagesc(errs_e)
% caxis([c1, C1]);
% set(gca,'FontSize', 18);
% set(gca,'FontName', 'Times');
% xticklabels(facts);
% yticklabels(rs);
% % colorbar
% title('(b)')
% 
% nexttile
% imagesc(errs_spec)
% caxis([c1, C1]);
% set(gca,'FontSize', 18);
% set(gca,'FontName', 'Times');
% xticklabels(facts);
% yticklabels(rs);
% title('(c)')
% 
% % nexttile
% % imagesc(errs_sdp)
% % caxis([c1, C1]);
% % set(gca,'FontSize', 18);
% % set(gca,'FontName', 'Times');
% % xticklabels(facts);
% % yticklabels(rs);
% 
% title('(d)')
% 
% cb = colorbar;
% cb.Layout.Tile = 'east';




