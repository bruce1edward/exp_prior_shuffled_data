load('block_BSD.mat');
error_est1 = zeros(8,100,5);
dev_est1 = zeros(8,100,5);
for i = 1:8
    error_est1(i,:,:) = error_est(1+(i-1)*100:100+(i-1)*100,:);
    dev_est1(i,:,:) = dev_est(1+(i-1)*100:100+(i-1)*100,:);
end

MSE_beta = reshape(mean(error_est1,2),[8,5]);
R_sq = reshape(mean(dev_est1,2),[8,5]);
K_n = 0.05:0.05:0.4;

std_error = reshape(std(error_est1,0,2),[8,5]);
std_error1 = reshape(std(dev_est1,0,2),[8,5]);

figure
hold on 
p1 = errorbar(K_n, MSE_beta(:,1), 3*std_error(:,1)./sqrt(100) ,'k-*', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'k');
p2 = errorbar(K_n, MSE_beta(:,2), 3*std_error(:,2)./sqrt(100) ,'b-s', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'b');
p3 = errorbar(K_n, MSE_beta(:,3), 3*std_error(:,3)./sqrt(100) ,'r-o', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'r');
p4 = errorbar(K_n, MSE_beta(:,4), 3*std_error(:,4)./sqrt(100), 'm-d', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'm');
p5 = errorbar(K_n, MSE_beta(:,5), 3*std_error(:,5)./sqrt(100), 'c-+', 'LineWidth', 1.5, 'MarkerSize',11, 'MarkerFaceColor', 'c');
xlabel('k*B/n', 'FontSize', 14)
%ylabel('Relative Estimation Errors', 'FontSize', 14)
xlim([0.1 max(K_n)])
xticks(K_n)
ylim([0 0.6])
lh = legend([p1 p2 p3 p4 p5], {'naive', 'EM', 'EMM', 'LL', 'C'}, 'location', [0.3738 0.7382 0.3000 0.1833]);
lh.NumColumns = 2;
lh.FontSize = 14;
lh.ItemTokenSize = [30 70];  
legend('boxoff')
pos_leg = get(lh,'Position');
ax=gca;
ax.FontSize = 14;
grid on
set(gcf, 'Color', 'w');
fig = gcf;
fig.Units               = 'centimeters';
fig.Position(3)         = 8;
fig.Position(4)         = 7;
hold off

