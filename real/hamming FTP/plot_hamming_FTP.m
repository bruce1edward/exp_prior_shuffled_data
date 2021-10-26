load('hamming_FTP.mat')
k_n = 0.1 : 0.05 : 0.4; nmethod = 4;iters = 100;
MSE_beta = reshape(mean(error_est,1),[numel(k_n),nmethod]);
std_beta = reshape(std(error_est,0,1),[numel(k_n),nmethod]);

figure
hold on 
p1 =  errorbar(k_n, MSE_beta(:,1), 3*std_beta(:,1)/sqrt(iters), 'k-*', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'k');
p2 = errorbar(k_n, MSE_beta(:,2), 3*std_beta(:,2)/sqrt(iters), 'b-^', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'b');
p3 = errorbar(k_n, MSE_beta(:,3), 3*std_beta(:,3)/sqrt(iters), 'r-s', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'r');
xlabel('k/n', 'FontSize', 14)
%ylabel('Relative Estiamtion error', 'FontSize', 14)
%yticks([0 : 0.25 :2])
ax=gca;
ax.FontSize = 14;
xlim([min(k_n) max(k_n)])
xticks(k_n)
lh = legend([p1 p2 p3], {'Naive','EM','EMH'}, 'location', [0.1708  0.7022 0.1768 0.1833]);
lh.NumColumns = 1;
lh.FontSize = 14;
lh.ItemTokenSize = [30 70];  
legend('boxoff')
pos_leg = get(lh,'Position');
grid on
set(gcf, 'Color', 'w');
fig = gcf;
fig.Units               = 'centimeters';
fig.Position(3)         = 8;
fig.Position(4)         = 7;
hold off

