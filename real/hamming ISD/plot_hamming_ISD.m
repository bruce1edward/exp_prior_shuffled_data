load('hamming_ISD.mat')
MSE_beta = reshape(mean(error_est,1),[numel(K_n),7]);
std_error = reshape(std(error_est,0,1),[numel(K_n),7]);

figure
hold on 
p1 = errorbar(K_n, MSE_beta(:,1), 3*std_error(:,1)/10, 'k-o', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'k');
p2 = errorbar(K_n, MSE_beta(:,2), 3*std_error(:,2)/10, 'm-^', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'm');
p4 = errorbar(K_n, MSE_beta(:,4), 3*std_error(:,4)/10, 'b-d', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'b');
p5 = errorbar(K_n, MSE_beta(:,5), 3*std_error(:,5)/10, 'r-*', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'r');
xlabel('k/n', 'FontSize', 14)
ylabel('Relative Estimation Errors', 'FontSize', 14)
ylim([0 0.8])
yticks([0 : 0.2 : 0.8])
xlim([min(K_n) max(K_n)])
xticks(K_n)
lh = legend([p1 p2 p4 p5], {'naive', 'robust', 'EM', 'EMH'}, 'location', [0.2488 0.7814 0.3500 0.1250]);
lh.NumColumns = 2;
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
