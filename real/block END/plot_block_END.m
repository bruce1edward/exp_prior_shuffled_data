load('block_END.mat')
MSE_beta = reshape(mean(error_est,1),[7,5]);
std_error = reshape(std(error_est,0,1),[7,5]);
K_n = 0.1:0.05:0.4;

figure
hold on 
p1 = errorbar(K_n, MSE_beta(:,1), 3*std_error(:,1)/10,'k-*', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'k');
p2 = errorbar(K_n, MSE_beta(:,2), 3*std_error(:,2)/10,'b-s', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'b');
p3 = errorbar(K_n, MSE_beta(:,3), 3*std_error(:,3)/10, 'r-o', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'r');
p4 = errorbar(K_n, MSE_beta(:,4), 3*std_error(:,4)/10, 'm-d', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'm');
p5 = errorbar(K_n, MSE_beta(:,5), 3*std_error(:,5)/10, 'c-+', 'LineWidth', 1.5, 'MarkerSize',11, 'MarkerFaceColor', 'c');
xlabel('k*B/n', 'FontSize', 14)
%ylabel('Relative Estimation Errors', 'FontSize', 14)
%yticks([0 : 0.25 :2])
ylim([0 2.5])
xlim([0.10 max(K_n)])
xticks(K_n)
lh = legend([p1 p2 p3 p4 p5], {'naive', 'EM', 'EMM', 'LL','C'}, 'location', [0.3696 0.6856 0.1696 0.1750]);
lh.NumColumns = 2;
lh.FontSize = 14;
lh.ItemTokenSize = [30 70];  
legend('boxoff')
ax=gca;
ax.FontSize = 14;
grid on
set(gcf, 'Color', 'w');
fig = gcf;
fig.Units               = 'centimeters';
fig.Position(3)         = 8;
fig.Position(4)         = 7;
hold off
