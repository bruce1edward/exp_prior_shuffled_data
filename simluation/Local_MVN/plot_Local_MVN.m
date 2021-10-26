load('Local_MVN.mat');
MSE_beta = reshape(mean(error_est,2),[numel(r_n),nmethod]);
std_beta = reshape(std(error_est,0,2),[numel(r_n),nmethod]);

figure
hold on 
p1 =  errorbar(r_n, MSE_beta(:,1), 3*std_beta(:,1)/sqrt(iters), 'k-*', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'k');
p2 = errorbar(r_n, MSE_beta(:,2), 3*std_beta(:,2)/sqrt(iters), 'g-^', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'g');
p3 = errorbar(r_n, MSE_beta(:,3), 3*std_beta(:,3)/sqrt(iters), 'b-s', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'b');
p4 = errorbar(r_n, MSE_beta(:,4), 3*std_beta(:,4)/sqrt(iters), 'r-o', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'r');
p5 = errorbar(r_n, MSE_beta(:,5), 3*std_beta(:,5)/sqrt(iters), 'm-+', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'b');
xlabel('r', 'FontSize', 14)
%ylabel('Relative Estiamtion error', 'FontSize', 14)
ax=gca;
ax.FontSize = 14;
xlim([min(r_n) max(r_n)])
xticks(r_n)
ylim([0 1.4])
yticks(0:0.2:1.4)
lh = legend([p1 p2 p3 p4 p5], {'naive','oracle','EM','EML','averaging'}, 'location', [0.416 0.7282 0.3661 0.1833]);
lh.NumColumns = 2;
lh.FontSize = 14;
lh.ItemTokenSize = [30 70];  
pos = get(lh,'Position');
legend('boxoff')
grid on
set(gcf, 'Color', 'w');
fig = gcf;
fig.Units               = 'centimeters';
fig.Position(3)         = 8;
fig.Position(4)         = 7;
hold off


