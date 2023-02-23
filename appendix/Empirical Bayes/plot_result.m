MSE_beta = reshape(mean(error_est,2),[numel(k_n),nmethod]);
R_sq = reshape(mean(r_square,2),[numel(k_n),nmethod]);
figure
hold on 
p1 =  plot(k_n, MSE_beta(:,1), 'k-*', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'k');
p2 = plot(k_n, MSE_beta(:,2), 'g-^', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'g');
p3 = plot(k_n, MSE_beta(:,3), 'b-s', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'b');
p4 = plot(k_n, MSE_beta(:,4), 'r-o', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'r');
%p6 = plot(k_n, MSE_beta(:,6), 'c-x', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'r');
xlabel('k/n', 'FontSize', 14)
ylabel('Relative Estiamtion error', 'FontSize', 14)
%yticks([0 : 0.25 :2])
xlim([min(k_n) max(k_n)])
xticks(k_n)
%title([' \beta_{0} = ' num2str(b0) ' , ' ' ||\beta^*||_{2} = ' num2str(b)])
lh = legend([p1 p2 p3 p4], {'Naive','EM','EMH','EMEB'}, 'location', 'best');
lh.NumColumns = 2;
lh.FontSize = 10;
lh.ItemTokenSize = [30 70];  
legend('boxoff')
grid on
set(gcf, 'Color', 'w');
fig = gcf;
fig.Units               = 'centimeters';
fig.Position(3)         = 8;
fig.Position(4)         = 7;
hold off

%export_fig('plot/lm_sparse_ree.pdf')
%export_fig('plot/glm_sparse_ree.pdf')
export_fig('lm_EB_sparse_ree.pdf')

figure
hold on 
p1 =  plot(k_n, R_sq(:,1), 'k-*', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'k');
p2 = plot(k_n, R_sq(:,2), 'g-^', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'g');
p3 = plot(k_n, R_sq(:,3), 'b-s', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'b');
p4 = plot(k_n, R_sq(:,4), 'r-o', 'LineWidth', 1.5, 'MarkerSize',9, 'MarkerFaceColor', 'r');
xlabel('k/n', 'FontSize', 14)
ylabel('R^2', 'FontSize', 14)
%yticks([0 : 0.25 :2])
xlim([min(k_n) max(k_n)])
xticks(k_n)
%title([' \beta_{0} = ' num2str(b0) ' , ' ' ||\beta^*||_{2} = ' num2str(b)])
lh = legend([p1 p2 p3 p4 p5], {'Naive','EM','EMH','EMEB'}, 'location', 'best');
lh.NumColumns = 1;
lh.FontSize = 10;
lh.ItemTokenSize = [30 70];  
legend('boxoff')
grid on
set(gcf, 'Color', 'w');
fig = gcf;
fig.Units               = 'centimeters';
fig.Position(3)         = 8;
fig.Position(4)         = 7;
hold off

export_fig('lm_EB_sparse_r_2.pdf')