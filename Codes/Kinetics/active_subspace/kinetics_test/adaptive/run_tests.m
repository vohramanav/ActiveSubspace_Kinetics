close all;
clear all;

% gradient based
active_subspace;

xi_grad = xi;
g_grad = G;

% gradient free
local_linear_approx;

xi_loclin = xi;
g_loclin = f;

% compare SSPs
%figure;
%hold on;
%g1 = V(:,1)'*xi_grad';
%g2 = W(:,1)'*xi_loclin';
%plot(g1, g_grad, 'ko', 'markerfacecolor', 'k');
%plot(g2, g_loclin, 'r*', 'markerfacecolor', 'r');
%xlabel('$$\mathrm{\eta^{\top}x}$$','interpreter','latex','fontsize',20);
%ylabel('$$\mathrm{G(\eta^{\top}x)}$$','interpreter','latex','fontsize',20);
%leg = legend('$\mathrm{Gradient~Based}$', '$\mathrm{Gradient~Free}$');
%set(leg,'interpreter','latex','fontsize',16,'location','NorthWest');
%set(gca,'TickLabelInterpreter','Latex','fontsize', 18);
%%title('comparing SSPs');
%box on;
%print -depsc comp_ssp.eps
%
%% compare eigenvalues
%figure;
%hold on;
%semilogy(1:length(lambda_grad),abs(lambda_grad)./lambda_grad(1),'-o','linewidth',2,'MarkerFaceColor','k');
%semilogy(1:length(lambda_loclin),abs(lambda_loclin)./lambda_loclin(1),'-*','linewidth',2,'MarkerFaceColor','r');
%xlabel('$$\mathrm{Index~(i)}$$','interpreter','latex','fontsize',20);
%ylabel('$$\mathrm{Eigenvalue~(\lambda_i)}$$','interpreter','latex','fontsize',20);
%leg = legend('$\mathrm{Gradient~Based}$', '$\mathrm{Gradient~Free}$');
%set(leg,'interpreter','latex','fontsize',16);
%set(gca,'TickLabelInterpreter','Latex','fontsize', 18);
%%title('comparing dominant eigenvalues');
%box on;
%print -depsc comp_eig.eps


% compare dominant eigenvectors
%vi = W(:,1);
%save('eigv.txt','vi','-ASCII');

% compute relative norm of the difference between e.v. components
ev1 = load('eigv_data_files/eigv1.txt');
ev2 = load('eigv_data_files/eigv2.txt');
ev3 = load('eigv_data_files/eigv3.txt');
ev4 = load('eigv_data_files/eigv4.txt');
ev5 = load('eigv_data_files/eigv5.txt');
ev6 = load('eigv_data_files/eigv6.txt');

rel_nor = zeros(6,1);
ev = [ev1 ev2 ev3 ev4 ev5 ev6];
ev_diff = [sqrt((ev1-ev2).^2) sqrt((ev2-ev3).^2) sqrt((ev3-ev4).^2) sqrt((ev4-ev5).^2) sqrt((ev5-ev6).^2)];
max_ev_diff = [max(ev_diff(:,1)) max(ev_diff(:,2)) max(ev_diff(:,3)) max(ev_diff(:,4)) max(ev_diff(:,5))];
global_max = max(max_ev_diff);
rel_nor = [sqrt(sum(ev_diff(:,1).^2)) sqrt(sum(ev_diff(:,2).^2)) sqrt(sum(ev_diff(:,3).^2))...
           sqrt(sum(ev_diff(:,4).^2)) sqrt(sum(ev_diff(:,5).^2))];
for i = 1:5
  rel_nor(i) = rel_nor(i)./global_max;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,1)
hold on;
% line them up: 
%V(:,1) = V(:,1)*sign(V(1,1));
%W(:,1) = W(:,1)*sign(W(1,1));

plot(V(:,1), '-ko', 'linewidth',2, 'markerfacecolor', 'k');
plot(W(:,1), '-r*', 'linewidth',2, 'markerfacecolor', 'r');
xlabel('$$\mathrm{Index~(i)}$$','interpreter','latex','fontsize',12);
ylabel('$$\mathrm{Eigenvector~Components~(v_i)}$$','interpreter','latex','fontsize',12);
leg = legend('$\mathrm{Gradient~Based}$', '$\mathrm{Gradient~Free}$','location','SouthEast');
set(leg,'interpreter','latex','fontsize',8);
set(gca,'TickLabelInterpreter','Latex','fontsize', 8);
%title('comparing dominant eigenvectors');
box on;
axis square;

subplot(1,2,2)
%plot(rel_nor(1:1),':ko','MarkerFaceColor','k');
%xticklabels({'$\mathrm{5-10}$','$\mathrm{10-15}$','$\mathrm{15-20}$','$\mathrm{20-25}$',...
%'$\mathrm{25-30}$'});
xticklabels({'$\mathrm{5-10}$'});
xlabel('$$\mathrm{N_{i} - N_{i+1}}$$','interpreter','latex','fontsize',12);
ylabel('$$\mathrm{\frac{\sqrt{\sum\limits_{j=1}^{N_p}(v_j^i - v_j^{i+1})^2}}{\|(|v^i - v^{i+1}|)\|_\infty}}$$','interpreter','latex','fontsize',12);
set(gca,'xtick',1:1,'fontsize',6,'TickLabelInterpreter','latex');
box on;
axis square;

print -depsc eigv1.eps

%% compare normalized eigenvalues
%l1 = lambda_grad./lambda_grad(1);
%l2 = lambda_loclin./lambda_loclin(1);
%figure;
%bar([abs(l1(:)) abs(l2(:))]);
%ylim([1e-16 1.5]);
%title('comparing eigenvalues');
%set(gca, 'yscale', 'log');
%set(gca, 'fontsize', 20);
%legend('alg 1.1', 'alg 1.2');

% activity scores (1.1) 
%ndim = 19;
%as_grad = zeros(ndim,1);
%for i = 1:ndim
%  for j=1:1
%    as_grad(i) = as_grad(i) + lambda_grad(j).*(V(i,j).^2);
%  end
%end
%as_grad = as_grad / sum(as_grad);
%
%% activity scores (1.1) 
%as_loclin = zeros(ndim,1);
%for i = 1:ndim
%  for j=1:1
%    as_loclin(i) = as_loclin(i) + lambda_loclin(j).*(W(i,j).^2);
%  end
%end
%as_loclin = as_loclin / sum(as_loclin);
%
%figure;
%bar([as_grad(:) as_loclin(:)]);
%xticklabels({'$\mathrm{A_1}$','$\mathrm{A_2}$','$\mathrm{A_3}$','$\mathrm{A_4}$',...
%             '$\mathrm{A_5}$','$\mathrm{A_6}$','$\mathrm{A_7}$','$\mathrm{A_8}$',...
%             '$\mathrm{A_9}$','$\mathrm{A_{10}}$','$\mathrm{A_{11}}$','$\mathrm{A_{12}}$',...
%             '$\mathrm{A_{13}}$','$\mathrm{A_{14}}$','$\mathrm{A_{15}}$','$\mathrm{A_{16}}$',...
%             '$\mathrm{A_{17}}$','$\mathrm{A_{18}}$','$\mathrm{A_{19}}$'});
%leg = legend('$\mathrm{Gradient~Based}$', '$\mathrm{Gradient~Free}$');
%set(leg,'interpreter','latex','fontsize',16);
%set(gca,'xtick',1:19,'fontsize',10,'TickLabelInterpreter','latex');
%ylabel('$$\mathrm{Activity~Scores}$$','interpreter','latex','fontsize',20);
%box on;
%print -depsc comp_as.eps

