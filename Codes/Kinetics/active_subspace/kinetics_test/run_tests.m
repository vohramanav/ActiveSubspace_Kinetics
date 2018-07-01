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

% compare eigenvalues
figure;
hold on;
g1 = V(:,1)'*xi_grad';
g2 = W(:,1)'*xi_loclin';
plot(g1, g_grad, 'ko', 'markerfacecolor', 'k');
plot(g2, g_loclin, 'r*', 'markerfacecolor', 'r');
legend('alg 1.1', 'alg 1.2');
set(gca, 'fontsize', 20);
title('comparing SSPs');
print -dpng comp_ssp.png

% compare SSPs
figure;
hold on;

semilogy(abs(lambda_grad)./lambda_grad(1),'-o','linewidth',2);
semilogy(abs(lambda_loclin)./lambda_loclin(1),'-*','linewidth',2);
legend('alg 1.1', 'alg 1.2');
set(gca, 'fontsize', 20);
title('comparing dominant eigenvalues');
print -dpng comp_eig.png


% compare dominant eigenvectors
figure;
hold on;

% line them up: 
%V(:,1) = V(:,1)*sign(V(1,1));
%W(:,1) = W(:,1)*sign(W(1,1));

plot(V(:,1), '-o', 'linewidth',2);
plot(W(:,1), '-*', 'linewidth',2);

legend('alg 1.1', 'alg 1.2');
set(gca, 'fontsize', 20);
title('comparing dominant eigenvectors');
print -dpng comp_eigv.png

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
ndim = 19;
as_grad = zeros(ndim,1);
for i = 1:ndim
  for j=1:1
    as_grad(i) = as_grad(i) + lambda_grad(j).*(V(i,j).^2);
  end
end
as_grad = as_grad / sum(as_grad);

% activity scores (1.1) 
as_loclin = zeros(ndim,1);
for i = 1:ndim
  for j=1:1
    as_loclin(i) = as_loclin(i) + lambda_loclin(j).*(W(i,j).^2);
  end
end
as_loclin = as_loclin / sum(as_loclin);

figure;
bar([as_grad(:) as_loclin(:)]);
legend('alg 1.1', 'alg 1.2');
set(gca, 'fontsize', 20);
print -dpng comp_as.png

