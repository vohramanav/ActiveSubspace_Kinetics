close all
clear all

eigN = load('eig_N1000.txt');
eiga2 = load('eig_a2.txt');
eiga3 = load('eig_a3.txt');
eiga4 = load('eig_a4.txt');

eigvN = load('eigv_N1000.txt');
eigva2 = load('eigv_a2.txt');
eigva3 = load('eigv_a3.txt');
eigva4 = load('eigv_a4.txt');

E = zeros(3,1);
E(1) = norm(eigva2(:,1).^2 - eigvN(:,1).^2);
E(2) = norm(eigva3(:,1).^2 - eigvN(:,1).^2);
E(3) = norm(eigva4(:,1).^2 - eigvN(:,1).^2);

figure;
semilogy(eigN(1:19),'--o','MarkerFaceColor','k','linewidth',2);
hold on;
semilogy(eiga2(1:19),'--d','MarkerFaceColor','r','linewidth',2);
semilogy(eiga3(1:19),'--s','MarkerFaceColor','b','linewidth',2);
semilogy(eiga4(1:19),'--*','MarkerFaceColor','g','linewidth',2);
xlabel('$$\mathrm{Index~(i)}$$','interpreter','latex','fontsize',20);
ylabel('$$\mathrm{\|v_{1,n}^2 - v_{1,N}^2\|}$$','interpreter','latex','fontsize',20);
leg = legend('$\mathrm{N~=~1000}$', '$\mathrm{N~=~117~(\alpha:2)}$',...
      '$\mathrm{N~=~176~(\alpha:3)}$', '$\mathrm{N~=~235~(\alpha:4)}$');
set(leg,'interpreter','latex','fontsize',16);
set(gca,'TickLabelInterpreter','Latex','fontsize', 18);
%title('comparing dominant eigenvalues');
box on;
print -depsc eig_comp.eps

figure;
plot(E,'--o','MarkerFaceColor','k');
xlabel('$$\mathrm{\alpha}$$','interpreter','latex','fontsize',20);
ylabel('$$\mathrm{\|v_{1,n}^2 - v_{1,N}^2\|}$$','interpreter','latex','fontsize',20);
%leg = legend('$\mathrm{N~=~117~(\alpha:2)}$','$\mathrm{N~=~176~(\alpha:3)}$',...
%             '$\mathrm{N~=~235~(\alpha:4)}$');
xticklabels({'$\mathrm{2}$','$\mathrm{3}$','$\mathrm{4}$'});
%xlim([1,5]);
%set(leg,'interpreter','latex','fontsize',16);
set(gca,'xtick',1:3,'TickLabelInterpreter','Latex','fontsize', 18);
%title('comparing dominant eigenvalues');
box on;
print -depsc eigv_comp.eps
