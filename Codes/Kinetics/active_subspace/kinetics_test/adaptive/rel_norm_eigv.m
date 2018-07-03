close all
clear all

ev1 = load('eigv1.txt');
ev2 = load('eigv2.txt');
ev3 = load('eigv3.txt');
ev4 = load('eigv4.txt');
ev5 = load('eigv5.txt');
ev6 = load('eigv6.txt');

rel_nor = zeros(6,1);
ev = [ev1 ev2 ev3 ev4 ev5 ev6];
ev_diff = [ev1-ev2 ev2-ev3 ev3-ev4 ev4-ev5 ev5-ev6];
rel_nor = [sqrt(sum(ev_diff(:,1).^2)) sqrt(sum(ev_diff(:,2).^2)) sqrt(sum(ev_diff(:,3).^2))...
           sqrt(sum(ev_diff(:,4).^2)) sqrt(sum(ev_diff(:,5).^2))];
rel_nor = rel_nor./max(rel_nor);

figure;
%plot(rel_nor(1:1),':ko','MarkerFaceColor','k');
%xticklabels({'$\mathrm{5-10}$','$\mathrm{10-15}$','$\mathrm{15-20}$','$\mathrm{20-25}$',...
%'$\mathrm{25-30}$'});
xticklabels({'$\mathrm{5-10}$'});
xlabel('$$\mathrm{N_{i} - N_{i+1}}$$','interpreter','latex','fontsize',20);
ylabel('$$\mathrm{\frac{\sqrt{\sum\limits_{j=1}^{N_p}(v_j^i - v_j^{i+1})^2}}{\|(|v_j^i - v_j^{i+1}|)\|_\infty}}$$','interpreter','latex','fontsize',20);
set(gca,'xtick',1:1,'fontsize',12,'TickLabelInterpreter','latex');
box on;
print -depsc rel_nor1.eps
