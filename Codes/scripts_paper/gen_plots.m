close all
clear all

%[rel_norm_per,rel_norm_reg] = eig_conv();
%[Y_surr_per,Y_surr_reg] = SSP_plot();
%[rel_l2_per,rel_l2_reg] = pdf_plot(Y_surr_per,Y_surr_reg);
[as_per,as_reg] = sensitivity();

function [rel_norm_per,rel_norm_reg] = eig_conv()
% compute relative norm of the difference between e.v. components for regression
ev1 = load('data_files_p2/eigv_gradfreeN185.txt');
ev2 = load('data_files_p2/eigv_gradfreeN370.txt');
ev3 = load('data_files_p2/eigv_gradfreeN555.txt');
ev4 = load('data_files_p2/eigv_gradfreeN740.txt');
ev5 = load('data_files_p2/eigv_gradfreeN925.txt');
ev6 = load('data_files_p2/eigv_gradfreeN1110.txt');
ev7 = load('data_files_p2/eigv_gradfreeN1295.txt');
ev8 = load('data_files_p2/eigv_gradfreeN1480.txt');
ev9 = load('data_files_p2/eigv_gradfreeN1665.txt');
ev10 = load('data_files_p2/eigv_gradfreeN1850.txt');
W = ev10;

ev = [ev1(:,1) ev2(:,1) ev3(:,1) ev4(:,1) ev5(:,1) ev6(:,1) ev7(:,1) ev8(:,1) ev9(:,1) ev10(:,1)];
rel_norm_reg = zeros(size(ev,2)-1,1);
for i = 1:length(rel_norm_reg)
  rel_norm_reg(i,1) = norm(ev(:,i).^2-ev(:,i+1).^2)./norm(ev(:,i).^2);
end

% compute relative norm of the difference between e.v. components for perturbation
ev1 = load('../gradient_based/data_files_p2/eigv_gradN05.txt');
ev2 = load('../gradient_based/data_files_p2/eigv_gradN10.txt');
ev3 = load('../gradient_based/data_files_p2/eigv_gradN15.txt');
ev4 = load('../gradient_based/data_files_p2/eigv_gradN20.txt');
ev5 = load('../gradient_based/data_files_p2/eigv_gradN25.txt');
ev6 = load('../gradient_based/data_files_p2/eigv_gradN30.txt');
ev7 = load('../gradient_based/data_files_p2/eigv_gradN35.txt');
ev8 = load('../gradient_based/data_files_p2/eigv_gradN40.txt');
ev9 = load('../gradient_based/data_files_p2/eigv_gradN45.txt');
ev10 = load('../gradient_based/data_files_p2/eigv_gradN50.txt');
V = ev5;

% convergence of perturbation approach
ev = [ev1(:,1) ev2(:,1) ev3(:,1) ev4(:,1) ev5(:,1) ev6(:,1) ev7(:,1) ev8(:,1) ev9(:,1) ev10(:,1)];
rel_norm_per = zeros(size(ev,2)-1,1);
for i = 1:length(rel_norm_reg)
  rel_norm_per(i,1) = norm(ev(:,i).^2-ev(:,i+1).^2)./norm(ev(:,i).^2);
end

subplot(1,2,1)
hold on;
plot(V(:,1).^2, '--ko');
plot(W(:,1).^2, '--r*');
xlabel('$$\mathrm{Index~(i)}$$','interpreter','latex','fontsize',20);
ylabel('$$\mathrm{w_i^2}$$','interpreter','latex','fontsize',20);
ylim([-0.1,1.0]);
leg = legend({'$\mathrm{Perturbation}$','$\mathrm{Regression}$'});
set(leg,'interpreter','latex','fontsize',12);
set(gca,'TickLabelInterpreter','Latex','fontsize', 16);
box on;
axis square;

hold off;
subplot(1,2,2)
hold on;
plot(1:4,rel_norm_per(1:4),':ko');
plot(1:9,rel_norm_reg(1:9),':r*');
xlabel('$$\mathrm{Iterations}$$','interpreter','latex','fontsize',20);
ylabel('$$\mathrm{\max\limits_j(\delta \hat{\mathbf{W}}_{1,j}^{(i)})}$$','interpreter','latex','fontsize',20);
xlim([0,10]);
xticklabels({'$\mathrm{1}$','$\mathrm{2}$','$\mathrm{3}$','$\mathrm{4}$',...
     '$\mathrm{5}$','$\mathrm{6}$','$\mathrm{7}$','$\mathrm{8}$',...
     '$\mathrm{9}$'});
leg = legend({'$\mathrm{Perturbation}$','$\mathrm{Regression}$'});
set(leg,'interpreter','latex','fontsize',12);
set(gca,'xtick',1:9,'fontsize',14,'TickLabelInterpreter','latex');
box on;
axis square;

print -depsc eig_conv36Dp2.eps
end

function [Y_surr_per,Y_surr_reg] = SSP_plot()
W_reg = load('data_files_p2/eigv_gradfreeN1850.txt');
lam_reg = load('data_files_p2/lambda_gradfreeN1850.txt');
xi_reg = load('data_files_p2/xi_gradfreeN1850.txt');
G_reg = load('data_files_p2/f_gradfreeN1850.txt');

W_per = load('../gradient_based/data_files_p2/eigv_gradN25.txt');
lam_per = load('../gradient_based/data_files_p2/lambda_gradN25.txt');
xi_per = load('../gradient_based/data_files_p2/xi_gradN25.txt');
G_per = load('../gradient_based/data_files_p2/f_gradN25.txt');

g1 = -W_per(:,1)'*xi_per';
g2 = W_reg(:,1)'*xi_reg';

% fit a surrogate
ssp_per = zeros(length(g1),2);
ssp_per(:,1) = g1; ssp_per(:,2) = G_per;
ssp_per_sort = sortrows(ssp_per);
Y_surr_per = get_polyfit_surr(ssp_per_sort(:,1),ssp_per_sort(:,2),-W_per(:,1),2);

ssp_reg = zeros(length(g2),2);
ssp_reg(:,1) = g2; ssp_reg(:,2) = G_reg;
ssp_reg_sort = sortrows(ssp_reg);
Y_surr_reg = get_polyfit_surr(ssp_reg_sort(:,1),ssp_reg_sort(:,2),W_reg(:,1),3);

np1 = length(g1); np2 = length(g2);
Y_surr_comp_per = zeros(np1,1);
Y_surr_comp_reg = zeros(np2,1);
for i = 1:np1
  Y_surr_data_per(i,1) = Y_surr_per(xi_per(i,:)');
end
ssp_per_data = zeros(length(g1),2);
ssp_per_data(:,1) = g1; ssp_per_data(:,2) = Y_surr_data_per;
ssp_per_data_sort = sortrows(ssp_per_data);

for i = 1:np2
  Y_surr_data_reg(i,1) = Y_surr_reg(xi_reg(i,:)');
end
ssp_reg_data = zeros(length(g2),2);
ssp_reg_data(:,1) = g2; ssp_reg_data(:,2) = Y_surr_data_reg;
ssp_reg_data_sort = sortrows(ssp_reg_data);

%hold off;
%subplot(1,3,1)
%hold on;
%plot(g1,G_per,'bo','MarkerFaceColor','k','MarkerSize',5);
%plot(ssp_per_data_sort(:,1),ssp_per_data_sort(:,2),'--r','LineWidth',2);
%xlabel('$$\mathrm{\hat{\mathbf{W}}_1^{\top}\mathbf{\xi}}$$','interpreter','latex','fontsize',20);
%ylabel('$$\mathrm{G(\hat{\mathbf{W}}_1^{\top}\mathbf{\xi})}$$','interpreter','latex','fontsize',20);
%ylim([0,0.6]);
%leg = legend({'$\mathrm{G}$','$\mathrm{\hat{G}}$'});
%set(leg,'interpreter','latex','fontsize',12,'location','NorthWest');
%set(gca,'TickLabelInterpreter','Latex','fontsize', 16);
%box on;
%axis square;
%
%hold off;
%subplot(1,3,2)
%hold on;
%plot(g2,G_reg,'bo','MarkerFaceColor','k','MarkerSize',5);
%plot(ssp_reg_data_sort(:,1),ssp_reg_data_sort(:,2),'--r','LineWidth',2);
%xlabel('$$\mathrm{\hat{\mathbf{W}}_1^{\top}\mathbf{\xi}}$$','interpreter','latex','fontsize',20);
%ylabel('$$\mathrm{G(\hat{\mathbf{W}}_1^{\top}\mathbf{\xi})}$$','interpreter','latex','fontsize',20);
%leg = legend({'$\mathrm{G}$','$\mathrm{\hat{G}}$'});
%set(leg,'interpreter','latex','fontsize',12,'location','NorthWest');
%set(gca,'TickLabelInterpreter','Latex','fontsize', 16);
%box on;
%axis square;
%
%hold off;
%subplot(1,3,3)
%hold on;
%plot(ssp_per_data_sort(:,1),ssp_per_data_sort(:,2),'--r','LineWidth',3);
%plot(ssp_reg_data_sort(:,1),ssp_reg_data_sort(:,2),'--b','LineWidth',1.5);
%xlabel('$$\mathrm{\hat{\mathbf{W}}_1^{\top}\mathbf{\xi}}$$','interpreter','latex','fontsize',20);
%ylabel('$$\mathrm{G(\hat{\mathbf{W}}_1^{\top}\mathbf{\xi})}$$','interpreter','latex','fontsize',20);
%leg = legend({'$\mathrm{\hat{G}:~Perturbation}$','$\mathrm{\hat{G}:~Regression}$'});
%set(leg,'interpreter','latex','fontsize',12,'location','NorthWest');
%set(gca,'TickLabelInterpreter','Latex','fontsize', 16);
%box on;
%axis square;
%print -depsc SSP_plot36Dp2.eps

end

function [rel_l2_per,rel_l2_reg] = pdf_plot(Y_surr_per,Y_surr_reg)
xi_ver = load('data_files_p2/xi_N10000ver.txt');
np = 10000;
data_per = zeros(np,1); data_reg = zeros(np,1);

for i = 1:np
  data_per(i,1) = Y_surr_per(xi_ver(i,:)');
  data_reg(i,1) = Y_surr_reg(xi_ver(i,:)');
end

G = load('data_files_p2/record_id_N10000ver.txt');
mu = zeros(3,1); stdev = zeros(3,1);
mu(1,1) = mean(G); stdev(1,1) = std(G);
mu(2,1) = mean(data_per); stdev(2,1) = std(data_per);
mu(3,1) = mean(data_reg); stdev(3,1) = std(data_reg);
save('mu.txt','mu','-ASCII');
save('stdev.txt','stdev','-ASCII');
pdf_data = zeros(np,3);
pdf_data(:,1) = G; pdf_data(:,2) = data_per; pdf_data(:,3) = data_reg;
save('pdf_data','pdf_data','-ASCII'); 

[density,xmesh] = ksdensity(G);
[density_per,xmesh_per] = ksdensity(data_per);
[density_reg,xmesh_reg] = ksdensity(data_reg);

rel_l2_per = norm(G-data_per)./norm(G);
rel_l2_reg = norm(G-data_reg)./norm(G);

hold off;
figure;
hold on;
%nbins = 30;
plot(xmesh,density,'-b','LineWidth',2);
plot(xmesh_per,density_per,'--k','LineWidth',2);
plot(xmesh_reg,density_reg,'-.r','LineWidth',2);
xlabel('$$\mathrm{Ignition~Delay (s)}$$','interpreter','latex','fontsize',20);
ylabel('$$\mathrm{PDF}$$','interpreter','latex','fontsize',20);
%xlim([-0.1,1.0]);
leg=legend('$\mathrm{Model}$','$\mathrm{Perturbation}$','$\mathrm{Regression}$');
set(leg,'interpreter','latex','fontsize',14,'location','NorthEast');
set(gca,'fontsize',18,'TickLabelInterpreter','latex');
box on;
grid on;
print -depsc pdf_plot36Dp2_1e5.eps
end

function [as_grad,as_loclin] = sensitivity()

% load data
W_reg = load('data_files_p2/eigv_gradfreeN1850.txt');
lam_reg = load('data_files_p2/lambda_gradfreeN1850.txt');
W_per = load('../gradient_based/data_files_p2/eigv_gradN25.txt');
lam_per = load('../gradient_based/data_files_p2/lambda_gradN25.txt');
xi_ref = load('data_files_p2/xi_N10000ver.txt');
f = load('data_files_p2/record_id_N10000ver.txt');

%% Activity Scores
m = 36;
as_loclin = zeros(m,1);
as_grad = zeros(m,1);

for i = 1:m
  for j=1:1
    as_loclin(i) = as_loclin(i) + lam_reg(j).*(W_reg(i,j).^2);
    as_grad(i) = as_grad(i) + lam_per(j).*(W_per(i,j).^2);
  end
end
as_loclin = as_loclin / sum(as_loclin);
as_grad = as_grad / sum(as_grad);

% Computing total Sobol index using eigenvector from the perturbation approach
eta1g = W_per(:,1);
avg = eta1g'*xi_ref';
Y_surrg = get_polyfit_surr(avg,f,eta1g,1);

N = 1e5;
An = -1 + 2 * rand(N , m);
Bn = -1 + 2 * rand(N , m);

% constructing the modified matrices A(B)i and B(A)i
An_Bn_1 = An; An_Bn_1(:,1) = Bn(:,1);
An_Bn_2 = An; An_Bn_2(:,2) = Bn(:,2);
An_Bn_3 = An; An_Bn_3(:,3) = Bn(:,3);
An_Bn_4 = An; An_Bn_4(:,4) = Bn(:,4);
An_Bn_5 = An; An_Bn_5(:,5) = Bn(:,5);
An_Bn_6 = An; An_Bn_6(:,6) = Bn(:,6);
An_Bn_7 = An; An_Bn_7(:,7) = Bn(:,7);
An_Bn_8 = An; An_Bn_8(:,8) = Bn(:,8);
An_Bn_9 = An; An_Bn_9(:,9) = Bn(:,9);
An_Bn_10 = An; An_Bn_10(:,10) = Bn(:,10);
An_Bn_11 = An; An_Bn_11(:,11) = Bn(:,11);
An_Bn_12 = An; An_Bn_12(:,12) = Bn(:,12);
An_Bn_13 = An; An_Bn_13(:,13) = Bn(:,13);
An_Bn_14 = An; An_Bn_14(:,14) = Bn(:,14);
An_Bn_15 = An; An_Bn_15(:,15) = Bn(:,15);
An_Bn_16 = An; An_Bn_16(:,16) = Bn(:,16);
An_Bn_17 = An; An_Bn_17(:,17) = Bn(:,17);
An_Bn_18 = An; An_Bn_18(:,18) = Bn(:,18);
An_Bn_19 = An; An_Bn_19(:,19) = Bn(:,19);
An_Bn_20 = An; An_Bn_20(:,20) = Bn(:,20);
An_Bn_21 = An; An_Bn_21(:,21) = Bn(:,21);
An_Bn_22 = An; An_Bn_22(:,22) = Bn(:,22);
An_Bn_23 = An; An_Bn_23(:,23) = Bn(:,23);
An_Bn_24 = An; An_Bn_24(:,24) = Bn(:,24);
An_Bn_25 = An; An_Bn_25(:,25) = Bn(:,25);
An_Bn_26 = An; An_Bn_26(:,26) = Bn(:,26);
An_Bn_27 = An; An_Bn_27(:,27) = Bn(:,27);
An_Bn_28 = An; An_Bn_28(:,28) = Bn(:,28);
An_Bn_29 = An; An_Bn_29(:,29) = Bn(:,29);
An_Bn_30 = An; An_Bn_30(:,30) = Bn(:,30);
An_Bn_31 = An; An_Bn_31(:,31) = Bn(:,31);
An_Bn_32 = An; An_Bn_32(:,32) = Bn(:,32);
An_Bn_33 = An; An_Bn_33(:,33) = Bn(:,33);
An_Bn_34 = An; An_Bn_34(:,34) = Bn(:,34);
An_Bn_35 = An; An_Bn_35(:,35) = Bn(:,35);
An_Bn_36 = An; An_Bn_33(:,36) = Bn(:,36);

Bn_An_1 = Bn; Bn_An_1(:,1) = An(:,1);
Bn_An_2 = Bn; Bn_An_2(:,2) = An(:,2);
Bn_An_3 = Bn; Bn_An_3(:,3) = An(:,3);
Bn_An_4 = Bn; Bn_An_4(:,4) = An(:,4);
Bn_An_5 = Bn; Bn_An_5(:,5) = An(:,5);
Bn_An_6 = Bn; Bn_An_6(:,6) = An(:,6);
Bn_An_7 = Bn; Bn_An_7(:,7) = An(:,7);
Bn_An_8 = Bn; Bn_An_8(:,8) = An(:,8);                      
Bn_An_9 = Bn; Bn_An_9(:,9) = An(:,9);
Bn_An_10 = Bn; Bn_An_10(:,10) = An(:,10);
Bn_An_11 = Bn; Bn_An_11(:,11) = An(:,11);
Bn_An_12 = Bn; Bn_An_12(:,12) = An(:,12);
Bn_An_13 = Bn; Bn_An_13(:,13) = An(:,13);
Bn_An_14 = Bn; Bn_An_14(:,14) = An(:,14);
Bn_An_15 = Bn; Bn_An_15(:,15) = An(:,15);
Bn_An_16 = Bn; Bn_An_16(:,16) = An(:,16);
Bn_An_17 = Bn; Bn_An_17(:,17) = An(:,17);
Bn_An_18 = Bn; Bn_An_18(:,18) = An(:,18);
Bn_An_19 = Bn; Bn_An_19(:,19) = An(:,19);
Bn_An_20 = Bn; Bn_An_20(:,20) = An(:,20);
Bn_An_21 = Bn; Bn_An_21(:,21) = An(:,21);
Bn_An_22 = Bn; Bn_An_22(:,22) = An(:,22);
Bn_An_23 = Bn; Bn_An_23(:,23) = An(:,23);
Bn_An_24 = Bn; Bn_An_24(:,24) = An(:,24);
Bn_An_25 = Bn; Bn_An_25(:,25) = An(:,25);
Bn_An_26 = Bn; Bn_An_26(:,26) = An(:,26);
Bn_An_27 = Bn; Bn_An_27(:,27) = An(:,27);
Bn_An_28 = Bn; Bn_An_28(:,28) = An(:,28);
Bn_An_29 = Bn; Bn_An_29(:,29) = An(:,29);
Bn_An_30 = Bn; Bn_An_30(:,30) = An(:,30);
Bn_An_31 = Bn; Bn_An_31(:,31) = An(:,31);
Bn_An_32 = Bn; Bn_An_32(:,32) = An(:,32);
Bn_An_33 = Bn; Bn_An_33(:,33) = An(:,33);
Bn_An_34 = Bn; Bn_An_34(:,34) = An(:,34);
Bn_An_35 = Bn; Bn_An_35(:,35) = An(:,35);
Bn_An_36 = Bn; Bn_An_36(:,36) = An(:,36);

%surrogate evaluations
fAn = zeros(N,1);
fAn_Bn = zeros(N,m); fBn_An = zeros(N,m);
for i = 1:N
  fAn(i,1) = Y_surrg(An(i,:)');
  fAn_Bn(i,1) = Y_surrg(An_Bn_1(i,:)'); fBn_An(i,1) = Y_surrg(Bn_An_1(i,:)');
  fAn_Bn(i,2) = Y_surrg(An_Bn_2(i,:)'); fBn_An(i,2) = Y_surrg(Bn_An_2(i,:)');
  fAn_Bn(i,3) = Y_surrg(An_Bn_3(i,:)'); fBn_An(i,3) = Y_surrg(Bn_An_3(i,:)');
  fAn_Bn(i,4) = Y_surrg(An_Bn_4(i,:)'); fBn_An(i,4) = Y_surrg(Bn_An_4(i,:)');
  fAn_Bn(i,5) = Y_surrg(An_Bn_5(i,:)'); fBn_An(i,5) = Y_surrg(Bn_An_5(i,:)');
  fAn_Bn(i,6) = Y_surrg(An_Bn_6(i,:)'); fBn_An(i,6) = Y_surrg(Bn_An_6(i,:)');
  fAn_Bn(i,7) = Y_surrg(An_Bn_7(i,:)'); fBn_An(i,7) = Y_surrg(Bn_An_7(i,:)');
  fAn_Bn(i,8) = Y_surrg(An_Bn_8(i,:)'); fBn_An(i,8) = Y_surrg(Bn_An_8(i,:)');
  fAn_Bn(i,9) = Y_surrg(An_Bn_9(i,:)'); fBn_An(i,9) = Y_surrg(Bn_An_9(i,:)');
  fAn_Bn(i,10) = Y_surrg(An_Bn_10(i,:)'); fBn_An(i,10) = Y_surrg(Bn_An_10(i,:)');
  fAn_Bn(i,11) = Y_surrg(An_Bn_11(i,:)'); fBn_An(i,11) = Y_surrg(Bn_An_11(i,:)');
  fAn_Bn(i,12) = Y_surrg(An_Bn_12(i,:)'); fBn_An(i,12) = Y_surrg(Bn_An_12(i,:)');
  fAn_Bn(i,13) = Y_surrg(An_Bn_13(i,:)'); fBn_An(i,13) = Y_surrg(Bn_An_13(i,:)');
  fAn_Bn(i,14) = Y_surrg(An_Bn_14(i,:)'); fBn_An(i,14) = Y_surrg(Bn_An_14(i,:)');
  fAn_Bn(i,15) = Y_surrg(An_Bn_15(i,:)'); fBn_An(i,15) = Y_surrg(Bn_An_15(i,:)');
  fAn_Bn(i,16) = Y_surrg(An_Bn_16(i,:)'); fBn_An(i,16) = Y_surrg(Bn_An_16(i,:)');
  fAn_Bn(i,17) = Y_surrg(An_Bn_17(i,:)'); fBn_An(i,17) = Y_surrg(Bn_An_17(i,:)');
  fAn_Bn(i,18) = Y_surrg(An_Bn_18(i,:)'); fBn_An(i,18) = Y_surrg(Bn_An_18(i,:)');
  fAn_Bn(i,19) = Y_surrg(An_Bn_19(i,:)'); fBn_An(i,19) = Y_surrg(Bn_An_19(i,:)');
  fAn_Bn(i,20) = Y_surrg(An_Bn_20(i,:)'); fBn_An(i,20) = Y_surrg(Bn_An_20(i,:)');
  fAn_Bn(i,21) = Y_surrg(An_Bn_21(i,:)'); fBn_An(i,21) = Y_surrg(Bn_An_21(i,:)');
  fAn_Bn(i,22) = Y_surrg(An_Bn_22(i,:)'); fBn_An(i,22) = Y_surrg(Bn_An_22(i,:)');
  fAn_Bn(i,23) = Y_surrg(An_Bn_23(i,:)'); fBn_An(i,23) = Y_surrg(Bn_An_23(i,:)');
  fAn_Bn(i,24) = Y_surrg(An_Bn_24(i,:)'); fBn_An(i,24) = Y_surrg(Bn_An_24(i,:)');
  fAn_Bn(i,25) = Y_surrg(An_Bn_25(i,:)'); fBn_An(i,25) = Y_surrg(Bn_An_25(i,:)');
  fAn_Bn(i,26) = Y_surrg(An_Bn_26(i,:)'); fBn_An(i,26) = Y_surrg(Bn_An_26(i,:)');
  fAn_Bn(i,27) = Y_surrg(An_Bn_27(i,:)'); fBn_An(i,27) = Y_surrg(Bn_An_27(i,:)');
  fAn_Bn(i,28) = Y_surrg(An_Bn_28(i,:)'); fBn_An(i,28) = Y_surrg(Bn_An_28(i,:)');
  fAn_Bn(i,29) = Y_surrg(An_Bn_29(i,:)'); fBn_An(i,29) = Y_surrg(Bn_An_29(i,:)');
  fAn_Bn(i,30) = Y_surrg(An_Bn_30(i,:)'); fBn_An(i,30) = Y_surrg(Bn_An_30(i,:)');
  fAn_Bn(i,31) = Y_surrg(An_Bn_31(i,:)'); fBn_An(i,31) = Y_surrg(Bn_An_31(i,:)');
  fAn_Bn(i,32) = Y_surrg(An_Bn_32(i,:)'); fBn_An(i,32) = Y_surrg(Bn_An_32(i,:)');
  fAn_Bn(i,33) = Y_surrg(An_Bn_33(i,:)'); fBn_An(i,33) = Y_surrg(Bn_An_33(i,:)');
  fAn_Bn(i,34) = Y_surrg(An_Bn_34(i,:)'); fBn_An(i,34) = Y_surrg(Bn_An_34(i,:)');
  fAn_Bn(i,35) = Y_surrg(An_Bn_35(i,:)'); fBn_An(i,35) = Y_surrg(Bn_An_35(i,:)');
  fAn_Bn(i,36) = Y_surrg(An_Bn_36(i,:)'); fBn_An(i,36) = Y_surrg(Bn_An_36(i,:)');
end
dim=m;
f_total = zeros((dim+1).*N,1);
f_total(1:N,1) = fAn;
f_total(N+1:end,1) = reshape(fBn_An,N*dim,1);
f0_m = mean(f_total); Dr_m = var(f_total);
f_total(N+1:end,1) = reshape(fAn_Bn,N*dim,1);
f0_t = mean(f_total); Dr_t = var(f_total);

Nr_m = zeros(dim,1); Nr_t = zeros(dim,1); gsa_m = zeros(dim,1); gsa_t = zeros(dim,1);

for i = 1:dim
   Nr_m(i) = (1.0./N).*dot(fAn,fBn_An(:,i)) - f0_m.^2.0;
%   Nr_t(i) = (1.0./N).*dot(fAn,fAn_Bn(:,i)) - f0_t.^2.0;
   diff = fAn - fAn_Bn(:,i);
   Nr_t(i) = (1.0/N).*dot(fAn,diff);
end

gsa_m = Nr_m./Dr_m;
%gsa_t = 1.0 - Nr_t./Dr_t;
gsa_t = Nr_t./Dr_t;

% Plot1
%hold off;
%subplot(1,2,1)
%bar([as_loclin(1:19) as_grad(1:19) gsa_t(1:19)]);
%%bar([as_loclin(1:19) as_grad(1:19)]);
%xticklabels({'$\mathrm{A_1}$','','','',...
%             '','','','',...
%             '$\mathrm{A_9}$','','','',...
%             '','','$\mathrm{A_{15}}$','',...
%             '$\mathrm{A_{17}}$','',''});
%ylabel('$$\mathrm{Activity~Scores~(\tilde{\nu}_{i,r})}$$','interpreter','latex','fontsize',14);
%ylim([0 0.35]);
%leg=legend('$\mathrm{Regression}$','$\mathrm{Perturbation}$','$\mathrm{Total~Sobol~Index}$');
%set(leg,'interpreter','latex','fontsize',9,'location','NorthWest');
%set(gca,'xtick',1:19,'fontsize',12,'TickLabelInterpreter','latex');
%axis square;
%box on;
%
%subplot(1,2,2)
%bar([as_loclin(20:m) as_grad(20:m) gsa_t(20:m)]);
%%bar([as_loclin(20:m) as_grad(20:m)]);
%xticklabels({'','','','',...
%             '','','','',...
%             '','$\mathrm{E_{a,15}}$','','',...
%             '','','$\mathrm{P_0}$','$\mathrm{T_0}$','$\mathrm{\Phi_0}$'});
%%ylim([0 0.35]);
%leg=legend('$\mathrm{Regression}$','$\mathrm{Perturbation}$','$\mathrm{Total~Sobol~Index}$');
%set(leg,'interpreter','latex','fontsize',9,'location','NorthWest');
%ylim([0 0.35]);
%set(gca,'xtick',1:17,'fontsize',12,'TickLabelInterpreter','latex');
%axis square;
%box on;
%print -depsc as_36D.eps

% Plot2
hold off;
bar([as_loclin(1:m) as_grad(1:m) gsa_t(1:m)],1.2);
%bar([as_loclin(1:19) as_grad(1:19)]);
xticklabels({'$\mathrm{A_1}$','','','',...
             '','','','',...
             '$\mathrm{A_9}$','','','',...
             '','','$\mathrm{A_{15}}$','',...
             '$\mathrm{A_{17}}$','','','','','','',...
             '','','','',...
             '','$\mathrm{E_{a,15}}$','','',...
             '','','$\mathrm{P_0}$','$\mathrm{T_0}$','$\mathrm{\Phi_0}$'});
ylabel('$$\mathrm{Activity~Scores~(\tilde{\nu}_{i,r})}$$','interpreter','latex','fontsize',21);
ylim([0 0.35]);
leg=legend('$\mathrm{Regression}$','$\mathrm{Perturbation}$','$\mathrm{Total-Effect~Sobol~Index}$');
set(leg,'interpreter','latex','fontsize',12,'location','NorthEast');
set(gca,'xtick',1:36,'fontsize',15,'TickLabelInterpreter','latex');
box on;

print -depsc as_36Dp2_rev2.eps
savefig('as_36Dp2_rev2.fig');
end

function surr = get_polyfit_surr(y,G,eta,deg)
p = polyfit(y(:),G,deg);
surr = @(x)(polyval(p,eta'*x));
end
