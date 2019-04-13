rng(123);
m1 = 19;
m2 = 14;
m3 = 3; % uncertain P,T,Phi

L1 = zeros(1,m1); U1 = zeros(1,m1);
L2 = zeros(1,m2); U2 = zeros(1,m2);
Lp = zeros(1,3); Up = zeros(1,3); % for P,T,Phi

% Nominal values for pre-factors a1 to a19 for the 19 reactions

nom1 = [1.915e14,5.080e04,2.160e08,1.230e04,4.577e19,6.165e15,4.714e18,2.240e22,6.170e19,...
       6.630e13,1.690e14,1.810e13,1.450e16,3.020e12,1.202e17,1.000e13,4.820e13,9.550e06,...
       7.000e12];
nom2 = [1.644e4,6.290e3,3.430e3,-1.880e3,1.044e5,...
        2.130e3,8.740e2,-4.000e2,1.390e3,4.550e4,3.590e3,7.950e3,3.970e3,1.430e3];

nomp = [1.0,900,2.0];

nom1 = log(nom1); % perturb log(A) instead of A 

L1(1,:) = 0.98.*nom1(1,:); U1(1,:) = 1.02.*nom1(1,:);
L2(1,:) = 0.98.*nom2(1,:); U2(1,:) = 1.02.*nom2(1,:);
Lp(1,:) = 0.98.*nomp(1,:); Up(1,:) = 1.02.*nomp(1,:);

%Algorithm 1.2 from book
m = m1+m2+m3;
k = m + 1;
alpha = 5;

%N = alpha * m + 1;
N = 100000;
M = floor( alpha *k * log(m) );
%M = N-1;
%M=N;

%
% get the samples
%
xi = -1 + 2 * lhsdesign(N,m); 
y = -1 + 2 * lhsdesign(M,m); 

% FUNCTION CALLS
gen_samples(m,N,xi,U1,L1,U2,L2,Up,Lp);
%[lambda_loclin,W,f,as_loclin,rn1,rn2,Y_surrg] = compute_subspace(m,N,M,xi,y);
%[lambda_loclin,W] = compute_subspace(m,N,M,xi,y);
%command = 'mail -s "run complete" vohra.manav@gmail.com < /dev/null'
%system(command);
%[gsa_m,gsa_t] = sobol(m,xi);

function gen_data = gen_samples(m,N,xi,U1,L1,U2,L2,Up,Lp)
pts_xi = xi;
save('xi_N1e5ver.txt','xi','-ASCII');
pts_x = zeros(size(pts_xi,1),size(pts_xi,2)+6); % 5 E's with 0 nom value + 1 for MO2

% Project points to the physical space

for i = 1:size(pts_x,1)
  for j = 1:19
    pts_x(i,j) = L1(1,j) + 0.5*(U1(1,j)-L1(1,j)).*(pts_xi(i,j)+1);
  end
  for j = 20:24
    pts_x(i,j) = L2(1,j-19) + 0.5*(U2(1,j-19)-L2(1,j-19)).*(pts_xi(i,j)+1);
  end
  for j = 25:28
    pts_x(i,j) = 0.000e0;
  end
  for j = 29:31
    pts_x(i,j) = L2(1,j-23) + 0.5*(U2(1,j-23)-L2(1,j-23)).*(pts_xi(i,j-4)+1);
  end
  pts_x(i,32) = 0.000e0;
  for j = 33:38
    pts_x(i,j) = L2(1,j-24) + 0.5*(U2(1,j-24)-L2(1,j-24)).*(pts_xi(i,j-5)+1);
  end    
  for j = 39:41
    pts_x(i,j) = Lp(1,j-38) + 0.5*(Up(1,j-38)-Lp(1,j-38)).*(pts_xi(i,j-5)+1);
  end
  phi = pts_x(i,41);
  nO2 = 0.9./(2.*phi+1);
  nH2 = 0.9-nO2;
  pts_x(i,41) = nH2;
  pts_x(i,42) = nO2;
end

pts_x(:,1:19) = exp(pts_x(:,1:19));

% Save physical points to a file
save('pts_N1e5ver.txt','pts_x','-ASCII');

end

%function [lambda_loclin,W,f,as_loclin,rn1,rn2,Y_surrg] = compute_subspace(m,N,M,xi,y)
function [lambda_loclin,W] = compute_subspace(m,N,M,xi,y)
id = load('data_files_p2/record_id_gradfreeN1850.txt');
%gsa = load('gsa_t.txt');
%W_grad = load('eigv_gradN40.txt');
%W_grad(:,1) = -W_grad(:,1);
%lam_grad = load('lambda_gradN40.txt');
%f_grad = load('f_gradN40.txt');
f = zeros(N,1);
f(:,1) = id(1:N);
xi_ref = xi;
y_ref = y;

%
% find the nearest p points in N for each point in M
%
%p = N - 1;  %integer between m+1 and N
p = N-1;  %integer between m+1 and N

d_matrix = zeros(N,1);

for i=1:M
    
    for j=1:N
        
        d_matrix(j) = 0;
        
        for k=1:m
            
            d_matrix(j) = d_matrix(j) + norm(y_ref(i,k) - xi_ref(j,k));
               
        end
        
    end
    
    [z,index(i,:)] = sort(d_matrix);
    
    for j=1:p
        
        ip = (i-1)*p + j;
        
        points(ip,:) = xi(index(i,j),:);
        
    end
    
end

%
% formulate the least square problem
%
for np = 1 : M
    
    A = [1 points((np-1)*p+1,:)];
    
    for i = (np-1)*p+2 : np*p
        
        A = [A; 1 points(i,:)];
        
    end
    
    B = f(index(np,1));
    
    for i=2:p
        
        B = [B; f(index(np,i))];
        
    end
    
    z = A \ B;
    
    if np == 1
        
       b_matrix = z(2:m+1);
       
    else
        
       b_matrix = [b_matrix z(2:m+1)];
       
    end
    
end

%construct the covariance matrix

C = 0;

for i=1 : M
    
    z = b_matrix(:,i);
    
    C = C + z * z';
    
end

C=C/M;

[W D] = eig(C);

[lambda_loclin, idx] = sort(diag(D), 'descend');

W = W(:,idx);
eta1 = W(:,1);
eta2 = W(:,2);
save('eigv_gradfreeN1850.txt','W','-ASCII');
save('lambda_gradfreeN1850.txt','lambda_loclin','-ASCII');
save('xiref_gradfreeN1850.txt','xi_ref','-ASCII');
save('f_gradfreeN1850.txt','f','-ASCII');

% Verification
%xi_test = load('pts_data_files/pts_xi_test_1e4.txt');
%xir_test = zeros(size(xi_test,1),m);
%xir_test(:,1:24) = xi_test(:,1:24);
%xir_test(:,25:27) = xi_test(:,29:31);
%xir_test(:,28:33) = xi_test(:,33:38);
%
%% 1D Surrogate
%% grad-free
%av = eta1'*xi_ref';
%Y_surr = get_polyfit_surr(av,f,eta1,1);
%% grad-based
%eta1g = W_grad(:,1);
%avg = eta1g'*xi_ref';
%Y_surrg = get_polyfit_surr(avg,f,eta1g,1);
%
%% Verification of both surrogates
%np = 1e4;
%Y_surr_data = zeros(np,1);
%Y_surr_datag = zeros(np,1);
%for i = 1:np
%  Y_surr_data(i,1) = Y_surr(xir_test(i,:)');
%  Y_surr_datag(i,1) = Y_surrg(xir_test(i,:)');
%end
%G = load('record_data_files/record_id_test_1e4.txt');
%G_mean = mean(G);
%G_var = var(G);
%Y_mean = mean(Y_surr_data);
%Y_var = var(Y_surr_data);
%Y_meang = mean(Y_surr_datag);
%Y_varg = var(Y_surr_datag);
%mean_var = zeros(3,2);
%mean_var(1,1) = G_mean; mean_var(1,2) = G_var;
%mean_var(2,1) = Y_mean; mean_var(2,2) = Y_var;
%mean_var(3,1) = Y_meang; mean_var(3,2) = Y_varg;
%save('mean_var.txt','mean_var','-ASCII');
%
%Nr = (sum((G-Y_surr_data).^2)).^(0.5);
%Dr = (sum((G).^2)).^(0.5);
%rn1 = Nr./Dr;
%Nrg = (sum((G-Y_surr_datag).^2)).^(0.5);
%Drg = (sum((G).^2)).^(0.5);
%rn1g = Nrg./Drg;
%rn = [rn1 rn1g];
%save('rel_norm_err_1e4.txt','rn','-ASCII');
%rel_err = [0.0850 0.0895 0.0904]; % fit analysis for the grad-free case
%figure;
%hold off;
%semilogy(rel_err,'ko','MarkerSize',5,'MarkerFaceColor','k');
%xlabel('$$\mathrm{Order~of~Polynomial~Fit}$$','interpreter','latex','fontsize',20);
%ylabel('$$\mathrm{\frac{\|G-\sim{G}\|_2}{\|G\|_2}}$$','interpreter','latex','fontsize',20);
%set(gca,'fontsize',18,'TickLabelInterpreter','latex');
%set(gcf,'color',[1,1,1]);
%box on;
%print -depsc rel_err.eps
%
% Load files for plotting

%W = load('eigv_gradfreeN1850.txt');
%lambda_loclin = load('lambda_gradfreeN1850.txt');
%xi_ref = load('xiref_gradfreeN1850.txt');
%f = load('f_gradfreeN1850.txt');
%
%% Activity Scores
%as_loclin = zeros(m,1);
%as_grad = zeros(m,1);
%for i = 1:m
%  for j=1:1
%    as_loclin(i) = as_loclin(i) + lambda_loclin(j).*(W(i,j).^2);
%    as_grad(i) = as_grad(i) + lam_grad(j).*(W_grad(i,j).^2);
%  end
%end
%as_loclin = as_loclin / sum(as_loclin);
%as_grad = as_grad / sum(as_grad);

% PLOTS

% Eigenvalue plot
%figure;
%semilogy(1:length(lambda_loclin),abs(lambda_loclin)./lambda_loclin(1),'-o','linewidth',2,'MarkerFaceColor','k');
%xlabel('$$\mathrm{Index~(i)}$$','interpreter','latex','fontsize',20);
%ylabel('$$\mathrm{Eigenvalue~(\lambda_i)}$$','interpreter','latex','fontsize',20);
%set(gca,'TickLabelInterpreter','Latex','fontsize', 18);
%%title('comparing dominant eigenvalues');
%box on;
%print -depsc eig.eps
%
%% Eigenvector plot
%figure;
%plot(W(:,1), '-bo', 'linewidth',2, 'markerfacecolor', 'k');
%xlabel('$$\mathrm{Index~(i)}$$','interpreter','latex','fontsize',20);
%ylabel('$$\mathrm{Eigenvector~Components~(w_i)}$$','interpreter','latex','fontsize',20);
%set(gca,'TickLabelInterpreter','Latex','fontsize', 18);
%%title('comparing dominant eigenvectors');
%box on;
%print -depsc eigv.eps
%
%% SSPs
%figure;
%hold on;
%g2 = -W(:,1)'*xi_ref';
%xig_ref = load('xiref_gradN40.txt');
%g2g = W_grad(:,1)'*xig_ref';
%g_loclin = f;
%np = length(g2);
%Y_surr_comp = zeros(np,1);
%Y_surr_compg = zeros(np,1);
%for i = 1:np
%  Y_surr_comp(i,1) = Y_surr(xi_ref(i,:)');
%  Y_surr_compg(i,1) = Y_surrg(xi_ref(i,:)');
%end
%plot(g2,g_loclin, 'ko', 'markerfacecolor', 'k','MarkerSize',5);
%hold on
%plot(g2g,f_grad, 'r*', 'markerfacecolor', 'r','MarkerSize',5);
%plot(g2, Y_surr_comp,'--k','LineWidth',2);
%plot(g2g, Y_surr_compg,'--r','LineWidth',2);
%xlabel('$$\mathrm{\hat{\mathbf{W}}_1^{\top}\mathbf{\xi}}$$',...
%'interpreter','latex','fontsize',20);
%ylabel('$$\mathrm{G(\hat{\mathbf{W}}_1^\top\mathbf{\xi})}$$','interpreter','latex','fontsize',20);
%leg = legend({'$\mathrm{G: Regression-based}$','$\mathrm{G: Perturbation-based}$','$\mathrm{\tilde{G}: Regression-based}$'...
%             ,'$\mathrm{\tilde{G}: Perturbation-based}$'});
%set(leg,'Interpreter','latex','fontsize',11,'location','SouthWest');
%set(gca,'TickLabelInterpreter','Latex','fontsize', 18);
%box on;
%print -depsc ssp_36D.eps
%%
%% Activity Scores
%subplot(1,2,1)
%bar([as_loclin(1:19) as_grad(1:19) gsa(1:19)]);
%bar([as_loclin(1:19) as_grad(1:19)]);

%xticklabels({'$\mathrm{A_1}$','$\mathrm{A_2}$','$\mathrm{A_3}$','$\mathrm{A_4}$',...
%             '$\mathrm{A_5}$','$\mathrm{A_6}$','$\mathrm{A_7}$','$\mathrm{A_8}$',...
%             '$\mathrm{A_9}$','$\mathrm{A_{10}}$','$\mathrm{A_{11}}$','$\mathrm{A_{12}}$',...
%             '$\mathrm{A_{13}}$','$\mathrm{A_{14}}$','$\mathrm{A_{15}}$','$\mathrm{A_{16}}$',...
%             '$\mathrm{A_{17}}$','$\mathrm{A_{18}}$','$\mathrm{A_{19}}$'});
%xticklabels({'$\mathrm{A_1}$','','','',...
%             '','','','',...
%             '$\mathrm{A_9}$','','','',...
%             '','','$\mathrm{A_{15}}$','',...
%             '$\mathrm{A_{17}}$','',''});
%ylabel('$$\mathrm{Activity~Scores~(\tilde{\nu_{i,r}})}$$','interpreter','latex','fontsize',14);
%ylim([0 0.35]);
%leg=legend('$\mathrm{Regression-based}$','$\mathrm{Perturbation-based}$','$\mathrm{Total~Sobol~Index}$');
%set(leg,'interpreter','latex','fontsize',11);
%set(gca,'xtick',1:19,'fontsize',12,'TickLabelInterpreter','latex');
%axis square;
%box on;
%
%subplot(1,2,2)
%bar([as_loclin(20:m) as_grad(20:m) gsa(20:m)]);
%bar([as_loclin(20:m) as_grad(20:m)]);
%xticklabels({'','','','',...
%             '','','','',...
%             '','$\mathrm{E_{a,15}}$','','',...
%             '','','$\mathrm{P_i}$','$\mathrm{T_i}$','$\mathrm{\Phi}$'});
%ylabel('$$\mathrm{Activity~Scores~(\nu_{i,r})}$$','interpreter','latex','fontsize',20);
%ylim([0 0.35]);
%leg=legend('$\mathrm{Regression-based}$','$\mathrm{Perturbation-based}$','$\mathrm{Total~Sobol~Index}$');
%set(leg,'interpreter','latex','fontsize',11,'location','NorthWest');
%set(gca,'xtick',1:17,'fontsize',12,'TickLabelInterpreter','latex');
%axis square;
%box on;
%print -depsc as_36D_new.eps
%
%% 2D Surrogate
%av1 = eta1'*xi_ref';
%av2 = eta2'*xi_ref';
%sf = fit([av1',av2'],f,'poly12');
%np = 1e4;
%Y_surr_2D = zeros(np,1);
%for i = 1:np
%  av1t = eta1'*xir_test(i,:)';
%  av2t = eta2'*xir_test(i,:)';
%  Y_surr_2D(i,1) = sf(av1t,av2t);
%end
%
%Nr = (sum((G-Y_surr_2D).^2)).^(0.5);
%Dr = (sum((G).^2)).^(0.5);
%rn2 = Nr./Dr;
%
%
%% PDF comparison
%step_d = (max(G) - min(G)).*(1e-4);
%step = (max(Y_surr_data) - min(Y_surr_data)).*(1e-4);
%stepg = (max(Y_surr_datag) - min(Y_surr_datag)).*(1e-4);
%%step2D = (max(Y_surr_2D) - min(Y_surr_2D)).*(1e-4);
%pts_surr_d = 0.9.*min(G):step_d:1.1.*max(G);
%pts_surr = 0.9.*min(Y_surr_data):step:1.1.*max(Y_surr_data);
%pts_surrg = 0.9.*min(Y_surr_datag):stepg:1.1.*max(Y_surr_datag);
%%pts_surr2D = 0.9.*min(Y_surr_2D):step2D:1.1.*max(Y_surr_2D);
%[density_surr_d,xmesh_surr_d] = ksdensity(G,pts_surr_d);
%[density_surr,xmesh_surr] = ksdensity(Y_surr_data,pts_surr);
%[density_surrg,xmesh_surrg] = ksdensity(Y_surr_datag,pts_surrg);
%%[density_surr2D,xmesh_surr2D] = ksdensity(Y_surr_2D,pts_surr2D);
%
%figure;
%hold on
%nbins = 30;
%%histogram(G,nbins,'Normalization','pdf','FaceAlpha',0.2);
%plot(xmesh_surr_d,density_surr_d,'Linewidth',2,'color','b');
%plot(xmesh_surr,density_surr,'--','Linewidth',2,'color','k');
%plot(xmesh_surrg,density_surrg,'-.','Linewidth',2,'color','r');
%%plot(xmesh_surr2D,density_surr2D,'--','Linewidth',2,'color','b');
%xlabel('$$\mathrm{Ignition~Delay (s)}$$','interpreter','latex','fontsize',20);
%ylabel('$$\mathrm{PDF}$$','interpreter','latex','fontsize',20);
%set(gca,'fontsize',18,'TickLabelInterpreter','latex');
%set(gcf,'color',[1,1,1]);
%leg = legend({'$\mathrm{Model}$','$\mathrm{\tilde{G}~(Regression-based)}$','$\mathrm{\tilde{G}~(Perturbation-based)}$'});
%set(leg,'Interpreter','latex','fontsize',14);
%box on;
%print -depsc pdf_comp_id_1e4.eps

end

function surr = get_polyfit_surr(y,G,eta,deg)
p = polyfit(y(:),G,deg);
surr = @(x)(polyval(p,eta'*x));
end

function [gsa_m,gsa_t] = sobol(m,xi)
W_grad = load('eigv_grad.txt');
xi_ref = zeros(size(xi,1),m);
xi_ref(:,1:24) = xi(:,1:24);
xi_ref(:,25:27) = xi(:,29:31);
xi_ref(:,28:33) = xi(:,33:38);
id = load('record_data_files/record_id_grad_free_n1000.txt');
f = zeros(1000,1);
f(:,1) = id(1:1000);
% grad-based
eta1g = W_grad(:,1);
avg = eta1g'*xi_ref';
Y_surrg = get_polyfit_surr(avg,f,eta1g,1);

N = 1e6;
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

end
