rng(123);
m1 = 19; % number of uncertain A's
m2 = 14; % number of uncertain E's
m3 = 3; % uncertain P,T,Phi
m = m1+m2+m3;
N = 105;

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

% locations of 0E's: 6(25),7(26),8(27),9(28),13(32)

nom1 = log(nom1); % perturb log(A) instead of A 

L1(1,:) = 0.98.*nom1(1,:); U1(1,:) = 1.02.*nom1(1,:);
L2(1,:) = 0.98.*nom2(1,:); U2(1,:) = 1.02.*nom2(1,:);
Lp(1,:) = 0.98.*nomp(1,:); Up(1,:) = 1.02.*nomp(1,:);

% draw samples in [-1,1]
xi = -1 + 2 * lhsdesign(N,m);
pts_xi = zeros(N*(m+1),m);
pts_xi(1:N,:) = xi;

% perturb the points
dxi = 2e-5;
for i = 1:m
  pts_xi(i*N+1:(i+1)*N,:) = pts_xi(1:N,:);
  pts_xi(i*N+1:(i+1)*N,i) = pts_xi(i*N+1:(i+1)*N,i) + dxi;
end

gen_samples(N,pts_xi,L1,L2,Lp,U1,U2,Up);
%acsub(N,dxi,m,xi);

%% Project points to the physical space
function gs = gen_samples(N,pts_xi,L1,L2,Lp,U1,U2,Up)
pts_x = zeros(size(pts_xi,1),size(pts_xi,2)+6); % 5 E's with 0 nom value + 1 for MO2

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
save('pts_gradN105.txt','pts_x','-ASCII');
end

function as = acsub(N,dxi,m,xi)
xi_ref = xi;
save('xi_gradN50.txt','xi_ref','-ASCII');

id = load('data_files_p2/record_id_gradN50.txt');
G = zeros(N,1); Gdxi = zeros(N,m);
G(:,1) = id(1:N,1);

dGdxi = zeros(N,m);
%
for j = 1:m
  Gdxi(:,j) = id(j*N+1:(j+1)*N,1);
  dGdxi(:,j) = (Gdxi(:,j) - G(:,1))./dxi;
end
%
C = zeros(m,m);
grad_f = zeros(m,1);

for i = 1:N
  grad_f(:,1) = dGdxi(i,:);
  C = C + grad_f*transpose(grad_f);
end

C = (1./N).*(C);
[V,D] = eig(C);

[lambda_grad, idx] = sort(diag(D), 'descend');
V = V(:,idx);
save('lambda_gradN50.txt','lambda_grad','-ASCII');
save('eigv_gradN50.txt','V','-ASCII');
%
%% Computing activity scores
%
%as = zeros(m,1);
%
%for i = 1:m
%  for j=1:1
%    as(i) = as(i) + lambda_grad(j).*(V(i,j).^2);
%  end
%end
%
%as = as./sum(as);
%
% Eigenvalue Plot
%figure
%semilogy(1:length(lambda_grad),abs(lambda_grad)./lambda_grad(1),'ko','linewidth',2,'MarkerFaceColor','k');
%xlabel('$$\mathrm{Index~(i)}$$','interpreter','latex','fontsize',20);
%ylabel('$$\mathrm{Eigenvalue~(\lambda_i)}$$','interpreter','latex','fontsize',20);
%set(gca,'TickLabelInterpreter','Latex','fontsize', 18);
%%title('comparing dominant eigenvalues');
%box on;
%print -depsc eig_33D.eps

%
eta1 = -V(:,1);
%eta2 = V(:,2);
save('f_gradN50.txt','G','-ASCII');
%
% univariate
%
figure;
g1 = eta1'*xi_ref';
plot(g1, G, 'ko', 'markerfacecolor', 'k')
set(gca, 'fontsize', 20);
xlabel('<eta1, x>');
ylabel('f(x)');
title('grad based SSP');
print -dpng ssp_grad.png

end


