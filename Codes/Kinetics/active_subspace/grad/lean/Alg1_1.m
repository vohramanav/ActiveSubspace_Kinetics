function [dGdx, lambda1, W1] = Alg1_1(dim);

%
%
% This script performs active subspace testing using Algorithm 1.1 
%
%
%   Input: 
%     dim --- dimension of input parameter space
%
%  Output:
%     dGdx --- gradient
%     lambda1 --- eigenvalues
%     W1 ---- eigenvectors

rng(100);
set(0,'DefaultFigureVisible','off');

dX = zeros(dim,1); 
L = zeros(dim,1); U = zeros(dim,1);
id_data = load('record_id_as.txt');
params = load('pts_as.txt');

% Nominal values of the parameters
nom = [1.915e14,5.080e04,2.160e08,1.230e04,4.577e19,6.165e15,4.714e18,2.240e22,6.170e19,...
       6.630e13,1.690e14,1.810e13,1.450e16,3.020e12,1.202e17,1.000e13,4.820e13,9.550e06,...
       7.000e12];

L(:,1) = 0.9.*nom(:,1); U(:,1) = 1.1.*nom(:,1);

% Project params in [-1,1]

nrows = size(params,1); % number of rows
ncols = size(params,2); % number of cols

xp = zeros(nrows,ncols);

for i = 1:nrows
  for j = 1:ncols
    xp(i,j) = 2.0.*(params(i,j)-L(j))./(U(j)-L(j)) - 1.0;
  end
end

nsams = 60;


for j = 1:dim
  dX(j) = 1e-3.*(U(j,1)-L(j,1));
end

G = zeros(nsams,1); Gdx = zeros(nsams,dim);
G(:,1) = id_data(1:nsams,1);
dGdx = zeros(nsams,dim);

for j = 1:dim
  Gdx(:,j) = id_data(j*nsams+1:(j+1)*nsams,1);
  dGdx(:,j) = (Gdx(:,j) - G(:,1))./dX(j);
end

% Computing the C matrix (Eq. 1.2)

C = zeros(dim,dim);
grad_f = zeros(dim,1);

for i = 1:nsams
  grad_f(:,1) = dGdx(i,:);
  C = C + grad_f*transpose(grad_f);
end

C = (1./nsams).*(C);
[W1,D] = eig(C);

[lambda1, idx] = sort(diag(D), 'descend');
W1 = W1(:,idx);

