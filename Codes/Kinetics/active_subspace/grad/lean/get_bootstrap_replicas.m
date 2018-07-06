function get_bootstrap_replicas(Mboot, df, W1, lambda1) 
%
%
%  This is the bootstrapping routine; first need to run main script file
%  Alg1_1.m
% 
%  Input: 
%     Mboot --- number of bootstrap replicates to perform
%     dGdx --- gradient
%     lambda1 --- eigenvalues
%     W1 ---- eigenvectors
%
%




m = 19; %dimension of the input parameter space
k = m + 1;

nsamples = size(df, 1);

all_lambda = zeros(m, Mboot);

K = 0;
for p = 1 : Mboot
       disp(p)
       % pick random numbers between 1 and the number of samples 
       j = randi([1, nsamples], nsamples, 1);

        for k = 1 : nsamples

            gfk = df(j(k),:)';

            K = K + gfk * gfk';
        end

        K = K / nsamples;

        %compute eigenvalues and eigenvectors
        [W D] = eig(K);
        [lambda idx] =sort(diag(D), 'descend');
        
        %saves all the values of lambda for each bootstrap replicate
        all_lambda(:, p) = lambda;

        W = W(:,idx);

        %first column of ev matrix for each replicate

        %eta1new(:,p) = W(:,1);
        
        %compute mean, lower and upper bound of the eigenvalues
        for i=1:m
           z(i) = mean(all_lambda(i,:));
           lower(i) = min(all_lambda(i,:));
           upper(i) = max(all_lambda(i,:));
        end

        %calculating the error between subspaces
        for j=1:8
          errornew(p,j)=norm(W1(:,j)*(W1(:,j)')-W(:,j)*(W(:,j)'));
        end
        Z1 = mean(errornew);
        Z2 = max(errornew);
        Z3 = min(errornew);

	%produce a plot of the eigenvalues and bootstrap intervals
	figure(1);
	y1 = lower ./ lower(1);
	y2 = upper ./ upper(1);
	leng = 1 : 19;
	X = [leng, fliplr(leng)];
	Y = [y1, fliplr(y2)];  
	fill(X, Y, 1, 'facecolor','blue', ...
              'edgecolor','none', ...
              'facealpha', 0.3); 
        set(gca,'YScale','log')
        hold on
end

% plot the lambdas obtained from active subspace calculation
semilogy(lambda1./lambda1(1),'-o','Color','k');
xlabel('Index');
ylabel('Eigenvalues with bootstrap intervals');
title('Estimated Eigenvalues and bootstrap intervals');
 
% subspace errors
figure(2);
x1 = Z3;
x2 = Z2;
leng = 1 : 8;
X1 = [leng, fliplr(leng)];
Y1 = [x1, fliplr(x2)];  
fill(X1, Y1, 1,'facecolor','blue', ...
      'edgecolor', 'none', ...
      'facealpha', 0.3); 
set(gca,'YScale','log')
ylim([1e-5 1.1]);
hold on
semilogy(Z1,'-o', 'Color', 'k');
xlabel('Subspace Dimension');
ylabel('Subspace Error');