


A = [
    0   0.89287 1.54987
    0   0.69933 1.81970
    0   0.84022 1.28598
    0   0.79523 1.16012
    0   0.61266 1.12835
    0   0.39950 0.37942
    0   0.54807 1.66173
    0   0.50882 1.43175
    0   0.68840 1.58589
    0   0.59572 1.29311
    1   1.00787 1.09905
    1   1.23724 0.98834
    1   1.02175 0.67245
    1   0.88458 0.36003
    1   0.66582 1.22097
    1   1.24408 0.59735
    1   1.03421 0.88595
    1   1.66279 0.84183
];

gscatter(A(:,2),A(:,3),A(:,1))





%# generate data
num = 50;
X = [ mvnrnd([0.5 1.5], [0.025 0.03 ; 0.03 0.16], num) ; ...
      mvnrnd([1 1], [0.09 -0.01 ; -0.01 0.08], num)   ];
G = [1*ones(num,1) ; 2*ones(num,1)];

gscatter(X(:,1), X(:,2), G)
axis equal, hold on

for k=1:2
    %# indices of points in this group
    idx = ( G == k );

    %# substract mean
    Mu = mean( X(idx,:) );
    X0 = bsxfun(@minus, X(idx,:), Mu);

    %# eigen decomposition [sorted by eigen values]
    [V D] = eig( X0'*X0 ./ (sum(idx)-1) );     %#' cov(X0)
    [D order] = sort(diag(D), 'descend');
    D = diag(D);
    V = V(:, order);

    t = linspace(0,2*pi,100);
    e = [cos(t) ; sin(t)];        %# unit circle
    VV = V*sqrt(D);               %# scale eigenvectors
    e = bsxfun(@plus, VV*e, Mu'); %#' project circle back to orig space

    %# plot cov and major/minor axes
    plot(e(1,:), e(2,:), 'Color','k');
    %#quiver(Mu(1),Mu(2), VV(1,1),VV(2,1), 'Color','k')
    %#quiver(Mu(1),Mu(2), VV(1,2),VV(2,2), 'Color','k')
end






STD = 2;                     %# 2 standard deviations
conf = 2*normcdf(STD)-1;     %# covers around 95% of population
scale = chi2inv(conf,2);     %# inverse chi-squared with dof=#dimensions

Cov = cov(X0) * scale;
[V D] = eig(Cov);


