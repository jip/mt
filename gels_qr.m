m = 9
n = 7
A = [1 6.37526 40.6439 259.115 1651.93 10531.4 67140.7;1 _1.27002 _72.1821 279.114 4734.14 _33072.3 _272987;0 1 _2.54003 _216.546 1116.46 23670.7 _198434;0 8.5904 _21.8199 _592.362 3150.01 36667.6 _330673;0 0 17.1808 _65.4597 _2369.45 15750.1 220006;1 _6.83522 46.7203 _319.344 2182.78 _14919.8 101980;1 0 0 0 0 0 0;0 1 0 0 0 0 0;0 0 2 0 0 0 0]
b = [587.136 _0.188615 _0.188615 0.208058 0.208058 0.00107523 1 1 1]'

Q,R,perm] = qr(A,0)
i = 1; while (abs(R(i,i)) > abs(R(1,1))*m*eps ) i=i+1; if (i==n+1), break; end; end;
rankA =i-1
R = R(1:rankA,1:n)
Q = Q(:,1:rankA)
if (rankA < n), 
   [Z,R] = qr( R(1:rankA,1:n)', 0 ); R = R'; Z=Z';
else
   Z = eye(n);
end
z = Q'*b
x(perm,1) = Z' *(R \ z)
mse = (norm(b)^2 - norm(z)^2) / (m-rankA)
stdx = zeros(n,1)
for i=1:n, 
   stdx(perm(i)) = norm( R\Z(:,i) ) * sqrt( mse ); 
end

fprintf('rankA = %d\n',rankA);
fprintf('x    = [ %+f %+f %+f ]''\n',x(1),x(2),x(3));
fprintf('mse  = %f\n',mse);
fprintf('stdx = [ %+f %+f %+f ]''\n',stdx(1),stdx(2),stdx(3));
fprintf('\n');
   