function [ Prox_lambda_nuclear_norm ] = Prox_lambda_nuclear_norm(A,lambda)
% iter=10000;
% i=1;
% A = rand(2,3);
% while i < iter
% % A = [1,2,3;
% %     -1,-2,-3];
%lambda = 1;


[s,u,v] = svd(A);
S_lambda = sign(u).*max(abs(u) - lambda,0);
Prox_lambda_nuclear_norm = s*S_lambda*v';
% B = rand(2,3);
% [s1,u1,v1] = svd(Prox_lambda);
% [s2,u2,v2] = svd(B);
% 
% nuclear_norm_Prox_lambda=sum(diag(u1));
% nuclear_norm_B=sum(diag(u2));
% A1 = 1/(2*lambda)*sum(diag((Prox_lambda-A)*(Prox_lambda-A)')) + nuclear_norm_Prox_lambda;
% B1 = 1/(2*lambda)*sum(diag((B-A)*(B-A)')) + nuclear_norm_B;
% if A1 > B1
%     break
% else 
%     i = i+1;
% end
%end

