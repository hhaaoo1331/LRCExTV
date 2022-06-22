function [ x,k,SNR,PSNR,SIM,tg,PSNR2]= Box_Dong_ADMM_constrained_impulsive(f_true,b,psf,iter,epsilon,known,rho)

% 2022.03.16
% This function implements ADMM method for  solving the total variation
% image denosing problem with impulse noise
% 
% \min_{x}  || x ||_{TV}+delta_C + \mu ||x||_{*}
% Kx \in C'
% equalently
%% Input 
% f_true ----------------- the original clear image
% b ---------------------- the noisy image 
% psf ---------------------- the blurring kernel
% iter -------------------- the total number of iterations
% epsilon ----------------- the stopping criterion
% known ------------------- the location of noise free images
% rho --------------------- the parameter of ADMM

%% Output 
% x ----------------------- the restored images
% k ----------------------- the number of iterations
% SNR --------------------- the SNR values
% PSNR -------------------- the PSNR values
% SIM  -------------------- the SSIM index
% tg  --------------------- the CPU computation time
% PSNR2 ------------------- the PSNR sequences
%
%

	[H,W]=size(f_true);
	opD = @(x) cat(3,[diff(x,1,1);zeros(1,W)],[diff(x,1,2) zeros(H,1)]);
	opDadj = @(u) [-u(1,:,1);-diff(u(1:end-1,:,1),1,1);u(end-1,:,1)]+...
		[-u(:,1,2) -diff(u(:,1:end-1,2),1,2) u(:,end-1,2)];


    
    % data fidelity
A_dir  = @(x) imfilter(x, psf,'circular');
A_adj  = @(x) imfilter(x, rot90(psf,2),'circular');  % WARNING: 'psf' must be a (2n+1)-by-(2n+1) matrix

 [m,n]=size(f_true);
% y =zeros(m,2*n);

 % C.eigsDtD = abs(psf2otf([1,-1],[m,n])).^2 + abs(psf2otf([1;-1],[m,n])).^2;

 % initial value
y = zeros([m,n,2]);
w = zeros(m,n);
z = zeros(m,n);
x = zeros(m,n);
q = zeros(m,n);

lam1 = zeros([m,n,2]);
lam2 = zeros(m,n);

lam4 = zeros(m,n);


lambda = 0.01;

% rho = 10;


SNR = [];
PSNR = [];
SIM = [];

%  SNR(1)    = 20*log10(norm(f_true(:))/norm(f_true(:)-f(:)));
%  % PSNR(1)   = 20*log10(sqrt(m*n)*255/norm(f_true(:)-f(:))) ; % [0,255]
%    PSNR(1)   = 20*log10(sqrt(m*n)/norm(f_true(:)-f(:))) ; % [0,1]
%   SIM(1) = ssim(f,f_true);

k =1;

tic;
t1 = clock;
while k <= iter
    
    
    x_update = x - lambda*( opDadj(lam1) + rho*opDadj(opD(x)-y) + lam2 + rho*(x-z)  + A_adj(lam4)  + rho*(A_adj(A_dir(x)-q)) );
    
    yp = opD(x_update) + 1/rho * lam1;
    y_update = max(abs(yp)-1/rho,0).*sign(yp);
    
    zp = x_update + 1/rho * lam2;
    z_update = proj_bound(zp,0,1);
    

    
    qp = A_dir(x_update) + 1/rho * lam4;
    q_update = b + qp.*(1-known);
    
    
   lam1_update = lam1 + rho*(opD(x_update)-y_update);
   
   lam2_update = lam2 + rho*(x_update - z_update);
   
   
   lam4_update = lam4 + rho*(A_dir(x_update)-q_update);

    
   
 
  %   x_update = omega1*x1_update + omega2*x2_update + omega3 * x3_update;
   
       %   e = omega1*e1 + omega2*e2 + omega3*e3;
       SNR    = 20*log10(norm(f_true(:))/norm(f_true(:)- x_update(:)));
%      PSNR(k+1)   = 20*log10(sqrt(m*n)*255/norm(f_true(:)-e(:))) ;
        PSNR = psnr(x_update,f_true); 
        PSNR2(k) = psnr(x_update,f_true); 
        SIM = ssim( x_update,f_true);
 %       fval(k) = 0.5*norm(x_update-f)^2 + lambda*norm(D*x_update,1);
    
tg(k) = etime(clock,t1);
 
 
    if   norm(x_update(:)-x(:))/norm(x(:)) <= epsilon
        break;
    else
        x = x_update;
    lam1 = lam1_update;
    lam2 = lam2_update;

    lam4 = lam4_update;
    
    y = y_update;
    z = z_update;

    q = q_update;
    
       k = k+1;
    end
    
 
% u=f+lamda*div(y1_update,y2_update,1);

% mse(k) = norm(f_true - u,'fro')^2/(m*n);

end
% tg = toc;
end

