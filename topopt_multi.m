%95 LINES MATLAB CODE MULTIPHASE MINIMUM COMPLIANCE TOPOLOGY OPTIMIZATION
function topopt_multi
clear all; clc; nx = 96; ny = 48; q = 3; p = 4;
epsi = 1; alpha = 0.5; zeta = 0.125; maxit = 1000;
nu = 0.3; e = [4 2 1 1e-9]'; v = [0.15 0.15 0.15 0.55]';
%INIT DESIGN VECTOR AND RHS OF VOLUME CONSTRAINTS
x = ones(nx*ny,p)*repmat(v',length(v),1)/length(v); sx = nx*ny*v';
%PREPARE FINITE ELEMENT TOOLS FOR ELASTICITY
A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4];
B12 = [ 2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nx)*(1+ny),1+ny,1+nx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nx*ny,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*ny+[2 3 0 1] -2 -1],nx*ny,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nx*ny,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nx*ny,1);
%PREPARE FINITE ELEMENT TOOLS FOR HELMHOLTZ SOLVER
D = alpha*epsi*zeta;
KEF = D*[4 -1 -2 -1; -1 4 -1 -2; -2 -1 4 -1; -1 -2 -1 4]/6 + ...
    [4 2 1 2; 2 4 2 1; 1 2 4 2; 2 1 2 4]/36;
edofVecF = reshape(nodenrs(1:end-1,1:end-1),nx*ny,1);
edofMatF = repmat(edofVecF,1,4)+repmat([0 ny+[1:2] 1],nx*ny,1);
iKF = reshape(kron(edofMatF,ones(4,1))',16*nx*ny,1);
jKF = reshape(kron(edofMatF,ones(1,4))',16*nx*ny,1);
sKF = reshape(KEF(:)*ones(1,nx*ny),16*nx*ny,1);
KF = sparse(iKF,jKF,sKF); LF = chol(KF,'lower');
iTF = reshape(edofMatF,4*nx*ny,1); jTF = reshape(repmat([1:nx*ny],4,1)',4*nx*ny,1);
sTF = repmat(1/4,4*nx*ny,1);
TF = sparse(iTF,jTF,sTF);
%DEFINE LOADS AND SUPPORTS FOR CANTILEVER BEAM
F = sparse(2*(ny+1)*(nx+1),1,-1,2*(ny+1)*(nx+1),1);
fixeddofs = [1:2*(ny+1)];
U = zeros(2*(ny+1)*(nx+1),1);
alldofs = [1:2*(ny+1)*(nx+1)];
freedofs = setdiff(alldofs,fixeddofs);
% OPTIMIZATION LOOP
for iter = 1:maxit
    [obj objc g] = eval_fg (nx,ny,epsi,zeta,p,q,e,KE,U,edofMat,freedofs,iK,jK,F,x);
    x = x + alpha * (proj (sx, x - g) - x);
    for i = 1:p-1
        x(:,i) = TF'*(LF'\(LF\(TF*x(:,i))));
    end
    x(:,p) = 1 - sum(x(:,1:p-1),2);
    fprintf('It:%5i c:%8.4f obj:%8.4f\n',iter, objc, obj);
    % VISUALIZATION OF CURRENT TOPOLOGY
    if mod(iter,10) == 0
        I = make_bitmap (p,nx,ny,x); image(I), axis image off, drawnow;
    end
end

%COMPUTE OBJECTIVE FUNCTION AND ITS GRADIENT
function [f fc g] = eval_fg (nx,ny,epsi,zeta,p,q,e,KE,U,edofMat,freedofs,iK,jK,F,x)
[fp, gp] = eval_fg_perimeter (epsi,nx,ny,p,x);
[fc, gc] = eval_fg_compliance (nx,ny,p,q,e,KE,U,edofMat,freedofs,iK,jK,F,x);
f = zeta * fp + fc;
g = zeta * gp + gc;
function [f, g] = eval_fg_compliance (nx,ny,p,q,e,KE,U,edofMat,freedofs,iK,jK,F,x)
E = (x.^q)*e;
sK = reshape(KE(:)*E(:)',64*nx*ny,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs) = K(freedofs,freedofs)\F(freedofs);
ce = sum((U(edofMat)*KE).*U(edofMat),2);
f = sum(sum(E.*ce));
g = (-0.5*q*x.^(q-1)).*(ce*e');
g = g - repmat(g(:,p),1,p);
function [f g] = eval_fg_perimeter (epsi,nx,ny,p,x)
i = 2:ny+1; j = 2:nx+1; im = i-1; ip = i+1; jm = j-1; jp = j+1;
v = zeros(ny+2,nx+2,p-1); v(2:ny+1,2:nx+1,:) = reshape(x(:,1:p-1),ny,nx,p-1);
v(1,:,:) = v(2,:,:); v(ny+2,:,:) = v(ny+1,:,:);
v(:,1,:) = v(:,2,:); v(:,nx+2,:) = v(:,nx+1,:);
f1 = (v(ip,j,:)-v(im,j,:)).^2 + (v(i,jp,:)-v(i,jm,:)).^2;
f2 = (x.^2).*(x - 1).^2;
f = (epsi/8)*sum(f1(:))+(1./epsi)*sum(f2(:));
g = (4*x.^3 - 6*x.^2 + 2*x)/epsi; g(:,p) = 0;

% ORTHOGONAL PROJECTION OF A ATRIAL POINT ONTO THE CONTROL SPACE
function x = proj (b, x)
m = size(x,2)-1; n = size(x,1); res = inf; mu = zeros(m,1);
while res > 1.e-10
    x_old = x; diff = 1;
    while diff > 1.e-4
        lambda_u = (1/m)*max(sum(x(:,1:m), 2)- 1 - sum(mu), 0);
        lambda_l = (1/m)*min(sum(x(:,1:m), 2) - sum(mu), 0);
        lambda = lambda_u + lambda_l;
        mu_new = (sum(x(:,1:m)-repmat(lambda,1,m))-b(1:m))'/n;
        diff = norm(mu - mu_new, inf); mu = mu_new;
    end
    x(:,1:m) = max(0, min (1, x(:,1:m) - repmat(lambda,1,m) - repmat(mu',n,1)));
    res = norm(x(:) - x_old(:), inf);
end

% MAKE BITMAP IMAGE OF MULTIPHASE TOPOLOGY
function I = make_bitmap(p,nx,ny,x)
color = [0 0 0; 1 0 0; 0 0 1; 0 1 0; 1 1 0; 1 .1 .7; 0.44 0.86 0.86; 1 1 1];
color(p,:) = [1 1 1];
I = imresize(reshape(min(max(x*color(1:p,:),0),1),ny,nx,3),10,'bilinear');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is from the paper                                               %
% "Multimaterial topology optimization by volume constrained               %
%  Allen–Cahn system and regularized projected steepest descent method"    %
% by Rouhollah Tavakoli   (2014),                                          %
% Computer methods in applied mechanics and engineering,                   %
% Vol 276, pp.534-565.  DOI:10.1016/j.cma.2014.04.005                      %
%                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%