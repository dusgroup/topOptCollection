%{
This code is written by 
Gao J, Luo Z, Xia L, Gao L.
Concurrent topology optimization of multiscale composite structures in Matlab,
Structural and Multidisciplinary Optimization, 2019, 60: 2621-2651.
DOI: 10.1007/s00158-019-02323-6

Commented by
DU Run and his team
%}
function ConTop2D(Macro_struct,Micro_struct,penal,rmin)
tic; % 开始计时
% USER-DEFINED DESIGN AND STRUCTURAL PARMETERS
Macro_struct = [10, 5, 100, 50, 0.4];
Micro_struct = [0.1, 0.1, 50, 50, 0.5]; penal = 3; rmin = 2;
maxloop = 200; E0 = 1; Emin = 1e-9; nu = 0.3;
Macro.length = Macro_struct(1); Macro.width = Macro_struct(2);
Micro.length = Micro_struct(1); Micro.width = Micro_struct(2);
Macro.nelx = Macro_struct(3); Macro.nely = Macro_struct(4);
Micro.nelx = Micro_struct(3); Micro.nely = Micro_struct(4);
Macro.Vol  = Macro_struct(5); Micro.Vol  = Micro_struct(5);
Macro.Elex = Macro.length/Macro.nelx; Macro.Eley = Macro.width/Macro.nely;
Macro.nele = Macro.nelx*Macro.nely; Micro.nele = Micro.nelx*Micro.nely;
Macro.ndof =2*(Macro.nelx+1)*(Macro.nely+1);
% PREPRE FINITE ELEMENT ANALYSIS
[load_x,load_y] = meshgrid(Macro.nelx, Macro.nely/2);
loadnid = load_x*(Macro.nely+1)+(Macro.nely+1-load_y);
F = sparse(2*loadnid(:),1,-1,2*(Macro.nelx+1)*(Macro.nely+1),1);
U = zeros(Macro.ndof, 1);
[fixed_x, fixed_y] = meshgrid(0,0:Macro.nely);
fixednid = fixed_x*(Macro.nely+1)+(Macro.nely+1-fixed_y);
fixeddofs = [2*fixednid(:);2*fixednid(:)-1];
% F = sparse(2,1,-1,2*(Macro.nely+1)*(Macro.nelx+1),1);
% fixeddofs = union(1:2:2*(Macro.nely+1),2*(Macro.nelx+1)*(Macro.nely+1));
freedofs = setdiff(1:Macro.ndof,fixeddofs);
nodenrs = reshape(1:(Macro.nely+1)*(Macro.nelx+1),1+Macro.nely,1+Macro.nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,Macro.nele,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*Macro.nely+[2 3 0 1] -2 -1],Macro.nele,1);
iK = reshape(kron(edofMat,ones(8,1))',64*Macro.nele,1);
jK = reshape(kron(edofMat,ones(1,8))',64*Macro.nele,1);
% PREPARE FILTER
[Macro.H,Macro.Hs] = filtering2d(Macro.nelx,Macro.nely,Macro.nele,rmin);
[Micro.H,Micro.Hs] = filtering2d(Micro.nelx,Micro.nely,Micro.nele,rmin);
% INITIALIZATIONS AT TWO SCALES 
Macro.x = repmat(Macro.Vol,Macro.nely,Macro.nelx);
Micro.x = ones(Micro.nely, Micro.nelx);
for i = 1:Micro.nelx
   for j= 1:Micro.nely
       if sqrt((i-Micro.nelx/2-0.5)^2+(j-Micro.nely/2-0.5)^2)<min(Micro.nelx,Micro.nely)/3
           Micro.x(j,i)= 0;
       end
   end
end
beta =1;
Macro.xTilde = Macro.x; Micro.xTilde = Micro.x;
Macro.xPhys = 1-exp(-beta*Macro.xTilde)+Macro.xTilde*exp(-beta);
Micro.xPhys = 1-exp(-beta*Micro.xTilde)+Micro.xTilde*exp(-beta);
loopbeta = 0; loop = 0; Macro.change = 1; Micro.change = 1;
while loop < maxloop || Macro.change > 0.01 || Micro.change >0.01
    loop = loop+1; loopbeta=loopbeta+1;
    % FE-ANALYSIS AT TWO SCALES
    [DH,dDH] = EBHM2D(Micro.xPhys,Micro.length,Micro.width,E0,Emin,nu,penal);
    Ke = elementMatVec2D(Macro.Elex/2,Macro.Eley/2,DH);
    sK = reshape(Ke(:)*(Emin+Macro.xPhys(:)'.^penal*(1-Emin)),64*Macro.nele,1);
    K = sparse(iK,jK,sK); K=(K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*Ke).*U(edofMat),2),Macro.nely,Macro.nelx);
    c = sum(sum((Emin+Macro.xPhys.^penal*(1-Emin)).*ce));
    Macro.dc = -penal*(1-Emin)*Macro.xPhys.^(penal-1).*ce;
    Macro.dv = ones(Macro.nely,Macro.nelx);
    Micro.dc = zeros(Micro.nely,Micro.nelx);
    for i = 1:Micro.nele
        dDHe = [dDH{1,1}(i) dDH{1,2}(i) dDH{1,3}(i);
                dDH{2,1}(i) dDH{2,2}(i) dDH{2,3}(i);
                dDH{3,1}(i) dDH{3,2}(i) dDH{3,3}(i)];
        [dKE] = elementMatVec2D(Macro.Elex,Macro.Eley,dDHe);
        dce = reshape(sum((U(edofMat)*dKE).*U(edofMat),2),Macro.nely,Macro.nelx);
        Micro.dc(i) = -sum(sum((Emin+Macro.xPhys.^penal*(1-Emin)).*dce));
    end
    Micro.dv = ones(Micro.nely,Micro.nelx);
    %FILTERING ANDMODIFICATION FSESITIVITIES
    Macro.dx = beta*exp(-beta*Macro.xTilde)+exp(-beta);
    Micro.dx = beta*exp(-beta*Micro.xTilde)+exp(-beta);
    Macro.dc(:) = Macro.H*(Macro.dc(:).*Macro.dx(:)./Macro.Hs);
    Macro.dv(:) = Macro.H*(Macro.dv(:).*Macro.dx(:)./Macro.Hs);
    Micro.dc(:) = Micro.H*(Micro.dc(:).*Micro.dx(:)./Micro.Hs);
    Micro.dv(:) = Micro.H*(Micro.dv(:).*Micro.dx(:)./Micro.Hs);
    % OPTIMALITY CRITERIA UPDATE FOR MACRO AND MICRO ELELMENT DENSITIES
    [Macro.x,Macro.xPhys,Macro.change] = OC(Macro.x, Macro.dc, Macro.dv, Macro.H, Macro.Hs, Macro.Vol, Macro.nele, 0.2, beta);
    [Micro.x,Micro.xPhys,Micro.change] = OC(Micro.x, Micro.dc, Micro.dv, Micro.H, Micro.Hs, Micro.Vol, Micro.nele, 0.2, beta);
    Macro.xPhys = reshape(Macro.xPhys, Macro.nely, Macro.nelx); 
    Micro.xPhys = reshape(Micro.xPhys, Micro.nely, Micro.nelx);
    % PRINT RESULTS
    fprintf('It.:%5i Obj.:%11.4f Macro_Vol.:%7.3f Micro_Vol.:%7.3f Macro_ch.:%7.3f Micro_ch.:%7.3fn', ...
        loop,c,mean(Macro.xPhys(:)),mean(Micro.xPhys(:)), Macro.change, Micro.change);
    figure(1); colormap(gray); imagesc(1-Macro.xPhys); caxis([0 1]); axis equal; axis off; drawnow;
    figure(2); colormap(gray); imagesc(1-Micro.xPhys); caxis([0 1]); axis equal; axis off; drawnow;
    % UPDATE HEAVISIDE REURATIONRTR
%     if beta <512 && (loopbeta >= 50 || Macro.change <= 0.01 || Micro.change <= 0.01)
%         beta=2*beta;loopbeta=0;Macro.change=1;Micro.change=1;
%         fprintf('Parameter beta increased to %g.\n',beta);
%     end
end
elapsedTime = toc; % 结束计时并显示时间  
fprintf('程序运行时间：%.6f 秒\n', elapsedTime);
end

function Ke = elementMatVec2D(a, b, DH)
GaussNodes = [-1/sqrt(3); 1/sqrt(3)]; GaussWeigh = [1 1];
L = [1 0 0 0; 0 0 0 1; 0 1 1 0]; Ke = zeros(8,8);
for i = 1:length(GaussNodes)
    for j = 1:length(GaussNodes)
        GN_x = GaussNodes(i); GN_y = GaussNodes(j);
        dN_x = 1/4*[-(1-GN_x) (1-GN_x) (1+GN_x) -(1+GN_x)];
        dN_y = 1/4*[-(1-GN_y) -(1+GN_y) (1+GN_y) (1-GN_y)];
        J = [dN_x; dN_y]*[ -a a a -a;-b -b b b]';
        G = [inv(J) zeros(size(J)); zeros(size(J)) inv(J)];
        dN(1,1:2:8) = dN_x; dN(2,1:2:8) = dN_y;
        dN(3,2:2:8) = dN_x; dN(4,2:2:8) = dN_y;
        Be = L*G*dN;
        Ke = Ke + GaussWeigh(i)*GaussWeigh(j)*det(J)*Be'*DH*Be;
    end
end
end

function [H,Hs] = filtering2d(nelx,nely,nele,rmin)
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
end

function [DH,dDH] = EBHM2D(xPhys,lx,ly,E0,Emin,nu,penal)
[nely,nelx] = size(xPhys); cellVolume = lx*ly;
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);

e0 = eye(3);
ufixed=zeros(8,3);
U=zeros(2*(nely+1)*(nelx+1),3);
alldofs=(1:2*(nely+1)*(nelx+1));
% 此处由5.1节详细描述，n1即为4个角处的节点编号，
% d1即通过乘2减1和乘2的方法计算n1节点对应的自由度
n1 = [nodenrs(end,[1,end]),nodenrs(1,[end,1])];
d1 = reshape([(2*n1-1);2*n1],1,8);
% 同样n3，d3用同样的方法表示左边界和下边界的节点和节点自由度
n3 = [nodenrs(2:end-1,1)',nodenrs(end,2:end-1)];
d3 = reshape([(2*n3-1);2*n3],1,2*(nelx+nely-2));
% n4，d4则是右边界和上边界的节点和节点自由度
n4 = [nodenrs(2:end-1,end)',nodenrs(1,2:end-1)];
d4 = reshape([(2*n4-1);2*n4],1,2*(nelx+nely-2));
% d2为剩下的中间点自由度
d2 = setdiff(alldofs,[d1,d3,d4]);
% 此处计算对应的3个单元实验应变场，计算的结果组装成的wfixed即为公式（13）（19）中的w
for j= 1:3
    ufixed(3:4,j)=[e0(1,j),e0(3,j)/2;e0(3,j)/2,e0(2,j)]*[nelx;0];
    ufixed(7:8,j)=[e0(1,j),e0(3,j)/2;e0(3,j)/2,e0(2,j)]*[0;nely];
    ufixed(5:6,j)= ufixed(3:4,j)+ufixed(7:8,j);
end
wfixed =[repmat(ufixed(3:4,:),nely-1,1);repmat(ufixed(7:8,:),nelx-1,1)];


qe = cell(3,3);
Q = zeros(3,3);
dQ = cell(3,3);

sK= reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
K=sparse(iK,jK,sK);K=(K+K')/2;
% 此处的计算为4.2节均质化方程的数值解中式（19）的表示
Kr =[K(d2,d2),K(d2,d3)+K(d2,d4);K(d3,d2)+K(d4,d2),K(d3,d3)+K(d4,d3)+K(d3,d4)+K(d4,d4)];
U(d1,:)= ufixed;
U([d2,d3],:)=Kr\(-[K(d2,d1);K(d3,d1)+K(d4,d1)]*ufixed-[K(d2,d4);K(d3,d4)+K(d4,d4)]*wfixed);
U(d4,:)=U(d3,:)+wfixed;
%% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
% 此处的qe为式（9）中的元素互能和，此处将通过遍历循环来将3*3的小单元组合，见式（7）
% Q为式（8）中的Q
for i = 1:3
    for j= 1:3
      U1 = U(:,i);U2 = U(:,j);
      qe{i,j}= reshape(sum((U1(edofMat)*KE).*U2(edofMat),2),nely,nelx);
      DH(i,j) = 1/cellVolume*sum(sum((Emin+xPhys.^penal*(E0-Emin)).*qe{i,j}));
      dDH{i,j} = 1/cellVolume*(penal*(E0-Emin)*xPhys.^(penal-1).*qe{i,j});
    end
end
end

function [x,xPhys,change] = OC(x,dc, dv, H, Hs, Vol, nele, changes, beta)
l1=0;l2=1e9;move =0.2;
while(l2-l1>1e-9)
lmid =0.5*(l2+l1);
xnew = max(0,max(x-move,min(1,min(x+move ,x.*sqrt(-dc./dv/lmid)))));
% if ft == 1
% xPhys = xnew;
% elseif ft == 2
xPhys(:)=(H*xnew(:))./Hs;
% end
  if mean(xPhys(:))>Vol,l1= lmid;
    else 
      l2 = lmid;
  end
end
change=max(abs(xnew(:)-x(:)));
x=xnew;
end