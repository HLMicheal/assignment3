% ELEC4700 Assignment 3
% By Huanyu Liu 100986552
% Part2

clear
clc
L=40;
W=60;
sigma1=1;
sigma2=0.01;
k=L*W;
G=sparse(k,k);
Z = zeros(k, 1);
Vo=1; % voltage drop across the whole area = the beginning voltage
S = zeros(L, W);    % sigma
for x = 1 : L
    for y = 1 : W
        if x >= 0.4*L && x <= 0.6*L && (y <= 0.4*W || y >= 0.6*W)
            % area in the blocks
            S(x, y) = sigma2;
        else
            % area outside the blocks
            S(x, y) = sigma1;
        end
    end
end
for i = 1:L
    for j = 1:W
        n = j + (i-1)*W;
        nxm = j + (i-2)*W;
        nxp = j + i*W;
        nym = j-1+(i-1)*W;
        nyp = j+1+(i-1)*W;
        if i == 1
            G(n, n) = 1;
            % assume the current flows from left to right
            Z(n) = Vo;
        elseif i == L
            G(n, n) = 1;
            % by default Z(n)=0 here
        elseif j == 1  % lower bound
            if i > 0.4*L && i < 0.6*L % inside the blocks
                G(n, n) = -3;
                G(n, nyp) = sigma2;
                G(n, nxp) = sigma2;
                G(n, nxm) = sigma2;
            else
                G(n, n) = -3;
                G(n, nyp) = sigma1;
                G(n, nxp) = sigma1;
                G(n, nxm) = sigma1;
            end
        elseif j == W % upper bound
            if i > 0.4*L && i < 0.6*L % inside the block
                G(n, n) = -3;
                G(n, nym) = sigma2;
                G(n, nxp) = sigma2;
                G(n, nxm) = sigma2;
            else
                G(n, n) = -3;
                G(n, nym) = sigma1;
                G(n, nxp) = sigma1;
                G(n, nxm) = sigma1;
            end
        else 
            if i > 0.4*L && i < 0.6*L && (j < 0.4*W||j > 0.6*W)
                % inside the blocks
                G(n, n) = -4;
                G(n, nyp) = sigma2;
                G(n, nym) = sigma2;
                G(n, nxp) = sigma2;
                G(n, nxm) = sigma2;
            else
                G(n, n) = -4;
                G(n, nyp) = sigma1;
                G(n, nym) = sigma1;
                G(n, nxp) = sigma1;
                G(n, nxm) = sigma1;
            end
        end
    end
end
% G*V=Z
V3 = G\Z;
V4=reshape(V3,L,W);
[Ex, Ey] = gradient(V4);
Jx = S.*Ex;
Jy = S.*Ey;
J = sqrt(Jx.^2 + Jy.^2);
% assume R=1
V=zeros(L,W);
V(:,1)=1;
for x=1:L
    for y=2:W
    V(x,y)=V(x,y-1)-J(x,y)/L;
    end
end
[X,Y]=meshgrid(1:1:W,1:1:L);
figure(1)
surf(X,Y,V);
title('surface plot of V(x,y)');
[E1,E2]=gradient(V);
figure(2)
quiver(X,Y,E1,E2);
title('electric field  vector plot');
 