function [ rho1, rho2, rho3, rho4, rho5, rho6, rho7, rho8, L, P2] = lmi2( k,epsilon_M,sigma,lamda,PHI,F1,F2,A,C,D,S,U,V,E,M,H,H1)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
    %% 初始化
    epsilon_M2=epsilon_M^2;
    lamdasigma = lamda*sigma(k);
    sigma_k1 = sigma(k+1);
    U_inv = inv(U);
    V_inv = inv(V);
    M_inv = inv(M);
    fcscf = F2'*C'*S*C*F2;
    escf = E'*S*C*F2;
    ese = E'*S*E;
    scf = S*C*F2;
    se = S*E;
    hcf = H*C*F2;
    h1cf= H1*C*F2;
    he = H * E;
    phif1=PHI*F1;
    
    setlmis([]);%初始化lmi系统 注意只描述上三角就可以了
    %% 设定矩阵变量
    rho1 = lmivar(1,[1,0]);
    rho2 = lmivar(1,[1,0]);
    rho3 = lmivar(1,[1,0]);
    rho4 = lmivar(1,[1,0]);
    rho5 = lmivar(1,[1,0]);
    rho6 = lmivar(1,[1,0]);
    rho7 = lmivar(1,[1,0]);
    rho8 = lmivar(1,[1,0]);
    L = lmivar(2,[2,1]);  %2*1阶矩阵
    P2 = lmivar(1,[2,1]); %2阶对称矩阵
    
    %% 描述线性矩阵不等式
    %(1,1)
    lmiterm([1,1,1,0],-1);%-1
    lmiterm([1,1,1,rho1],1,1);
    lmiterm([1,1,1,rho2],1,1);
    lmiterm([1,1,1,rho3],1,1);
    lmiterm([1,1,1,rho4],1,1);
    lmiterm([1,1,1,rho5],1,1);
    lmiterm([1,1,1,rho6],1,epsilon_M2);
    lmiterm([1,1,1,rho7],1,sigma_k1);
    lmiterm([1,1,1,rho7],-1,lamdasigma);
    
    %(2,2)
    lmiterm([1,2,2,rho1],-1,1);
    
    %(3,3)
    lmiterm([1,3,3,rho2],-1,1);
    lmiterm([1,3,3,rho7],-1,fcscf);
    
    %(4,4)
    lmiterm([1,4,4,rho3],-1,U_inv);
    
    %(5,3)
    lmiterm([1,5,3,rho7],-1,escf);
    
    %(5,5)
    lmiterm([1,5,5,rho4],-1,V_inv);
    lmiterm([1,5,5,rho7],-1,ese);
    
    %(6,3)
    lmiterm([1,6,3,rho7],-1,scf);
    
    %(6,5)
    lmiterm([1,6,5,rho7],-1,se);
    
    %(6,6)
    lmiterm([1,6,6,rho5],-1,M_inv);
    lmiterm([1,6,6,rho7],-1,S);
    
    %(7,3)
    lmiterm([1,7,3,rho8],0.5,hcf);
    
    %(7,5)
    lmiterm([1,7,5,rho8],0.5,he);
    
    %(7,6)
    lmiterm([1,7,6,rho8],0.5,H);
    
    %(7,7)
    lmiterm([1,7,7,rho8],-1,1);
    
    %(8,8)
    lmiterm([1,8,8,rho6],-1,1);
    
    %(9,2)
    lmiterm([1,9,2,0],phif1);
    
    %(9,3)
    lmiterm([1,9,3,0],A*F2);
    lmiterm([1,9,3,L],-1, h1cf);
    
    %(9,4)
    lmiterm([1,9,4,0],D);
    
    %(9,5)
    lmiterm([1,9,5,L],-1, H1*E);
    
    %(9,6)
    lmiterm([1,9,6,L],-1, H1);
    
    %(9,7)
    lmiterm([1,9,7,L],-1, 1);
    
    %(9,8)
    lmiterm([1,9,8,0],1);
    
    %(9,9)
    lmiterm([1,9,9,P2],-1, 1);
 
    
    
    %% 正定约束
    %rho1, rho2, rho3, rho4, rho5, rho6, rho7, rho8,P2>0
    lmiterm([-2 1 1 rho1], 1, 1);
    lmiterm([-3 1 1 rho2], 1, 1);
    lmiterm([-4 1 1 rho3], 1, 1);
    lmiterm([-5 1 1 rho4], 1, 1);
    lmiterm([-6 1 1 rho5], 1, 1);
    lmiterm([-7 1 1 rho6], 1, 1);
    lmiterm([-8 1 1 rho7], 1, 1);
    lmiterm([-9 1 1 rho8], 1, 1);
    lmiterm([-10 1 1 P2], 1, 1);

    

    %% 求可行解
    lmisys = getlmis;%获取lmi的信息
    
    [tmin, xfeas] = feasp(lmisys);
    
    rho1 = dec2mat(lmisys, xfeas, rho1);
    rho2 = dec2mat(lmisys, xfeas, rho2);
    rho3 = dec2mat(lmisys, xfeas, rho3);
    rho4 = dec2mat(lmisys, xfeas, rho4);
    rho5 = dec2mat(lmisys, xfeas, rho5);
    rho6 = dec2mat(lmisys, xfeas, rho6);
    rho7 = dec2mat(lmisys, xfeas, rho7);
    rho8 = dec2mat(lmisys, xfeas, rho8);
    L = dec2mat(lmisys, xfeas, L);
    P2 = dec2mat(lmisys, xfeas, P2);

end

