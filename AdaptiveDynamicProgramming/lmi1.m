function [ beta1, beta2, k1, k2, P1_bar ] = lmi1( F1,C, phi, Phi, Wm )
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
    
    Wm2=Wm^2;
    a=kron(C'*phi, Phi);

    setlmis([]);%初始化lmi系统 注意只描述上三角就可以了
    %% 设定矩阵变量
    k1 = lmivar(1,[1,0]);
    k2 = lmivar(1,[1,0]);
    beta1 = lmivar(1,[1,0]);
    beta2 = lmivar(1,[1,0]);
    P1_bar= lmivar(1,[4,1]);
    
    %% 描述线性矩阵不等式
    %(1,1)
    lmiterm([1,1,1,k1],1,1);%k1
    lmiterm([1,1,1,k1],1,Wm2);%k2*Wm^2
    lmiterm([1,1,1,0],-1);%-1
    %(2,2)
    lmiterm([1,2,2,k2],-1,1)%-k2*I
    %(3,3)
    lmiterm([1,3,3,k1],-1,1)%-k1*I
    %(4,1)
    lmiterm([1,4,1,beta2],-1,a);
    %(4,2)
    lmiterm([1,4,2,beta1],-1,1);%-beta1*I
    lmiterm([1,4,2,0],1);%I
    %(4,3)
    lmiterm([1,4,3,beta1],1,F1);
    %(4,4)
    lmiterm([1,4,4,P1_bar],-1,1);
    
    %% 正定约束
    %k1,k2,beta1,beta2,P1_bar>0
    lmiterm([-2 1 1 k1], 1, 1);
    lmiterm([-3 1 1 k2], 1, 1);
    lmiterm([-4 1 1 beta1], 1, 1);
    lmiterm([-5 1 1 beta2], 1, 1);
    lmiterm([-6 1 1 P1_bar], 1, 1);

    lmisys = getlmis;%获取lmi的信息

    %% 求可行解
    [tmin, xfeas] = feasp(lmisys);
    
    beta1 = dec2mat(lmisys, xfeas, beta1);
    beta2 = dec2mat(lmisys, xfeas, beta2);
    k1 = dec2mat(lmisys, xfeas, k1);
    k2 = dec2mat(lmisys, xfeas, k2);
    P1_bar = dec2mat(lmisys, xfeas, P1_bar);

end

