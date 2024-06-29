clear
clc

%% 初始化
T=500;
A=[1 0.05;-0.0005 1];
B=[0;0.05];
%D=[0.3;0.2];
%C=[0.2 0.3];
Q=eye(2);
R=0.5;
alpha_bar = 0.6;
alpha=binornd(1,alpha_bar,T,1);
rho_w = 0.01 ;
bata_a = 0.5;
bata_c = 0.9;
bata_d = 0.5;
V = zeros(1,T);
x=zeros(2,T);
x(:,1)=[0.1;0.1];
u = zeros(1,T);
u=-0.5*x(1,1)-2.5*x(2,1);
mu = zeros(1,T);
W_do = zeros(6,T);
W_do(:,1)=[0.0,0.1,0.2,0.3,0.4,0.5];
W_di = ones(6,2);
W_ao = zeros(6,T);
W_ao(:,1)=[0.0,0.1,0.2,0.3,0.4,0.5];
W_ai = ones(6,2);
W_co = zeros(8,T);
W_ci = ones(8,2);
e_d = zeros(1,T);
e_c = zeros(1,T);
e_a = zeros(1,T);

%% sys
for k=1:T-1
    %C=[0.2+0.01*cos(k) 0.3+0.02*sin(k)];
    D=[0.04*sin(k)+0.2;0.3+0.02*cos(k)];
    x(:,k+1)=A*x(:,k)+ [0;-0.0335*x(1,k)^3]+B*(u(k)+alpha(k)*mu(k))+D*0.3*cos(1.5*k);
    
    %V更新
    deltaphi =  tanh(W_ci*x(:,k+1)) - tanh(W_ci*x(:,k));
    r =x(:,k)'*Q*x(:,k) + u(k)'*R*u(k) - alpha_bar*mu(k)'*mu(k);
    e_c(k) = W_co(:,k)'*deltaphi + r;
    W_co(:,k+1) = W_co(:,k) - bata_c*deltaphi*e_c(k); %梯度下降
    V(k+1) = W_co(:,k+1)'*tanh(W_ci*x(:,k+1));

    %控制更新
    u_t = -0.5*R^(-1)*B'*W_ci'*calculate_tanh_jacobian(tanh(W_ci*x(:,k+1)))*W_co(:,k+1);
    e_d(k) = W_do(:,k)'*tanh(W_di*x(:,k)) - u_t;
    W_do(:,k+1) = W_do(:,k) - bata_d*tanh(W_di*x(:,k))*e_d(k);
    u(k+1) = W_do(:,k+1)' * tanh(W_di*x(:,k+1));

    %攻击更新
    mu_t = 0.5*B'*W_ci'*calculate_tanh_jacobian(tanh(W_ci*x(:,k+1)))*W_co(:,k+1);
    e_a(k) = W_ao(:,k)'*tanh(W_ai*x(:,k)) - mu_t;
    W_ao(:,k+1) = W_ao(:,k) - bata_a*tanh(W_ai*x(:,k))*e_a(k);
    mu(k+1) = W_ao(:,k+1)' *  tanh(W_ai*x(:,k+1));

    
end

figure(1);
t=1:1:T;
plot(t,x(1,:),'b');
ylabel('x1')
xlabel('k')
ylim([-5,5])
figure(2);
t=1:1:T;
plot(t,x(2,:),'b');
ylabel('x2')
xlabel('k')
%ylim([-0.1,0.1])
figure(3);
t=1:1:T;
plot(t,W_do(1,:),'b',t,W_do(2,:),'b');
ylabel('x2')
xlabel('k')
%ylim([-0.1,0.1])
