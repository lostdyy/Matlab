clear
clc

%% 初始化
T=500;
A=[0.1 0.3;0.3 0.1];
B=[0;0.5];
%D=[0.3;0.2];
%C=[0.2 0.3];
E=0.1;
Q=eye(2);
R=0.5;
alpha_bar = 0.3;
alpha=binornd(1,alpha_bar,T,1);
H1 = 0.3;
H2 = 1.1;
H = H2-H1;
U = 0.5;
V = 0.3;
M = 0.3;
S = 0.01;
lamda = 0.98;
Wm = 1; 
epsilon_M= 1.7 ; 
x=zeros(2,T);
x(:,1)=[0.6;-0.6];
x_hat=zeros(2,T);
x_hat(:,1)=[0.1;-0.2];
Wf_hat=[0.1,0;0.2,-0.2];%eye(2);
Wf =zeros(4,T);
P1 = 100*eye(2);
P2 = 100*eye(2);
P1_bar = kron(P1,eye(2));
P1_bar_inv = inv(P1_bar);
y=zeros(1,T);
sigma=zeros(1,T);
sigma(1)=0.3;
bata_a = 0.7;
bata_c = 0.4;
bata_d = 0.5;
Va = zeros(1,T);
u = zeros(1,T);
u(1) = 0.1;
mu = zeros(1,T);
W_do = zeros(6,T);
W_do(:,1)=[0.0,0.1,0.2,0.3,0.4,0.5];
W_di = ones(6,2);
W_ao = zeros(6,T);
W_ao(:,1)=[0.0,0.1,0.2,0.3,0.4,0.5];
W_ai = ones(6,2);
W_co = zeros(8,T);
W_co(:,1)=[-0.1,-0.2, 0.0 ,0.1,0.2,0.3,0.4,0.5];
W_ci = ones(8,2);
e_d = zeros(1,T);
e_c = zeros(1,T);
e_a = zeros(1,T);

%% sys
for k=1:T-1
    Ck=[0.2+0.01*cos(k) 0.3+0.02*sin(k)];
    Ckp1=[0.2+0.01*cos(k+1) 0.3+0.02*sin(k+1)];
    Dk=[0.04*sin(k)+0.2; 0.3+0.02*cos(k)];
    x(:,k+1)=A*x(:,k)+ f(x(:,k))+B*(u(k)+alpha(k)*mu(k))+Dk*0.01*cos(1.5*k);
    y(k)=Ck*x(:,k)+E*0.3*cos(k);
    y(k+1)=Ck*x(:,k+1)+E*0.3*cos(k+1);
    sigma(k+1)=lamda*sigma(k)+S*(y(k)-Ck*x_hat(:,k))^2;
    
    Phi=tanh(x_hat(:,k));
    PHI=[Phi', 0, 0; 0, 0, Phi'];
    F1=chol(P1_bar)';
    F2=chol(P2)';
     
     %solve L
    [ rho1, rho2, rho3, rho4, rho5, rho6, rho7, rho8, L, P2] = lmi2(k,epsilon_M,sigma,lamda,PHI,F1,F2,A,Ck,Dk,S,U,V,E,M,H,H1);
    
    %估计器
    x_hat(:,k+1)=A*x_hat(:,k)+Wf_hat'*tanh(x_hat(:,k))+B*(u(k)+alpha(k)*mu(k))+L*Sat((y(k)-Ck*x_hat(:,k)),sigma(k));
    innovation = y(k+1)-Ck*x_hat(:,k+1); %计算新息
    phi=Sat(innovation,sigma(k));
    
    %solve beta1k, beta2k
    [ beta1, beta2, k1, k2, P1_bar ] = lmi1( F1, Ckp1, phi, Phi, Wm );
   
    %更新权重
    Wf_hat=beta1*Wf_hat+beta2*Ckp1'*Sat((y(k+1)-Ckp1*x_hat(:,k+1)),sigma(k+1))*Phi'; 
    Wf(:,k+1) =  reshape(Wf_hat, 4, 1);
    
     %V更新
    deltaphi =  tanh(W_ci*x_hat(:,k+1)) - tanh(W_ci*x_hat(:,k));
    r =x_hat(:,k)'*Q*x_hat(:,k) + u(k)'*R*u(k) - alpha_bar*mu(k)'*mu(k);
    e_c(k) = W_co(:,k)'*deltaphi + r;
    W_co(:,k+1) = W_co(:,k) - bata_c*deltaphi*e_c(k); %梯度下降
    Va(k+1) = W_co(:,k+1)'*tanh(W_ci*x_hat(:,k+1));

    %控制更新
    u_t = -0.5*R^(-1)*B'*W_ci'*calculate_tanh_jacobian(tanh(W_ci*x_hat(:,k+1)))*W_co(:,k+1);
    e_d(k) = W_do(:,k)'*tanh(W_di*x_hat(:,k)) - u_t;
    W_do(:,k+1) = W_do(:,k) - bata_d*tanh(W_di*x_hat(:,k))*e_d(k);
    u(k+1) = W_do(:,k+1)' * tanh(W_di*x_hat(:,k+1));

    %攻击更新
    mu_t = 0.5*B'*W_ci'*calculate_tanh_jacobian(tanh(W_ci*x_hat(:,k+1)))*W_co(:,k+1);
    e_a(k) = W_ao(:,k)'*tanh(W_ai*x_hat(:,k)) - mu_t;
    W_ao(:,k+1) = W_ao(:,k) - bata_a*tanh(W_ai*x_hat(:,k))*e_a(k);
    mu(k+1) = W_ao(:,k+1)' *  tanh(W_ai*x_hat(:,k+1));
   
end

%% plot
figure(1);
t=1:1:T;
plot(t,x(1,:),'b--hexagram',t,x_hat(1,:),'r--*');
xlabel('k')
i=legend('$x_1$','$\hat{x}_1$');
set(i,'Interpreter','latex');
%ylim([-0.7,0.3]);
figure(2);
t=1:1:T;
plot(t,x(2,:),'b--hexagram',t,x_hat(2,:),'r--*');
i=legend('$x_2$','$\hat{x}_2$');
set(i,'Interpreter','latex');
xlabel('k')
%ylim([-0.7,0.3])
figure(3);
t=1:1:T;
plot(t,x(1,:)-x_hat(1,:),'b--hexagram',t,x(2,:)-x_hat(2,:),'r--*');
xlabel('k')
%ylim([-2,2])
legend('e_1','e_2')
figure(4);
t=1:1:T;
plot(t,u,'b',t,mu,'r');
xlabel('k')
%ylim([-0.1,0.1])
legend('u','\mu')
figure(5);
t=1:1:T;
plot(t,W_do(1,:),'b',t,W_do(2,:),'r',t,W_do(3,:),'g',t,W_do(4,:),'c',t,W_do(5,:),'m',t,W_do(6,:),'k');
xlabel('k')
%ylim([-0.1,0.1])
figure(6);
t=1:1:T;
plot(t,W_ao(1,:),'b',t,W_ao(2,:),'r',t,W_ao(3,:),'g',t,W_ao(4,:),'c',t,W_ao(5,:),'m',t,W_ao(6,:),'k');
xlabel('k')
%ylim([-0.1,0.1])
figure(7);
t=1:1:T;
plot(t,W_co(1,:),'b',t,W_co(2,:),'r',t,W_co(3,:),'g',t,W_co(4,:),'c',t,W_co(5,:),'m',t,W_co(6,:),'k',t,W_co(7,:),'y',t,W_co(8,:),'g');
xlabel('k')
%ylim([-0.1,0.1])
figure(8);
t=1:1:T;
plot(t,Wf(1,:),'g',t,Wf(2,:),'c',t,Wf(3,:),'b',t,Wf(4,:),'r');
xlabel('k')
%ylim([-0.1,0.1])