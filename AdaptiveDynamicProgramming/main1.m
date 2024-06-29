clear
clc

%% 初始化
T=50;
A=[0.1 0.3;0.3 0.1];
B=[0;0.5];
%D=[0.3;0.2];
%C=[0.2 0.3];
E=0.1;
u=0.1;
H1 = 0.9;
H2 = 1;
H = H2-H1;
S = 0.01; 
U = 0.5;
V = 0.3;
M = 0.3;
lamda = 0.99;
Wm = 1; 
epsilon_M= 1.7 ; 
x=zeros(2,T);
x(:,1)=[0.2; -0.1];
x_hat=zeros(2,T);
x_hat(:,1)=[0.1;-0.2];
Wf_hat=eye(2);
P1 = 100*eye(2);
P2 = 100*eye(2);
P1_bar = kron(P1,eye(2));
P1_bar_inv = inv(P1_bar);
y=zeros(1,T);
sigma=zeros(1,T);
sigma(1)=0.3;

%% sys
for k=1:T-1
    C=[0.2+0.01*cos(k) 0.3+0.02*sin(k)];
    Ck=[0.2+0.01*cos(k+1) 0.1+0.02*sin(k+1)];
    D=[0.04*sin(k)+0.2;0.3+0.02*cos(k)];
    x(:,k+1)=A*x(:,k)+ f(x(:,k))+B*u+D*0.3*cos(1.5*k);
    y(k)=C*x(:,k)+E*0.3*cos(k)+0.1*sin(k);
    y(k+1)=C*x(:,k+1)+E*0.3*cos(k+1)+0.1*sin(k+1);
    sigma(k+1)=lamda*sigma(k)+S*(y(k)-C*x_hat(:,k))^2;
    
    Phi=tanh(x_hat(:,k));
    PHI=[Phi', 0, 0; 0, 0, Phi'];
    F1=chol(P1_bar)';
    F2=chol(P2)';
    
     %solve L
    [ rho1, rho2, rho3, rho4, rho5, rho6, rho7, rho8, L, P2] = lmi2(k,epsilon_M,sigma,lamda,PHI,F1,F2,A,C,D,S,U,V,E,M,H,H1);
    
    x_hat(:,k+1)=A*x_hat(:,k)+Wf_hat'*tanh(x_hat(:,k))+B*u+L*Sat((y(k)-C*x_hat(:,k)),sigma(k));
    innovation = y(k+1)-Ck*x_hat(:,k+1); %计算新息
    phi=Sat(innovation,sigma(k));
    
    %solve beta1k, beta2k
    [ beta1, beta2, k1, k2, P1_bar ] = lmi1( F1, C, phi, Phi, Wm );
   
    Wf_hat=beta1*Wf_hat+beta2*Ck'*Sat((y(k+1)-Ck*x_hat(:,k+1)),sigma(k+1))*Phi'; %更新权重

    
end

%% plot
figure(1);
t=1:1:T;
plot(t,x(1,:),'b--hexagram',t,x_hat(1,:),'r--*');
xlabel('k')
i=legend('$x_1$','$\hat{x}_1$');
set(i,'Interpreter','latex');
%ylim([-0.5,0.5]);
figure(2);
t=1:1:T;
plot(t,x(2,:),'b--hexagram',t,x_hat(2,:),'r--*');
i=legend('$x_2$','$\hat{x}_2$');
set(i,'Interpreter','latex');
xlabel('k')
%ylim([-0.5,0.5])
figure(3);
t=1:1:T;
plot(t,x(1,:)-x_hat(1,:),'b--hexagram',t,x(2,:)-x_hat(2,:),'r--*');
xlabel('k')
%ylim([-2,2])
legend('e_1','e_2')