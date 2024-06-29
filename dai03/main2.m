clear
clc

%% 初始化
T=150;
A=[1 0.05;-0.0005 1];
B=[0;0.05];
%D=[0.3;0.2];
%C=[0.2 0.3];
E=0.1;
H1 = 0.3;
H2 = 1.1;
H = H2-H1;
U = 0.5;
V = 0.4;
M = 0.2;
S = 0.5; 
lamda = 0.5;
Wm = 0.1; 
epsilon_M= 0.8 ; 
x=zeros(2,T);
x(:,1)=[0.1; 0.1];
%u=-0.5*x(1,1)-2.5*x(2,1);
x_hat=zeros(2,T);
x_hat(:,1)=[0.1; -0.2];
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
    C=[0.5+0.01*cos(k) 0.6+0.01*sin(k)];
    D=[0.04*sin(k)+0.2;0.3+0.02*cos(k)];
    u=-0.5*x(1,k)-2.5*x(2,k);
    x(:,k+1)=A*x(:,k)+ [0;-0.0335*x(1,k)^3]+B*u+D*0.3*cos(k);
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
    innovation = y(k+1)-C*x_hat(:,k+1); %计算新息
    phi=Sat(innovation,sigma(k));
    
    %solve beta1k, beta2k
    [ beta1, beta2, k1, k2, P1_bar ] = lmi1( F1, C, phi, Phi, Wm );
   
    Wf_hat=beta1*Wf_hat+beta2*C*Sat((y(k+1)-C*x_hat(:,k+1)),sigma(k+1)); %更新权重

    
end

%% plot
figure(1);
t=1:1:T;
plot(t,x(1,:),'b--^',t,x_hat(1,:),'r--*');
ylabel('x1')
xlabel('k')
%ylim([-0.4,0.4])
figure(2);
t=1:1:T;
plot(t,x(2,:),'b--^',t,x_hat(2,:),'r--*');
ylabel('x2')
xlabel('k')
%ylim([-0.4,0.4])
figure(3);
t=1:1:T;
plot(t,x(1,:)-x_hat(1,:),'b--^');
ylabel('e')
xlabel('k')
ylim([-2,2])
figure(4);
t=1:1:T;
plot(t,x(2,:)-x_hat(2,:),'b--^');
ylabel('e')
xlabel('k')
ylim([-2,2])