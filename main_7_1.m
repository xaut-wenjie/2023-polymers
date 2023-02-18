function main_7_1
%---this code was edited by wenjie sun on Jan 6,2023
%---偏置电压的影响
format long
clc,clear;
global L_0 B A lambda_a_p lambda_p_p fre mu J_lim H_0 phi_ac xita M epsilon phi_dc
%---结构与材料参数
A=15e-3;        %A表示圆形质量块的半径；
B=35e-3;        %B表示环形约束框架的内径；
L_0=28e-3;      %L_0表示环形约束框架的间距；
lambda_a_p=1.5;
lambda_p_p=1.5;
M=15.16e-3; %M为中间质量块的总质量；
%---材料参数来自于文献(SMS-2019,AIS-2020)
H_0=224e-6;
epsilon=2.74*8.85e-12;
mu=620e3;
J_lim=538;


%---主动区域施加的电压
for i=1:121
    if i<=31
        fre=20+(i-1)*0.1;
    end
    if (i>31 && i<=61)
        fre=23.4+(i-31)*0.01;
    end
    if (i>61)
        fre=24+(i-61)*0.1;
    end
    phi_ac=2000;
    phi_dc=3000;
    xita=0.122932264219044*fre.^2-5.783882047130196*fre+68.172890656419824;

    x_0=0;
    x_end=60;
    h=0.001;
    N_step=(x_end-x_0)/h;
    
    %计算初值
    phi=phi_dc;
    fun=@(D)f(D,A,B,L_0,H_0,lambda_a_p,lambda_p_p,phi,epsilon,mu,J_lim);
    initial=L_0/2;
    y_initial=fzero(fun,initial);

    y_0=[y_initial;0];
    RR=rk4(@odefun,x_0,x_end,y_0,N_step);
    time=RR(:,1);
    dis=(RR(:,2)-y_initial)*1000;
    NN=length(dis);
    max_1=max(dis(NN-200:NN));
    min_1=min(dis(NN-200:NN));
    Amp(i)=(max_1-min_1)/2;
    FF(i)=fre;
    i=i+1;
end
plot(FF,Amp);
data_theory=[FF' Amp'];
csvwrite('phi_dc_3000_sun.csv',data_theory);








function dy=odefun(x,y)

%F_p*sin(theta_p)-F_a*sin(theta_a)-xita*(dD/dt)=M*(ddD/dt^2)
%D=y(1);dD/dt=y(2)

%F_p=2*pi*A*h_p*sigma_p_r;
%F_a=2*pi*A*h_a*sigma_a_r;
%sigma_p_r=mu*(lambda_p_r^2-lambda_p_r^-2*lambda_p_p^-2)/...
%             (1-(lambda_p_r^2+lambda_p_p^2+lambda_p_r^-2*lambda_p_p^-2-3)/J_lim);
%sigma_a_r=mu*(lambda_a_r^2-lambda_a_r^-2*lambda_a_p^-2)/...
%             (1-(lambda_a_r^2+lambda_a_p^2+lambda_a_r^-2*lambda_a_p^-2-3)/J_lim)-...
%             (epsilon*phi^2/H_0^2)*lambda_a_r^2*lambda_a_p^2;
%lambda_p_r=sqrt((L_0-y(1))^2+(B-A)^2)/(B-A)*lambda_p_p;
%lambda_a_r=sqrt(y(1)^2+(B-A)^2)/(B-A)*lambda_a_p;

global L_0 B A lambda_a_p lambda_p_p fre mu J_lim H_0 phi_ac xita M epsilon phi_dc
dy=zeros(2,1);

%---sin(theta_p)、sin(theta_a)
sin_theta_p=(L_0-y(1))/sqrt((L_0-y(1))^2+(B-A)^2);
sin_theta_a=y(1)/sqrt(y(1)^2+(B-A)^2);
%---lambda_p_r、lambda_a_r
lambda_p_r=sqrt((L_0-y(1))^2+(B-A)^2)/(B-A)*lambda_p_p;
lambda_a_r=sqrt(y(1)^2+(B-A)^2)/(B-A)*lambda_a_p;
%---h_p、h_a
h_p=H_0*lambda_p_r^-1*lambda_p_p^-1;
h_a=H_0*lambda_a_r^-1*lambda_a_p^-1;
%---phi
phi=phi_dc+phi_ac*sin(2*pi*fre*x);
%---sigma_p_r、sigma_a_r
sigma_p_r=mu*(lambda_p_r^2-lambda_p_r^-2*lambda_p_p^-2)/...
             (1-(lambda_p_r^2+lambda_p_p^2+lambda_p_r^-2*lambda_p_p^-2-3)/J_lim);
sigma_a_r=mu*(lambda_a_r^2-lambda_a_r^-2*lambda_a_p^-2)/...
             (1-(lambda_a_r^2+lambda_a_p^2+lambda_a_r^-2*lambda_a_p^-2-3)/J_lim)-...
             (epsilon*phi^2/H_0^2)*lambda_a_r^2*lambda_a_p^2;
%---F_p、F_a
F_p=2*pi*A*h_p*sigma_p_r;
F_a=2*pi*A*h_a*sigma_a_r;

%动力学方程；
dy(1)=y(2);
dy(2)=-xita*y(2)/M+F_p*sin_theta_p/M-F_a*sin_theta_a/M;
         
         
         

%------------四阶龙格库塔法函数
function R=rk4(odefun,x_0,x_end,y_0,N_step)
h=(x_end-x_0)/N_step;
x=x_0:h:x_end;
N=length(x);
x0=x_0;
y0=y_0;
for j=1:N
    K1=h*odefun(x0,y0);
    K2=h*odefun(x0+h/2,y0+K1/2);
    K3=h*odefun(x0+h/2,y0+K2/2);
    K4=h*odefun(x0+h,y0+K3);
    y1=y0+(K1+2*K2+2*K3+K4)/6;
    tt(j)=x0;
    yy1(j)=y1(1);
    yy2(j)=y1(2);
    x0=x0+h;
    y0=y1;
    j=j+1;
end
R=[tt' yy1' yy2'];





function y=f(D,A,B,L_0,H_0,lambda_a_p,lambda_p_p,phi,epsilon,mu,J_lim)

%---sin(theta_p)、sin(theta_a)
sin_theta_p=(L_0-D)/sqrt((L_0-D)^2+(B-A)^2);
sin_theta_a=D/sqrt(D^2+(B-A)^2);
%---lambda_p_r、lambda_a_r
lambda_p_r=sqrt((L_0-D)^2+(B-A)^2)/(B-A)*lambda_p_p;
lambda_a_r=sqrt(D^2+(B-A)^2)/(B-A)*lambda_a_p;
%---h_p、h_a
h_p=H_0*lambda_p_r^-1*lambda_p_p^-1;
h_a=H_0*lambda_a_r^-1*lambda_a_p^-1;
%---sigma_p_r、sigma_a_r
sigma_p_r=mu*(lambda_p_r^2-lambda_p_r^-2*lambda_p_p^-2)/...
             (1-(lambda_p_r^2+lambda_p_p^2+lambda_p_r^-2*lambda_p_p^-2-3)/J_lim);
sigma_a_r=mu*(lambda_a_r^2-lambda_a_r^-2*lambda_a_p^-2)/...
             (1-(lambda_a_r^2+lambda_a_p^2+lambda_a_r^-2*lambda_a_p^-2-3)/J_lim)-...
             (epsilon*phi^2/H_0^2)*lambda_a_r^2*lambda_a_p^2;
%---F_p、F_a
F_p=2*pi*A*h_p*sigma_p_r;
F_a=2*pi*A*h_a*sigma_a_r;

%---
y=F_p*sin_theta_p-F_a*sin_theta_a;

    























