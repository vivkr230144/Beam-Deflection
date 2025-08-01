clear all;
close all;
L=10; %beam lenght, in m
b=1;  %beam width, in m
h=0.323; %slab thickness, in m 
E=3.2837*10^7; %Young's modulus of concrete, in kPa (kN/m^2)
q=8.5; %UDL, in kN/m
I=b*h^3/12; %moment of inertia
syms x
y=q*x^4/(24*E*I)-5*q*L*x^3/(48*E*I)+q*L^2*x^2/(16*E*I);
k=64*10^7; %soil springs stiffness, in kPa (kN/m^2)
n=200; %number of steps
Dx=L/n; %step-size
RHS=zeros(n+1, 1);
RHS(3:n-1,1)=Dx^4*q/(E*I);
C=zeros(n+1, n+1);
C(1,1)=1;
C(n+1,n+1)=1;
C(2,1)=-3;
C(2,2)=4;
C(2,3)=-1;
C(n, n-2)=-1;
C(n, n-1)=4;
C(n, n)=-5;
C(n, n+1)=2;
for i=3:n-1
C(i,i-2)=1;
C(i,i-1)=-4;
C(i,i)=6+k*Dx^4/(E*I);
C(i,i+1)=-4;
C(i,i+2)=1;
end
y_num=linsolve(C,RHS);
x_num=(0:Dx:10)';
fplot(x,-y, [0,10]);
title ("x (m) vs -y(m)")
xlabel("x (m) ")
ylabel ("-y(m)")
hold on
scatter(x_num, -y_num); %in red 
p=k*y_num;
for i=1:n
    pi=Dx*(p(i,1)+p(i+1,1))/2;
    mi=pi*(i-0.5)*Dx;
    sum_pi=sum(pi,"all");
    sum_mi=sum(mi,"all");
end
x_bar=sum_mi/sum_pi;


