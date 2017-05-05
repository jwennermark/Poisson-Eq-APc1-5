clear all;
clc;
close all;
%%
%Setting up the number of XY grid points
stepper=[1 2 4 5 10];
for i=1:length(stepper)
NXinc=100/stepper(i);
NYinc=NXinc;
Step=1/NXinc;

%%
%Setting up length of X and Y regions
Ax=-pi;
Bx=-Ax;
Ay=Ax;
By=Bx;

%Defining Li, X&Y domain lengths:
Lx=2*pi;
Ly=Lx;

%%
%Set up Min x&y coordinates
MinX1=Ay;
MinY1=MinX1;

%Set up Max x&y coordinates
MaxX1=By;
MaxY1=By;


%%
%Determining # elements  on X&Y domains and mesh
x=linspace(MinX1,MaxX1,NXinc);
y=linspace(MinY1,MaxY1,NYinc);
[xx,yy]=meshgrid(x,y);

Hx=x(2)-x(1);

%Dirichlet boundary condition inputs
psi=cos((pi.*(y-Ay))).*cosh(By-y); 
phi=((y-Ay).^2).*(sin((pi./2)*((x-Ax)/(By-Ay)).*2+1));
%%
%Defining my Boundary Conditions and Initial Conditions with pre populated
%zeros matrix to reduce run time
U=zeros(NYinc,NXinc);
U(NYinc,2:NXinc-1)=psi(2:NXinc-1);
U(1,2:NXinc-1)=phi(2:NXinc-1);

%%
%Related Forces initial matrix population to reduce run time
F=zeros(NYinc,NXinc);

%%
%Initializing error for Gauss-Seidel and bound for application of Nuemann
%boundary condition.
error=1;
bound=1;

%Setting Eset (acceptable limit of error at system convergence)
Eset=10^-3;

Lcount=0; %Initializing counter for # iterations during Gauss-Seidel

%Defining initial Tcheck value for checkpointing set time period of
%the Gauss-Seidel iteration loop:
Tinc = 60; %checkpoint every 60 seconds

%Defining initial timer variable value to be greater than start timer initial value to start timer loop in Gauss-Seidel loop:
Tint = 0;
Tcount = 0; %Stores time values for combined loop times, if loop takes less than a minute to run

save('Var.mat') %Retains variables if code fails

%%
%checkpointing, rerun
load('Var.mat')

%While error is less than Eset, save variables at Tinc
while error>Eset;
            tic;
        if (Tcount >= Tinc) %Saves at least every 60 seconds or every Gauss-Seidel loop iteration, if the iterations take longer
            Tcount = 0; %Resets time counter to 0
            save('Var') %Saves variables for checkpointing for stop/restart capability
        end

%Defining previous "U" for later use in SOR implementation        
Up=U;

%Force values defined by Dirichlet Boundary Conditions(constant) populating
%existing zeros matrix
for i=2:NXinc-1;
for j=2:NYinc-1;
% F(i,j)=sin(pi.*(x(i)-Ay)/(By-Ay)).*cos((pi/2).*(2.*((y(j) - Ay)./(By-Ay))+1));

%%
%Nuemann boundary conditions applied
if bound==1;
U(2:NYinc-1,NXinc)=(1/4)*(2*U(2:NYinc-1,NYinc-1)+U([2: NYinc-1]-1,NYinc)+U((2:NYinc-1)+1,NYinc)+( Hx^2)*F((2:NYinc-1),NYinc));
U(2:NYinc-1,NXinc)=(U(2:NYinc-1,NYinc-1)+U([2: NYinc-1]-1,NYinc)+U((2:NYinc-1)+1,NYinc)+(Hx^2)*F((2:NYinc-1),NYinc))/4;
bound=bound+1;
end;
U(i,j)=0.25*( U(i+1,j)+U(i-1,j)+ U(i,j+1)+ U(i,j-1) + F(i,j));
end
end

%Defining foward "U" value for implmentation of SOR & Lambda (lam)
Uf=U;

%%
%Defining Error post iterations
E=U-Up;
NX=NXinc-1;
error=mean(mean(E(2:NX,2:NX).^2));
Lcount=Lcount+1; %# of loop iterations

%SOR Implmentation (Lambda=lam and ranges 1:1.59
lam=1.5;
U=Uf*lam+(1-lam)*Up;

%Time elapsed retrieved
Tint=toc; %Determines how long the Gauss-Seidel iteration took
Tcount=Tcount+Tint; %time elapsed from checkpoint save
bound=1;
end

%%
%Plotting U
figure(i)
surf(xx,yy,U');
rotate3d
%Setting parameters of plot view, labels, and format
h.Label.String='U[Deviation]';
view(-37.5,30);
h=colorbar;
xlabel('x');ylabel('y');
axis tight
zlabel('U');
title('Gauss-Seidel Grid','fontweight','nor mal');
set(gca,'fontsize',18);
box on
end



