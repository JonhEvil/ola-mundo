%--------------------------------------------------------------------------
% Ola Simple FINITE ELEMENT CODE FOR a 2node (linear) 1D 2nd order diff. equation
%
% version 04 M.M. Neves - It is not a fully tested code. 
% This version * Mar 2020 (Octave, but runs in Matlab)
% Solve 
%       a*d2u/dx2+c*u=f     In P7 a=-1, c=1 and f=1  
%       u(0)=u(1)=0         Only Dirichelet (essential) boundary conditions
%
%
%------Prob 7 * This is the short script _AP3_FEA_1D_Prob7s.m 
clear all
clc
dt2pause=0
%
% Input data including Physical properties (pre-processing)
Nelem=4                 % Total Number of Finite Elements
LL=1                    % Total Lenght of the Domain (in m)
he=LL/Nelem             % Element lenght (m)
fixdofs=[1;Nelem+1]     % Nodal value of u fixed at first node
%
format rat
% Element stiffness matrix construction:
Ke=1/he*[1 -1; -1 1] + he*[1/3 1/6; 1/6 1/3]
% Distributed load vector construction (Only uniform load f=1 case):
Fe=he/2*[1 ; 1]
disp ("wait please...");
pause (dt2pause);
%
% Definition of global matrices:
K=sparse(Nelem+1,Nelem+1);
F=sparse(Nelem+1,1);
%
for ie=1:Nelem
edof=[ie,ie+1];
K(edof,edof)=K(edof,edof)+Ke;
F(edof,1)=F(edof,1)+Fe;
 full(K)
 full(F)
disp ("wait please...");
pause (dt2pause);
end
%
% Dirichelet Boundary conditions:
alldofs=1:Nelem+1
freedofs=setdiff(alldofs,fixdofs)
%
 full(K(freedofs,freedofs))
 full(F(freedofs,1))
%
% Solving Ku=F (processing)
U(freedofs,1)=K(freedofs,freedofs)\F(freedofs,1)
U(fixdofs,1)=0;
full(U)
%
format long
full(U)
disp ("wait please...");
pause (dt2pause);
%
% plots (postprocessing)
% 
clf
figure(1)
x=0:he/100:1; y=-0.2689414*exp(x)-0.7310585*exp(-x)+1; plot(x,y,'b-'); hold on;
x=0:he:1; plot(x,U,'r:o','LineWidth',2);
%x=0:he:1; plot(x,U,'g:o','LineWidth',2);
%x=0:he:1; plot(x,U,'k:o','LineWidth',2);
xlabel('x'); ylabel('u'); title('Solution of Prob. 7 with linear FE'); hold on;
legend('exact u(x)','FEM nodal w/ plot');
title('Solution of Prob. 7 with linear FE')

figure(2)
x=0:he/100:1; dy=-0.2689414*exp(x)+0.7310585*exp(-x); plot(x,dy,'b-'); hold on
for ie=1:Nelem
dfe=-(1/he)*U(ie,1)+(1/he)*U(ie+1,1);
xe=(ie-1)*he:he:ie*he; plot(xe,[dfe dfe],'m-o','LineWidth',2); hold on
%xe=(ie-1)*he:he:ie*he; plot(xe,[dfe dfe],'g-o','LineWidth',2); hold on
%xe=(ie-1)*he:he:ie*he; plot(xe,[dfe dfe],'k-o','LineWidth',2); hold on
end
xlabel('x'); ylabel('du'); 
legend('exact du(x)','FEM nodal du w/ plot')
title('Derivative of Solution of Prob. 7 with linear FE')




