clc
clear all
close all

ind=1;
NumberOfClasses = 2;
Atoms = 64;
len=32;


A = (randn(len,Atoms));
for i=1:Atoms
    A(:,i) = A(:,i)/norm(A(:,i),2);
end
Dictionary = A;
AtA = A'*A;

TrueSparsityLevel = 10;
Noise_Std = 0.01;

x0 = (randn(Atoms,1));
x0( randperm(Atoms, Atoms - TrueSparsityLevel) ) = 0;
x0 = x0/max(abs(x0))/2;
% x0=x0/norm(x0);
TrueMajorElements = find(x0~=0);
TotalAtoms_True = length(TrueMajorElements);

n0 =   randn(len,1)*Noise_Std;
y0= [];
y0 = [y0  (A*x0+n0) / norm((A*x0+n0),2) ];
y = y0;

Lambda =0.0002; 
Sigma = 0.018;
Kappa = 0.47 * ones(Atoms,1);
Rho = Sigma^2 * log ( ((2*pi*Sigma^2)/Lambda)  * ((1-Kappa )./Kappa ).^2 );

%%
[x,Avg_x,Time_ICR,Objective] = ICR_Func(y,A,AtA,'lambda',Lambda,'rho',Rho,'GroundTruth',x0, ...
    'algorithm', 1, 'verbose', 1);

x_j = x;
Final_gamma =  abs(x_j)./abs(Avg_x);
MajorElements = find( abs(Final_gamma - 1 ) <0.5 );
Objective 
Time_ICR
