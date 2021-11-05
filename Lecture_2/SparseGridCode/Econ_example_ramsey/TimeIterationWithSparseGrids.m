%**************************************************************************
% Solving 2-Dimensional Ramsey Model Using Time-Iteration with Sparse Grids 
%**************************************************************************
%
% This script solves a 2-dimensional Ramsey model with a time iteration 
% algorithm that interpolates the policy function using sparse grids.
%
%**************************************************************************
%
% Uses the spinter Sparse Grid Interpolation Toolbox from Andreas Klimke:
% http://www.ians.uni-stuttgart.de/spinterp/
%
%**************************************************************************
%
% Calls FOCs.m
% 
%**************************************************************************
% By Johannes Brumm, Simon Scheidegger 07/2018
%**************************************************************************

%% Initialization:
clear all;

% Load Sparse Grid Library
addpath('../spinterp_v5.1.1')  

tic % start clock

% Economic paramters to be chosen: 
beta  = 0.95; % discount factor
gamma = 2; % risk aversion / inverse of intertemporal elasticity
delta = 0.1;  % rate of depreciation
alpha = 0.3;  % capital share
sigma = 0.02; % std of productivity shocks
pers=0.7; % persistence of productivity process

% We use a Gauss-Hermite quadrature rule (see Judd 1998, p.262):
wquad=pi^(-1/2)*[0.2954 1.1816 0.2954]; % quadrature weights
pquad=sqrt(2)*sigma*[-1.225 0.0 1.225]; % (scaled) quadrature points

% Computational parameters to be chosen:
tol_EE = 10^-8; % error tolerance for solution to the Euler equation
tol_it = 10^-6; % error tolerance of time-iteration algorithm

% set options for solver:
solver_op=optimset('Display','off','TolFun',tol_EE,'MaxIter',500,'TolX',10^-10);
% set options for SPARSE GRID:

% +++ SPARSE GRID HERE, choose from different options +++

% Test Global Polynomial Basis Functions: (smaller RelTol -> finer grid)
 grid_op =  spset('GridType', 'Chebyshev','DimensionAdaptive', 'off','RelTol',10^-1);

% Test Piecewise Linear Basis Functions:
% grid_op =  spset('GridType', 'Clenshaw-Curtis','RelTol',10^-1); warning off;
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Parameters to be determined:
kstar=((1/beta-(1-delta))/alpha)^(1/(alpha-1)); % steady state capital stock
k_min = kstar*0.3;   % minimal capital stock
k_max = kstar*1.5;   % maximum capital stock
a_min = 1 - 5*sigma; % minimal productivity level
a_max = 1 + 5*sigma; % maximum productivity level
error_it = 1; % initial value for error in timeiteration
test_it=[rand(100,1)*(k_max-k_min)+k_min rand(100,1)*(a_max-a_min)+a_min];   % define a set of test points for iteration
test_ee=[rand(2000,1)*(k_max-k_min)+k_min rand(2000,1)*(a_max-a_min)+a_min]; % define a set of test points for Euler errors

%initial guess for time iteration:
p = @(k,a) 0.95*k+0.05*(k_max+k_min)/2;
z = spvals(p,2,[k_min k_max; a_min a_max],grid_op);

%% Time Iteration:

% Run time-iteration algorithm until error below tolerance
tic
while error_it > tol_it
    
    % Use FOCs and policy function guess to get new policy function:
    pnew = @(k,a) fsolve( @(v) FOCs(v,k,a,alpha,beta,delta,gamma,pers,wquad,pquad,z),k,solver_op);
    
    % Get sparse grid representation of policy function:
    znew = spvals(pnew,2,[k_min k_max; a_min a_max],grid_op);
    % Calculate iteration error measure and update policy guess:
    error_it=max(abs(spinterp(znew,test_it)-spinterp(z,test_it)))/(max(spinterp(z,test_it)))
    z=znew; p=pnew;
    
    % plot policy along iteration if necessary:
    % ezmesh(@(k,a) spinterp(z,k,a),[k_min, k_max, a_min, a_max]); pause(1)
end
toc

%% Error Analysis:

k=test_ee(:,1); a=test_ee(:,2);
kp=spinterp(z,[k a]);
for i=1:size(test_ee(:,1),1)
    EulerErrors(i)=abs(FOCs(kp(i),k(i),a(i),alpha,beta,delta,gamma,pers,wquad,pquad,z));
end

% Compute errors
MaxEulerError=max(EulerErrors); % calculate the maximum Euler error 
AverageEulerError=mean(EulerErrors); % calculate the average Euler error 
display(['Method: Time Iteration With Sparse Grid of type ',num2str(z.gridType),' with ',num2str(z.nPoints),' Points.']) 
display(['Maximum EE: ',num2str(MaxEulerError),'. Average EE: ',num2str(AverageEulerError),'.'])
%**************************************************************************

