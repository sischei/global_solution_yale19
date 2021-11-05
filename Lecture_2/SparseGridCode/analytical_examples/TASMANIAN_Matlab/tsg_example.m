function tsg_example()

close all;
clear all;
clc;
addpath('../../TasmanianSparseGrids/InterfaceMATLAB');

plot_choice = 0; % 0=make no plots; 1=make plots


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% tsg_example()
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXAMPLE 1:
%
% interpolate: f(x,y) = exp( -x^2 ) * cos( y )   -- The so-called "Darth Vader" Function
% using a "Classical Sparse Grid"
%
% We compute the max error for 1000 random points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic Parameters for the sparse grid

dim =   2;        %dimension of problem
outs =  1;        %how many outputs
l_min = 5;        %refinement level of Sparse grid
which_basis = 1;  % 1 = linear basis functions 

disp(['----------------------------------------------------------------------------']);
disp(['    Example 1: interpolate: f(x,y) = exp( -x^2 ) * cos( y ) ']);
disp(['    using "classical" sparse grid with depth ',num2str(l_min)]);
disp(['    the error is estimated as the maximum from 1000 random points']);
disp([' ']); % clear all temporary used files


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build-up sparse grid

%%  "name of sparse grid"
grid_name = 'example_grid';   

%% Classical Sparse grid of max level l_min
[ lGrid, points ] = tsgMakeLocalPolynomial( grid_name, dim, outs, 'localp', l_min, which_basis );

%% analytical test function
vals = ( exp( -points(:,1).^2 ) .* cos( points(:,2) ) );  
tsgLoadValues( lGrid, vals );                                

%% Generate 1000 random points in (x,y)
pnts = [ -1 + 2 * rand( 1000, 2 ) ];              

%% Evaluate analytical function at the 1000 random points 
tres = exp( -pnts(:,1).^2 ) .* cos( pnts(:,2) );  

%% Evaluate interpolant at the same 1000 random points
[ res ] = tsgEvaluate( lGrid, pnts );             

%% number of points in Sparse Grid
disp(['   Number of points: ',num2str( size( points, 1 ) )]);                  
%% compute error
disp(['   Max Error for 1000 points : ',num2str( max( abs( res - tres ) ) )]); 


%% Example on how to evaluate sparse grid at a single point 

%% Generate 1 random point in (x,y)
pnts2 = [ -1 + 2 * rand( 1, 2 ) ];              

%% Evaluate analytical function at the 1 random point 
one_res_analyt = exp( -pnts2(:,1).^2 ) .* cos( pnts2(:,2) ); 
one_res_analyt;

%% Evaluate interpolant at the same 1 random point
tic
[ one_res ] = tsgEvaluate( lGrid, pnts2);            
toc
disp([' ']);
%% compute error
disp(['   Error for 1 point: ',num2str( max( abs( one_res_analyt - one_res) ) )]);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Results:

if plot_choice==1

%%% plot the sparse grid points
figure;
plot(points(:,1),points(:,2), 'Marker','x','LineStyle','none');
grid on;
box on;

%%%%%%%%%%%
%plot the sparse grid points and the corresponding function value
figure;
[ res ] = tsgEvaluate( lGrid, points );             % Evaluate interpolant at the grid points
plot3(points(:,1),points(:,2), res,'Marker','+','LineStyle','none');
grid on;
box on;
hold on;
title({'Darth Vader'});

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% here its how you could store the Sparse Grid  -> see documentation p.50ff
%filename = 'Sparse_grid.txt';  %filename
%save(filename, 'lGrid');       %store file there

tsgDeleteGrid( lGrid ); %you have to delete the grid










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXAMPLE 2:
%
% interpolate: f(x,y) = exp( -x^2 ) * cos( y )
% we refine the sparse grid adaptively
%
% We comute the max error for 1000 random points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic Parameters for the sparse grid

dim =  2;         % number of dimensions
outs = 1;         % how many outputs
l_min = 2;        % linear basis function
tol = 1.E-5;      % refinement criterion
which_basis = 1;  % 1 = linear basis functions 

disp(['----------------------------------------------------------------------------']);
disp(['    Example 2: interpolate: f(x,y) = exp( -x^2 ) * cos( y )  ']);
disp(['    the error is estimated as the maximum from 1000 random points']);
disp(['    tolerance is set at 1.E-5']);
disp([' ']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build-up sparse grid

%% "name of sparse grid"
grid_name_2 = 'example_grid_2';  

%% sparse grid of level l_min
[ lGrid1, points ] = tsgMakeLocalPolynomial(grid_name_2, dim, outs, 'localp', l_min, which_basis ); 

%% analytical test function
vals = ( exp( -points(:,1).^2 ) .* cos( points(:,2) ) ); 
tsgLoadValues( lGrid1, vals );      

%% generate 1000 random points in (x,y)
pnts = [ -1 + 2 * rand( 1000, 2 ) ]; 

%% evaluate analytical test function at (x,y)
tres = ( exp( -pnts(:,1).^2 ) .* cos( pnts(:,2) ) );

disp(['             Ordinary refinement']);
disp([' iteration   nodes       error  ']);

nump1 = size( points, 1 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build-up ADAPTIVE sparse grid

%% here we refine the grid up to level 7
for iI = 1:7
    tt = num2str(iI);
    ss = [blanks(6 - length(tt)),tt];
    %% ordinary refinement criterion    
    [ points ] = tsgRefineSurplus( lGrid1, tol, 'classic' );  
    %% analytical test function
    vals = ( exp( -points(:,1).^2 ) .* cos( points(:,2) ) );  
    tsgLoadValues( lGrid1, vals );
    nump1 = nump1 + size( points, 1 );
    tt = num2str(nump1);
    ss = [ss,'      ',blanks(5 - length(tt)),tt];
    [ res ] = tsgEvaluate( lGrid1, pnts );
    tt = num2str(max(abs( res - tres )),5);
    ss = [ss,' ',blanks(12 - length(tt)),tt];
    ss = [ss,blanks( 30 - length(ss) )];
    disp(ss);
end
disp([' ']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Results:

if plot_choice==1

%%%%%%%%%%%
% Plot adaptive sparse grid
figure;
grid on;
plot(points(:,1),points(:,2), 'Marker','x','LineStyle','none');
grid on;
box on;
ylim([-1 1]);
xlim([-1 1]);

%%%%%%%%%%%
% plot the sparse grid points and the corresponding function value
figure;
[ res ] = tsgEvaluate( lGrid1, points );             % Evaluate interpolant at the grid points
plot3(points(:,1),points(:,2), res,'Marker','+','LineStyle','none');
grid on;
box on;
hold on;
ylim([-1 1]);
xlim([-1 1]);

title({'Darth Vader with an adaptive sparse grid '});

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tsgDeleteGrid( lGrid1 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



