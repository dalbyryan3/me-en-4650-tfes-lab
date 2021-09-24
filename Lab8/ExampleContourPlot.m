% ExampleContourPlot.m
%
% This code demonstrates how to generate filled contour plots and line
% contour plots for a 2D function.
%
% M. Metzger
% 2-2018

%% User-specified parameters
Lx = 2*pi;     %length of domain in x-direction (cm)
Ly = 1;        %length of domain in y-direction (cm)
delta = 0.005; %spacing between points (cm)
Ncontour = 15; %number of colors/lines to use in contour plot

%% Create the 2D function to display

% Create vectors containing discrete x & y values in the domain
x=[0:delta:Lx];
y=[0:delta:Ly];

% create 2D grid of coordinates in domain: Xmat and Ymat are matrices
[Xmat,Ymat]=meshgrid(x,y);

% create an interesting 2D function representing height of a surface
Zmat = 0.5*Ymat.^2+0.1*cos(Xmat);

%% Generate a Color Contour Plot

% generate a filled contour plot with 'Ncontour' different colors
figure;
contourf(Xmat,Ymat,Zmat,Ncontour);
xlabel('X (cm)');
ylabel('Y (cm)');

% add a colorbar to the figure with appropriate label
cb=colorbar;
cb.Label.String='Z (cm)';

%% Generate a Line Contour Plot

% determine the values of the iso-lines (these will be lines of constant Z)
maxZ=max(max(Zmat));
Z_isolines=[0:maxZ/Ncontour:maxZ];

% generate a contour plot with iso-lines
figure;
[C,hnd]=contour(Xmat,Ymat,Zmat,Z_isolines);
xlabel('X (cm)');
ylabel('Y (cm)');

% set the linestyle of the isolines and turn on labels
set(hnd,'color','k','linestyle','-');
clabel(C,hnd);

