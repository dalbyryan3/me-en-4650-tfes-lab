%% ME EN 4650  Lab6:CFD    Ryan Dalby    
% 1f
close all;

% CD
openfig('NACA0012_CD.fig');
hold on;
plot(5, 0.01785, 'rx', 'DisplayName', 'Simulated at alpha = 5° and Re_c = 1.5e05', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
plot(12, 0.05359, 'bx', 'DisplayName', 'Simulated at alpha = 12° and Re_c = 1.5e05', 'MarkerSize', 10, 'LineWidth', 2);
title('Experimental and Simulated Coefficent of Drag Values');


openfig('NACA0012_CL.fig');
hold on;
plot(5, 0.4961, 'rx', 'DisplayName', 'Simulated at alpha = 5° and Re_c = 1.5e05', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
plot(12, 0.3097, 'bx', 'DisplayName', 'Simulated at alpha = 12° and Re_c = 1.5e05', 'MarkerSize', 10, 'LineWidth', 2);
title('Experimental and Simulated Coefficent of Lift Values');