% This is a demo of the multi-objective optimistic optimization algorithm
% (MO-SOO)
% The open-source code and data will be made available at the projects website.
% Author : Abdullah Al-Dujaili 2015

clear;clc;
close all;
% define a problem

n = 2; % decision space dimension
m = 2; % objective space dimensio
p = 0.5; % parameter of h_max
l = -ones(n,1)';% lower bound of the decision space
u = ones(n,1)';% upper bound of the decision space
f = @(x) [ ((x(1)-0.25).^2+(x(2)-0.66).^2) ; ((x(1)+0.25).^2+(x(2)-0.66).^2)]';
numEvaluations = 250;
% execute the algorithm
[PF,PS,fc]= MOSOO(f,l , u, numEvaluations, m, p);






% plot approximation set and pareto front:
% get a sampled set of the Pareto front
[x,y,z]= ndgrid(-1:0.01:1);
x = [x(:) y(:) z(:)]';
f = @(x) [ ((x(1,:)-0.25).^2+(x(2,:)-0.66).^2) ; ((x(1,:)+0.25).^2+(x(2,:)-0.66).^2)]';
y = f(x);

front = paretofront(y);
pf = y(front,:);
figure(1)
scatter(pf(:,1),pf(:,2),'b.')
figure(1), hold on; scatter(PF(:,1),PF(:,2),'r*')
legend('Pareto front','Approximation set')
