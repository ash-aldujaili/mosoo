


% General settings of the benchmark:
addpath(genpath(pwd))
resultDir= 'numerical_data';
benchmark = 'theory';
algorithms = {'MO-SOO'}; 
dimensions = [1,2]; % or can be used as a replcaement for the number of instances
functions = 1;
numRun = 4;
% run the benchmark for each of the algorithms above:
 for algIdx = 1 : length(algorithms)
    evaluateAlg;
 end
% copy the theoretical files to 
% compile and generate the tables:
system('python compileTheoreticalPlots.py');



