% file names
algorithm = algorithms{algIdx};
prefixFile = 'resultfile';
timeFile = fullfile(resultDir, 'timing.txt');
prefix = fullfile(resultDir, prefixFile);
% general settings
numDim = length(dimensions);
numFunc = length(functions);
timeInfo = zeros(numFunc,numDim);

% check if the folder exists
if (isdir(fullfile(resultDir)))
  delete(fullfile(resultDir,'resultfile*.txt')); % delete a prcd file from the folder of results
else
  mkdir(fullfile(resultDir))
end


if strcmp(benchmark,'theory')
   maxRange = 1;
   minRange = -1;
   samplingStep = 1;
   for dimIdx = 1: numDim    
        dimension = dimensions(dimIdx);
        numEvaluations = 10000 ;
        resLog = zeros(numFunc,numRun);  
        for func_num = 1:numFunc
            for runtime = 1: numRun
                % logging:
                global counter;
                global logX;
                global logY;
                global bestX;
                global bestY;
                global logIdx;
                logIdx = 0;
                bestY = inf;
                bestX = [];
                counter = 0;
                logX = [];
                logY = [];
                % MO-SOO
                if (strcmp(algorithm,'MO-SOO'))
                    n = dimension; % decision space dimension
                    m = 2; % objective space dimensio
                    p = 0.5; % parameter of h_max
                    l = zeros(n,1)';% lower bound of the decision space
                    u = ones(n,1)';% upper bound of the decision space
                    isMaximization = false;
                    isParetoRecord = true;
                    func = @(x,iter) moTheoryFunction(x, false,isParetoRecord, samplingStep, iter, runtime);
                    % execute the algorithm
                    [DEPTH{dimIdx,runtime}, fCount{dimIdx,runtime},paretofront, ~,~]=MOSOO(func,l , u, numEvaluations, m, p);
                else
                    error ('No such algorithm')
                end
                % print the log
                % Reporting results==========================================================
                % write down the results of the current evaluation:
                name = [ prefix '-f' num2str(func_num) '-d' num2str(dimension) '-r' num2str(runtime) ];
                count= samplingStep:samplingStep:numEvaluations;
                % extend when one algorithm terminates before the number of
                % evaluations
                lenX = length(logY);
                fidFx = fopen([name '-fx.txt'], 'w');
                fprintf(fidFx,'%%count\t x\n');
                A = logY;
                dlmwrite([name '-fx.txt'], A,'delimiter', '\t', '-append');
                fclose(fidFx);    
            end
        end
        % Write stats
        fid= fopen([prefix '-stat-d' num2str(dimension) '.txt'], 'w');
        fprintf(fid,'%%		Statistics result of \n%% Func.	Best	Worst	Median	Mean	Std\n');
        A = zeros(numFunc, 6);
        for i = 1 : numFunc
            A(i,:) = [ i min(resLog(i,:)) max(resLog(i,:)) median(resLog(i,:)) mean(resLog(i,:)) std(resLog(i,:))];
        end
        dlmwrite([prefix '-stat-d' num2str(dimension) '.txt'], A,'delimiter', '\t', '-append');
   end 
   fclose(fid);
else
  error('no such benchmark !');
end

% write depth info
for n = 1 : 2
    for i = 1 : 4
       dlmwrite(sprintf('hmax_data/hmax-n%d-mosoo-i%d.txt',n,i), [[1:numel(DEPTH{n,i})]' DEPTH{n,i}], 'delimiter', '\t') 
    end
end

% write timing info:
timeInfo = timeInfo./numRun;
timeInfo =[ [1:numFunc]' timeInfo];
dlmwrite(timeFile, timeInfo);
