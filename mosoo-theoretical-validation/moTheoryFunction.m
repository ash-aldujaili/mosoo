function y = moTheoryFunction(x,isMinus, isThining,samplingStep, numIteration,instance)
% motheoryFunction represent a wrapper and a logger for multi-objective
% theoretical function used to demonstrate MO-SOO loss bounds
% isMinus : true if the algorithm is designed for maximzing the function
% isIterationBased : true if you would like to record the samples based on
% the thining strategy, attaching the time stamp of the sample (reduce the
% size of the data)
% numIteration : the current timestamp of the sample
% samplingStep : recod it based on the value of the counter relative to the
% Sampling step, set it to one to record all the samplings
% while we are interested in minimizng it.
global counter;
%global bestX;
%global bestY;
%global logX;
global logY;
global logIdx;



counter = counter + 1;


% make sure x is a column vector
if (size(x,1) < size(x,2))
    x = x';
end

%y = feval(fhd,x,varargin{:});
alpha_1 = 1.2;
alpha_2 = 0.9;
opt1_vals = [0 0.21 0.47 0.57];
opt2_vals = [1 0.81 0.61 0.57];

opt_1 = opt1_vals(instance);
opt_2 = opt2_vals(instance);
y = [(max(abs(x-opt_1))^alpha_1), (max(abs(x-opt_2))^alpha_2)]; % objectives 

if (counter == 1)
    %bestY = y;
    %bestX = x';
    logY(1,:)= [1,y];
    logIdx =logIdx + 1;
    %logX(1,:)= bestX;
else

    %if (y < bestY)
    %    bestY = y;
    %bestX = x';
    %end

    if (mod(counter,samplingStep)==0) % 50 interval
        % perform thining (i.e if the new sample is non-dominated)
        if (isThining)
            front = paretofront([logY(:,2:end);y]);
            if (front(end))
                logIdx = logIdx + 1;
                logY(logIdx,:)= [numIteration, y];
            end
        else
            logIdx = logIdx + 1;
            logY(logIdx,:)= [numIteration, y];
            %numIteration
        end
    %logX(logIdx,:)= bestX;
    end
end

if isMinus % maximize/minimizaion problem
    y = -y; % minimize problem
end


end


