function [pf,ps,fc]= MOSOO(problem,l , u, numEvaluations, k)% temp func for vis
%% FUNCTION MODIRECT
%   A multi-objective optimizer from SOO
% problem : function handler problem : R^n -> R^k, gives a vector 1xk
% l : lower bound of the dimension space, 1xn
% u : upper bound of the dimension space, 1xn
% numEvaluations : evaluation budget
% k : number of functions
% p : value of h_max power

%-- Initialize the variables --------------------------------------%
% problem specific 
delta    = u - l;
n        = length(l);
%fCount   = 0;
func     = @(x) problem( delta.*x + l);
%keyboard
% tree structure
level = zeros(numEvaluations,1); % level of the node, determines its expansion and selection
c = zeros(numEvaluations,n); % centre point of the node, representative state
fc = zeros(numEvaluations, k); % f-value at that centre 1*k
%dimSeq = bsxfun(@times, 1 :  n,ones(numEvaluations, n)); % the sequence of the dimension splits which the node should follow in expansion
% other variables
I = eye(n);
pf = []; ps = [];
idx = [];
% variables for governing the tree depth
%H_MAX =  log(numEvaluations^3); fixed
H_MAX = @(x,y) ceil( x + log(2*(numEvaluations - y))/log(3)+n^(1.5)); % it was 10*x^p
%hCount = 0;
% -- Initialize --
front = 1;
c(1,:) = 0.5 * ones(1,n);
fc(1,:) = func(c(1,:));
pf(1,:) = fc(1,:);
fCount = 1;
iter = 0;
%prevfCount = 1; % to keep record of how many (original )
level(1) = 0;
candidateNodes = 1;
%dimSeq(1,:) = 1:n;
verbose = false;
if verbose && n ==2
    figure(2)
    rectangle('Position', [ 0 0 1 1])
    hold on
%     % level set:
%     [x,y]= ndgrid(-1:0.01:1);
%     x1 = [x(:) y(:)]';
%     %y = linspace(0,4*pi);
%     %[X,Y] = meshgrid(x,y);
%     Z = visFunc(x1);
% 
%     figure(2)
%     contour(x,y,reshape(Z(1,:),[201 201]))
%     hold on
%     figure(2)
%     contour(x,y,reshape(Z(2,:),[201 201]))
end
%attribute = 'EE'; % EE: explore-exploit , EED : explore-exploit-discover
splitLevel = 1;
% -- Main Loop ---
while (fCount < numEvaluations)
	%----------------------------------------------------------------
	% 1. expand the potentially optimal nodes, one by one:
	%----------------------------------------------------------------
	iter = iter + 1;
	for fId = 1 : length(front)
		i = front(fId);
		% sample 
		offset = (1/3)^(1 + floor(level(i)/n)); % compute the offset of the samples from the centre point
		% book-keeping variables
		splitDim = mod(level(i),n) + 1;
		numDim = 1;% (n - firstDim + 1);
		numSamples = 2 ;%* numDim;
		firstIdx = fCount + 1;
		lastIdx = fCount + numSamples;
		% check if you need more space for the coming samples:
		if (lastIdx > numEvaluations)
			level(firstIdx : lastIdx) = zeros(numSamples,1);
			c(firstIdx : lastIdx,:) = zeros(numSamples,n);
			fc(firstIdx : lastIdx,:) = zeros(numSamples, k);
		end
		c(firstIdx : lastIdx, :) = bsxfun(@plus, c(i, :), [offset .* I(splitDim,:); -offset .* I(splitDim,:)]);
		% evaluate samples
		for j = 1 : numDim
			fc(fCount + j,:) = func(c(fCount + j,:));
			fc(fCount + numDim + j, :)= func(c(fCount + numDim + j,:));
			level(fCount + j,:) = level(i) + 1;
			level(fCount + numDim + j, :)= level(i) + 1;
		end
		% update the level 
		level(i) = level(i) + numDim;
		% visualize
		if (verbose && (n == 2 || n ==3))
			figure(2)
            xlim([0 1]);
            ylim([0 1])
            set(gca, 'XTick', []);
            set(gca, 'YTick', []);
            if firstIdx==2
                firstIdx = 1;
            end
			if (n==2) % 2 dimension
				scatter(c(firstIdx : lastIdx,1),c(firstIdx : lastIdx,2),'.k')
                
                for j = 1 : numDim
					rectangle('Position',[c(fCount + j,1) - 3/2.*(1/3).^(1 + floor(level(fCount + j)./n) + mod(level(fCount + j),2)*(splitDim==1)),...
										c(fCount + j,2) - 3/2.*(1/3).^(1 + floor(level(fCount + j)./n) + mod(level(fCount + j),2)*(splitDim==2)),... 
										3.*(1/3).^(1 + floor(level(fCount + j)./n) + mod(level(fCount + j),2)*(splitDim==1)),...
										3.*(1/3).^(1 + floor(level(fCount + j)./n) + mod(level(fCount + j),2)*(splitDim==2))]);
					rectangle('Position',[c(fCount + numDim + j,1) - 3/2.*(1/3).^(1 + floor(level(fCount + numDim + j)./n) + mod(level(fCount + numDim + j),2)*(splitDim==1)),...
										c(fCount + numDim + j,2) - 3/2.*(1/3).^(1 + floor(level(fCount + numDim + j)./n) + mod(level(fCount + numDim + j),2)*(splitDim==2)),... 
										3.*(1/3).^(1 + floor(level(fCount + numDim + j)./n) + mod(level(fCount + numDim + j),2)*(splitDim==1)),...
										3.*(1/3).^(1 + floor(level(fCount + numDim + j)./n) + mod(level(fCount + numDim + j),2)*(splitDim==2))]);

				end
			else % 3 dimension
				scatter3(c(firstIdx : lastIdx,1),c(firstIdx : lastIdx,2),c(firstIdx : lastIdx,3),'.k')
			end
			hold on
            %fc(firstIdx : lastIdx,:)
            %c(firstIdx:lastIdx,:)
            %iter
			%keyboard
            pause(0.1)
        end
       	
		% update fCount
		fCount = fCount + numSamples;
		if (fCount >= numEvaluations)
			break;
		end
	end
	%----------------------------------------------------------------
	% 2. identify the potentially optimal nodes at the split level:
	%----------------------------------------------------------------
	% take the candidates fromp the previous round and the candidate at the current split level
	splitLevel = splitLevel + 1;
	minLevel = min(level(1:fCount));
	if (splitLevel > H_MAX(minLevel, fCount) || splitLevel > max(level(1:fCount)))
	  splitLevel = minLevel;
	  % check if the level is beyond
	  if (splitLevel > H_MAX(minLevel, fCount))
		disp('Tree exhaused !')
		%fCount
		%break;
	  end
	  % compare only with candidates from the current sweep of the tree
	  prevCandidates= [];
	  %prevfCount = fCount;
	  %disp('restarting');
	else
	  prevCandidates = candidateNodes;
	end
	%candidateNodes = find(level(1:prevfCount) == splitLevel);
	candidateNodes = find(level(1:fCount) == splitLevel);
	[~,ia]= intersect(prevCandidates,candidateNodes);
	prevCandidates(ia)= []; % to ensure only childs are considered.
	% get the selected node IDs
	front = (paretofront(fc([prevCandidates; candidateNodes],:)));
	front = candidateNodes(front(numel(prevCandidates)+1:end));
	%----------------------------------------------------------------
	% 3. Get the approximation set (no need at the time being)
	%frontSet = paretofront([pf; fc(fCount-numSamples:fCount,:)]);
	%pf = fc(frontSet,:);
	%----------------------------------------------------------------
end
% -- End of Main Loop --

% Get the approximation set
front = paretofront(fc);
pf = fc(front,:);
ps = bsxfun(@plus, bsxfun(@times,delta, c(front,:)),l);
end
