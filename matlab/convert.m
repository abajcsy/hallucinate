clear all
clc

originalBaseDir = '/home/eratner/Documents/hallucinate/bayes_two_goal';
newBaseDir = '/home/eratner/Documents/hallucinate/bayes_two_goal_new/';

parfor idx = 0:10
    b = idx * 0.1;
    bdir = sprintf('pb%0.2f', b);
    
    if ~isdir(fullfile(newBaseDir, bdir))
        mkdir(fullfile(newBaseDir, bdir));
    end
    
    for xval = -3:0.5:3
        for yval = -1.5:0.5:1.5
            name = sprintf('xval_%0.2f_yval_%0.2f_occuMap.mat', xval, yval);
            pathToNew = fullfile(newBaseDir, bdir, name);
            pathToOriginal = fullfile(originalBaseDir, bdir, name);
            convertData(pathToOriginal, pathToNew);
        end
    end
end

% convertData('/home/eratner/Documents/hallucinate/bayes_two_goal/pb0.00/xval_-1.00_yval_0.50_occuMap.mat', ...
%             '/home/eratner/Documents/hallucinate/bayes_two_goal_new/convered.mat');
% convertData('/home/eratner/Documents/hallucinate/bayes_two_goal/pb0.50/xval_-1.00_yval_0.50_occuMap.mat', ...
%             '/home/eratner/Documents/hallucinate/bayes_two_goal_new/convered2.mat');

function convertData(pathToOriginalFile, pathToNewFile)
    predTmin = 0;
    predTmax = 6.2;
    predDt = 0.05;

%     humanSpeed = 0.6;
%     triSideLength = 0.1;
%     oldDt = triSideLength / humanSpeed;

    oldData = load(pathToOriginalFile);
    
    [oldRows, oldCols] = size(oldData.predGridX);

    predTimes = linspace(predTmin, predTmax, ...
        1 + (predTmax - predTmin) / predDt);

    gridLow = [-4; -4; predTmin];
    gridUpp = [4; 4; predTmax];
    N = [81; 81; 125];
    predGrid = createGrid(gridLow, gridUpp, N);
    predictions = zeros(N(1), N(2), N(3));
    
    fprintf('Converting file %s...', pathToOriginalFile);

    tStart = cputime;
    
    oldDataPredTimes = oldData.predDt * linspace(0, 34, 35);
    
    for tIdx = 1:N(3)
        for xIdx = 1:N(1)
            for yIdx = 1:N(2)
                x = predGrid.xs{1}(xIdx, yIdx, tIdx);
                y = predGrid.xs{2}(xIdx, yIdx, tIdx);
                t = predGrid.xs{3}(xIdx, yIdx, tIdx);
                
                error = abs(oldDataPredTimes - t);
                [~, bestTimeIdx] = min(error(:));
                
                % Get data point in old data that is closest to this (x, y)
                error = sqrt((oldData.predGridX - x).^2 + ...
                    (oldData.predGridY - y).^2);
                [minErr, idx] = min(error(:));
                if minErr < 0.1
                    [bestRow, bestCol] = ind2sub(size(error), idx);
                
                    bestGridIdx = (bestRow - 1) * oldCols + bestCol;
                
                    predictions(xIdx, yIdx, tIdx) = oldData.predictions(bestTimeIdx, bestGridIdx);
                end
            end
        end
    end
    
    tEnd = cputime;
    fprintf('done. Conversion took %f s\n', tEnd - tStart);
    
    fprintf('Saving data to %s...', pathToNewFile);
    save(pathToNewFile, 'predictions', 'predGrid', 'predTmin', 'predTmax', 'predDt', 'predTimes');
    fprintf('done\n');
end