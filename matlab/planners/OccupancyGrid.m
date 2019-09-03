classdef OccupancyGrid < handle
    properties
        xDisc % Discretization in x (m / unit).
        yDisc % Discretization in y (m / unit).
        tDisc % Discretization in t (s / unit).
        
        xMin % Lower-bound on x (m).
        xMax % Upper-bound on x (m).
        yMin % Lower-bound on y (m).
        yMax % Upper-bound on y (m).
        tMin % Lower-bound on t (s).
        tMax % Upper-bound on t (s).
        
        data  % Cell array of grids (indexed by time).
        times % Array of times.
        
        figh % figure handle.
        hjiGrid;
        
        pathToFRSDir % Path to the directory containing the precomputed sets.
        
        xH;
    end
    
    methods
        function obj = OccupancyGrid(xDisc, yDisc, tDisc, xMin, xMax, yMin, ...
                yMax, tMin, tMax, hjiGrid, pathToFRSDir)
            obj.xDisc = xDisc;
            obj.yDisc = yDisc;
            obj.tDisc = tDisc;
            fprintf('Creating occupancy grid with disc (%f, %f, %f)\n', ...
                obj.xDisc, obj.yDisc, obj.tDisc);
            
            obj.xMin = xMin;
            obj.xMax = xMax;
            obj.yMin = yMin;
            obj.yMax = yMax;
            obj.tMin = tMin;
            obj.tMax = tMax;
            
            xExtent = round((xMax - xMin) / xDisc) + 1;
            yExtent = round((yMax - yMin) / yDisc) + 1;
            tExtent = round((tMax - tMin) / tDisc) + 1;
            
            obj.times = tMin:tDisc:tMax;
            obj.data = {};
            for t=1:length(obj.times)
                obj.data{t} = zeros(xExtent, yExtent);
            end
            
            obj.figh = [];
            if nargin > 9
                obj.hjiGrid = hjiGrid;
            else
                obj.hjiGrid = [];
            end
            
            if nargin > 10
                obj.pathToFRSDir = pathToFRSDir;
                fprintf('Path to FRS directory set to %s\n', ...
                    obj.pathToFRSDir);
            else
                fprintf('Warning: no path to FRS directory specified\n');
                obj.pathToFRSDir = '.';
            end
            
            obj.xH = [];
        end
        
        function fromValueFuns(obj, grid, valueFuns, times, t0)
            obj.times = zeros(1,length(times));
            obj.data = cell(1,length(times));
            obj.tMin = t0;
            obj.hjiGrid = grid;
            
            % We want to project all dimensions down into 2D.
            projDims = zeros(1,grid.dim);
            if grid.dim > 2
                projDims(3:end) = 1;
                
                for i=1:length(times)
                    % NOTE: this is hard-coded for 3D dyanamical systems.
                    vF = valueFuns(:,:,:,i);
                    
                    % Project to 2D by taking union over entire value function.
                    [~, data2D] = proj(grid, vF, projDims, 'min');
                    
                    % Convert to binary map.
                    data2D = 1*(data2D <= 0) + 0*(data2D > 0);
                    
                    % Store data and time-index based on starting time of
                    % prediction.
                    obj.data{i} = data2D;
                    obj.times(i) = t0 + times(i);
                end
            else
                for i=1:length(times)
                    data2D = valueFuns(:,:,i);
                    
                    % Convert to binary map.
                    data2D = 1*(data2D <= 0) + 0*(data2D > 0);
                    
                    obj.data{i} = data2D;
                    obj.times(i) = t0 + times(i);
                end
            end
        end
        
        function setRectangularObs(obj, lowXY, upXY, t)
            [lowi, lowj] = obj.xyToIndex(lowXY(1), lowXY(2));
            [upi, upj] = obj.xyToIndex(upXY(1), upXY(2));
            k = obj.timeToIndex(t);
            
            for i=lowi:upi
                for j=lowj:upj
                    obj.data{k}(i, j) = 1;
                end
            end
        end
        
        function setData(obj, x, y, t, val)
            [i, j] = obj.xyToIndex(x, y);
            k = obj.timeToIndex(t);
            
            obj.data{k}(i, j) = val;
        end
        
        function val = getData(obj, x, y, t, interpolateInTime)
            [i, j] = obj.xyToIndex(x, y);
            
            % If out-of-bounds, then assume val == 0.
            if i == 0 || j == 0 || ...
                    i > obj.hjiGrid.N(1) || j > obj.hjiGrid.N(2)
                fprintf('Coordinate (%f, %f, %f) is out of bounds!\n', ...
                    x, y, t);
                val = 0;
            else
                k = obj.timeToIndex(t);
                
                if nargin > 4 && interpolateInTime
                    % Interpolate in time.
                    [kBelow, kAbove] = obj.closestIndexToTime(t);
                    tBelow = obj.indexToTime(kBelow);
                    tAbove = obj.indexToTime(kAbove);
                    
                    alpha = (t - tBelow) / (tAbove - tBelow);
                    valBelow = obj.data{kBelow}(i, j);
                    valAbove = obj.data{kAbove}(i, j);
                    
                    val = valBelow + alpha * (valAbove - valBelow);
                else
                    if k > 0 && k <= length(obj.data)
                        val = obj.data{k}(i, j);
                    else
                        val = 0;
                    end
                end
            end
        end
        
        function [i, j] = xyToIndex(obj, x, y)
            if isempty(obj.hjiGrid)
                i = floor((x - obj.xMin) / obj.xDisc) + 1;
                j = floor((y - obj.yMin) / obj.yDisc) + 1;
            else
                % If the human state is set, then the center of the grid
                % should centered on this human state.
                xRel = x;
                yRel = y;
                
                if ~isempty(obj.xH)
                    xRel = x - obj.xH(1);
                    yRel = y - obj.xH(2);
                end
                
                % Check if the relative coordinate is out-of-bounds.
                if xRel < obj.hjiGrid.min(1) || xRel > obj.hjiGrid.max(1) || ...
                        yRel < obj.hjiGrid.min(2) || yRel > obj.hjiGrid.max(2)
                    i = 0;
                    j = 0;
                else
                    error = sqrt((obj.hjiGrid.xs{1} - xRel).^2 + ...
                        (obj.hjiGrid.xs{2} - yRel).^2);
                    [~,idx] = min(error(:));
                    [i, j] = ind2sub(size(error),idx);
                end
            end
        end
        
        function [x, y] = indexToXY(obj, i, j)
            x = obj.xMin + (i - 1) * obj.xDisc;
            y = obj.yMin + (j - 1) * obj.yDisc;
        end
        
        function k = timeToIndex(obj, t)
            k = floor((t - obj.tMin) / obj.tDisc) + 1;
        end
        
        function t = indexToTime(obj, k)
            t = obj.tMin + (k - 1) * obj.tDisc;
        end
        
        function [kBelow, kAbove] = closestIndexToTime(obj, t)
            idxLow = 1;
            idxHigh = length(obj.times);
            % Handle the case where t is out of bounds.
            if t <= obj.times(idxLow)
                kBelow = idxLow;
                kAbove = idxLow;
                return;
            elseif t >= obj.times(idxHigh)
                kBelow = idxHigh;
                kAbove = idxHigh;
                return;
            end
            
            kBelow = 0;
            kAbove = 0;
            
            idxMid = 0;
            
            % Binary search for the closest grids in time to the query time.
            while (idxHigh - idxLow) > 1
                idxMid = round(idxLow + (idxHigh - idxLow) / 2);
                tMid = obj.times(idxMid);
                fprintf("idxMid = %f, tMid = %f, t = %f\n", idxMid, tMid, t);
                if t < tMid
                    % Too high.
                    idxHigh = idxMid;
                elseif t > tMid
                    % Too low.
                    idxLow = idxMid;
                else
                    kBelow = idxMid;
                    kAbove = idxMid;
                    break;
                end
            end
            
            % Choose the indices of the grids before and after the query time.
            if kBelow < 1 && kAbove < 1
                if obj.times(idxMid) < t
                    kBelow = idxMid;
                    kAbove = min(idxMid + 1, length(obj.times));
                else
                    kBelow = max(1, idxMid - 1);
                    kAbove = idxMid;
                end
            end
        end
        
        function handles = draw(obj, color, prevHandles)
            numSlicesToSkip = 10;
            
            if nargin < 3
                for k = 1:numSlicesToSkip:length(obj.times)
                    handles{k} = scatter([], [], 30);
                end
            else
                handles = prevHandles;
            end
            
            linNums = linspace(0, 1, length(obj.times));
            
            xLB = obj.xMin;
            xUB = obj.xMax;
            yLB = obj.yMin;
            yUB = obj.yMax;
            
            if ~isempty(obj.xH)
                xLB = xLB + obj.xH(1);
                xUB = xUB + obj.xH(1);
                yLB = yLB + obj.xH(2);
                yUB = yUB + obj.xH(2);
            end
            
            for k = 1:numSlicesToSkip:length(obj.times)
                xs = [];
                ys = [];
                
                alpha = linNums(k) * 0.1 + (1 -linNums(k)) * 1;
                
                for x = xLB:obj.xDisc:xUB
                    for y = yLB:obj.yDisc:yUB
                        % Get the occupancy at this (x, y).
                        if obj.getData(x, y, obj.tMin + k * obj.tDisc) > 0
                            xs = [xs, x];
                            ys = [ys, y];
                        end
                    end
                end
                
                handles{k}.set('Xdata', xs);
                handles{k}.set('Ydata', ys);
                handles{k}.set('MarkerFaceColor', color);
                handles{k}.set('MarkerFaceAlpha', alpha);
                handles{k}.set('MarkerEdgeColor', 'none');
            end
            colormap('gray');
        end
        
        function setHumanState(obj, xH)
            obj.xH = xH;
        end
        
        function success = loadFromFile(obj, priorProb, horizon, currTime)
            success = true;
            
            % Load the data file.
            pathToFile = fullfile(obj.pathToFRSDir, ...
                strcat('pb', num2str(priorProb)), ...
                strcat('h', num2str(horizon)), ...
                'occuMap.mat');
            if isfile(pathToFile)
                fprintf('Loading data from %s\n', pathToFile);
                d = load(pathToFile);
            else
                fprintf('File %s does not exist!\n', pathToFile);
                d = {};
                success = false;
            end
            
            if success
                [~, ~, K] = size(d.predictions);
                obj.data = {};
                for k = 1:K
                    obj.data{k} = d.predictions(:, :, k);
                end
                
                obj.xDisc = (d.predGrid.max(1) - d.predGrid.min(1)) / ...
                    (d.predGrid.N(1) - 1);
                obj.yDisc = (d.predGrid.max(2) - d.predGrid.min(2)) / ...
                    (d.predGrid.N(2) - 1);
                
                obj.hjiGrid = d.predGrid;
                obj.times = currTime + d.predTimes;
                obj.tDisc = d.predDt;
                obj.tMin = currTime + d.predTmin;
                obj.tMax = currTime + d.predTmax;
            end
        end
    end
end