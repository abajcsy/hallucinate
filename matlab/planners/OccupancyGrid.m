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
           if i == 0 || j == 0
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
                   gridOrigin = [obj.xH(1) - (obj.hjiGrid.max(1) - obj.hjiGrid.min(1)) / 2;
                                 obj.xH(2) - (obj.hjiGrid.max(2) - obj.hjiGrid.min(2)) / 2];
                   xRel = x - gridOrigin(1);
                   yRel = y - gridOrigin(2);
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
       
        function draw(obj, grid2D, color)
            linNums = linspace(0, 1, length(obj.times));
            if ~isempty(obj.figh)
                for i=1:length(obj.figh)
                    delete(obj.figh{i});
                end
            end
            
            for i=1:10:length(obj.times)
                vF = obj.data{i};
                posIdx = find(vF >= 1);
                alpha = linNums(i)*0.1 + (1-linNums(i))*1;
                obj.figh{end+1} = scatter(grid2D.xs{1}(posIdx), grid2D.xs{2}(posIdx), 30, ...
                    'MarkerFaceColor', color, 'MarkerFaceAlpha', alpha, 'MarkerEdgeColor', 'none');
                %contourf(grid2D.xs{1}, grid2D.xs{2}, -vF, -1:1:1);
            end
            colormap('gray');
        end
        
        function updateHumanState(obj, xH)
            obj.xH = xH;
        end
        
        function success = loadFromFile(obj, priorProb, horizon)
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
                obj.data = d.preds;
                obj.hjiGrid = d.predGrid;
                obj.times = d.predTimes;
                obj.tDisc = d.predDt;
                obj.tMin = d.predTMin;
                obj.tMax = d.predTMax;
            end
        end
   end
end