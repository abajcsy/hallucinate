classdef OccupancyGrid < handle
   properties
       xDisc
       yDisc
       tDisc

       xMin
       xMax
       yMin
       yMax

       data
       times
   end

   methods
       function obj = OccupancyGrid(xDisc, yDisc, tDisc, xMin, xMax, yMin, yMax)
           obj.xDisc = xDisc;
           obj.yDisc = yDisc;
           obj.tDisc = tDisc;

           obj.xMin = xMin;
           obj.xMax = xMax;
           obj.yMin = yMin;
           obj.yMax = yMax;

           %xExtent = ceil((xMax - xMin) / xDisc);
           %yExtent = ceil((yMax - yMin) / yDisc);

           %obj.data = zeros(yExtent, xExtent);
       end

       function fromValueFuns(obj, grid, valueFuns, times, t0)
           obj.times = zeros(1,length(times));
           obj.data = cell(1,length(times));
           
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
               obj.data = valueFuns;
               for i=1:length(times)
                obj.times(i) = t0 + times(i);
               end
           end
       end
       
       function setData(obj, x, y, val)
           [i, j] = obj.xyToIndex(x, y);
           obj.data(i, j) = val;
       end

       function val = getData(obj, x, y)
           [i, j] = obj.xyToIndex(x, y);
           val = obj.data(i, j);
       end

       function [i, j] = xyToIndex(obj, x, y)
           i = floor((y - obj.yMin) / obj.yDisc) + 1;
           j = floor((x - obj.xMin) / obj.xDisc) + 1;
       end

       function [x, y] = indexToXY(obj, i, j)
           x = obj.xMin + (j - 1) * obj.xDisc;
           y = obj.yMin + (i - 1) * obj.yDisc;
       end
      
   end
end