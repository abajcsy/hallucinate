classdef OccupancyGrid < handle
   properties
       xDisc
       yDisc

       xMin
       xMax
       yMin
       yMax

       data
   end

   methods
       function obj = OccupancyGrid(xDisc, yDisc, xMin, xMax, yMin, yMax)
           obj.xDisc = xDisc;
           obj.yDisc = yDisc;

           obj.xMin = xMin;
           obj.xMax = xMax;
           obj.yMin = yMin;
           obj.yMax = yMax;

           xExtent = ceil((xMax - xMin) / xDisc);
           yExtent = ceil((yMax - yMin) / yDisc);

           obj.data = zeros(yExtent, xExtent);
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