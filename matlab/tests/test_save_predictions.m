%clear all
%clc

path = '/home/abajcsy/hybrid_ws/src/hallucinate/frs/';

filename = 'pb0.80_xval_0.00_yval_0.00_occuMap.mat';

% Grab:
%   predictions
%   predGrid
%   'predTimes', 'predDt', 'predTmin', 'predTmax
load(strcat(path, filename));

predT1 = predictions(:,:,1);
predT60 = predictions(:,:,60);
predT124 = predictions(:,:,124);

%figure(1)
%contourf(predGrid.xs{1}, predGrid.xs{2}, predT1, [0,0.1]);

%figure()
%contourf(predGrid.xs{1}, predGrid.xs{2}, predT60, [0,0.1]);

hold on
[M,c2] = contour(predGrid.xs{1}, predGrid.xs{2}, predT124, [0,0.1]);
c2.LineWidth = 2;
c2.Color = 'r';