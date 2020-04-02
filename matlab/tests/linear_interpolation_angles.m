clear all;

numCtrls = 4;

% State space is 2 --> x = 0 or x = 1. 
% lowerBound and upperBound is Nx=2 size
%   stores the lower and upper control bounds for each state.
lowerBound = [-pi/4, -3*pi/4];
upperBound = [pi/4, 3*pi/4];

likelyCtrls = cell(1, numCtrls);

% Minimum angular distance between the control bounds. 
diff = angdiff(lowerBound, upperBound);

% direction to traverse unit circle when linearly interpolating.
reverse_dir = sign(diff); 

% increment to add to control.
incr = abs(diff)/(numCtrls-1);

% generate controls \in [lowerBound, upperBound]
for i=1:numCtrls
    likelyCtrls{i} = ((i-1)*incr + lowerBound) .* (reverse_dir > 0) + ...
                        wrapToPi(lowerBound - (i-1)*incr) .* (reverse_dir < 0);
    
    fprintf(strcat("u",num2str(i),"(x=0)=",num2str(likelyCtrls{i}(1)),...
        ", u",num2str(i),"(x=1)=",num2str(likelyCtrls{i}(2)),"\n"));
end