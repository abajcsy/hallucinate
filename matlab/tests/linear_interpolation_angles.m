clear all;

numCtrls = 4;

% State space is 3 --> x = {0,1,2,3}. 
% lowerBound and upperBound is Nx=4 size
%   stores the lower and upper control bounds for each state.
lowerBound = [-pi/4, -3*pi/4, pi/4, -pi];
upperBound = [pi/4, 3*pi/4, pi/4, pi];

likelyCtrls = cell(1, numCtrls);

tol = 1e-6;

% Minimum angular distance between the control bounds. 
diff = angdiff(lowerBound, upperBound);

% need to account for sneaky boundary condition when lower=-pi and upper bound=pi 
boundary_cond_indicies = find(abs(lowerBound - (-pi)) < tol & abs(upperBound - pi) < tol);
diff(boundary_cond_indicies) = abs(upperBound(boundary_cond_indicies) - lowerBound(boundary_cond_indicies));

% direction to traverse unit circle when linearly interpolating.
dir = sign(diff);

% increment to add to control.
incr = abs(diff)/(numCtrls-1);

% generate controls \in [lowerBound, upperBound]
for i=1:numCtrls
    likelyCtrls{i} = ((i-1)*incr + lowerBound) .* (dir > 0) + ...
                        wrapToPi(lowerBound - (i-1)*incr) .* (dir < 0) + ...
                        lowerBound .* (dir == 0);
                    
    str1 = strcat("u",num2str(i),"(x=0)=",num2str(likelyCtrls{i}(1)),", ");
    str2 = strcat("u",num2str(i),"(x=1)=",num2str(likelyCtrls{i}(2)),", ");
    str3 = strcat("u",num2str(i),"(x=2)=",num2str(likelyCtrls{i}(3)),", ");
    str4 = strcat("u",num2str(i),"(x=3)=",num2str(likelyCtrls{i}(4)));
    fprintf(strcat(str1, str2, str3, str4, "\n"));
end