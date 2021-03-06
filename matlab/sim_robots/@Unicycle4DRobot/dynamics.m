function dx = dynamics(obj, ~, x, u)
    % dx = dynamics(obj,x, u)
    %         \dot{x}_1 = x_4 * cos(x_3) + d_1
    %         \dot{x}_2 = x_4 * sin(x_3) + d_1
    %         \dot{x}_3 = u_1 = u_1
    %         \dot{x}_4 = u_2 = u_2

    dx = cell(obj.nx, 1);

    returnVector = false;
    if ~iscell(x)
      returnVector = true;
      x = num2cell(x);
      u = num2cell(u);
    end

    for i = 1:length(obj.dims)
      dx{i} = dynamics_i(x, u, obj.vRange, obj.dims, obj.dims(i));
    end

    if returnVector
      dx = cell2mat(dx);
    end
end

function dx = dynamics_i(x, u, vRange, dims, dim)

    switch dim
      case 1
        dx = x{dims==4} .* cos(x{dims==3});
      case 2
        dx = x{dims==4} .* sin(x{dims==3});
      case 3
        dx = u{1};
      case 4
        dx = (x{dims==4} >= vRange(2)) .* min(u{2}, 0.0) + ...
          (x{dims==4} <= vRange(1)) .* max(u{2}, 0.0) + ...
          ((x{dims==4} > vRange(1)) .* (x{dims==4} < vRange(2))) .* u{2};
      otherwise
        error('Only dimension 1-4 are defined for dynamics of Plane4D!')    
    end
end
