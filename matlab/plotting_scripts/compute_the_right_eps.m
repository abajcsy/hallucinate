function [opt_eps, P, X, Y] = compute_the_right_eps(pred2D, predictor, preds, eps_type, delta_reachability)
    
    % Grid for Bayesian prediction
    [X, Y] = predictor.getLatticeMeshgrid();
    
    % Compute the epsilon heurictically
    if strcmp(eps_type, 'heuristic')
      opt_eps = (delta_reachability)^length(preds(:, 1));
    end
    
    % Compute the epsilon heurictically
    if strcmp(eps_type, 'same')
      opt_eps = delta_reachability;
    end
   
    % Compute the epsilon heurictically
    if strcmp(eps_type, 'discard_unlikely_states')
      valid_indices = find(preds(end, :) > 0);
      valid_data = preds(end, valid_indices);
      sorted_valid_data = sort(valid_data, 'descend');
      eps_index = find(cumsum(sorted_valid_data) > (1 - delta_reachability), 1, 'first');
      opt_eps = sorted_valid_data(eps_index);
    end 
    
    % Compute the epsilon based on same area
    if strcmp(eps_type, 'same_area')
      % Find the number of grid points in the reachability predictior
      num_valid_pts = length(find(pred2D <= 0));
      num_eps_searches_per_direction = 30;
      error_so_far = inf;
      for i=-num_eps_searches_per_direction:1:num_eps_searches_per_direction
        eps = (delta_reachability)^(length(preds(:, 1)) + i);
        num_pts = length(find(preds(end, :) >= eps));
        error = abs(num_pts - num_valid_pts);
        if  error < error_so_far
          error_so_far = error
          i
          opt_eps = eps;
        end
      end
    end

    % Compute the optimal predictions
    P = zeros(size(X));
    for k = 1:predictor.numRows
        for l = 1:predictor.numCols
            linIdx = (k-1)*predictor.numCols + l;
            P(k, l) = 1*(preds(end, linIdx) >= opt_eps) + 0*(preds(end, linIdx) < opt_eps);
        end
    end
end