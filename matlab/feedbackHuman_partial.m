function alpha = feedbackHuman_partial( ...
  t, data, derivMin, derivMax, schemeData, dim)

checkStructureFields(schemeData, 'grid');

g = schemeData.grid;
switch dim
  case 1
    % Control
    alpha = schemeData.dynSys.v;
    
  case 2
    % Control
    alpha = schemeData.dynSys.v;
    
  case 3
    % Control
    alpha = schemeData.dynSys.alpha;
end
end