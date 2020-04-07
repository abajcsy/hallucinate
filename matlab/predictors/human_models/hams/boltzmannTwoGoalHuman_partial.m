function alpha = boltzmannTwoGoalHuman_partial( ...
  t, data, derivMin, derivMax, schemeData, dim)

checkStructureFields(schemeData, 'dynSys');

switch dim
  case 1
    alpha = schemeData.dynSys.v;
    
  case 2
    alpha = schemeData.dynSys.v;
    
  case 3
    alpha = schemeData.dynSys.gamma;
end

end