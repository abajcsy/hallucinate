function alpha = feedbackHuman_partial( ...
  t, data, derivMin, derivMax, schemeData, dim)

checkStructureFields(schemeData, 'dynSys');

switch dim
  case 1
    alpha = schemeData.dynSys.v;
    
  case 2
    alpha = schemeData.dynSys.v;
    
  case 3
    alpha = schemeData.dynSys.alpha;
end
end