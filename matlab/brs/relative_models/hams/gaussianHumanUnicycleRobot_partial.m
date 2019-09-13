function alpha = gaussianHumanUnicycleRobot_partial( ...
  t, data, derivMin, derivMax, schemeData, dim)

checkStructureFields(schemeData, 'dynSys');

switch dim
  case 1
    alpha = schemeData.dynSys.v + schemeData.dynSys.vRRange(2);
    
  case 2
    alpha = schemeData.dynSys.v + schemeData.dynSys.vRRange(2);
    
  case 3
    alpha = schemeData.dynSys.gamma;
  
  case 4
    alpha = schemeData.dynSys.wRRange(2);  
end

end