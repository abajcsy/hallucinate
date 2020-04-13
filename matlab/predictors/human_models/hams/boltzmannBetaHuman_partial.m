function alpha = boltzmannBetaHuman_partial( ...
  t, data, derivMin, derivMax, schemeData, dim)

checkStructureFields(schemeData, 'dynSys');

if dim == 1
    alpha = schemeData.dynSys.v;
elseif dim == 2
    alpha = schemeData.dynSys.v;
else
    alpha = schemeData.dynSys.gamma;
end

end