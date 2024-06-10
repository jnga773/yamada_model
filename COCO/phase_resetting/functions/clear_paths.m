function clear_paths()
  % Clear all added folders from path
  rmpath('./boundary_conditions/');
  rmpath('./boundary_conditions/symbolic/');

  rmpath('./continuation_scripts/');
  rmpath('./continuation_scripts/floquet_bundle/');
  rmpath('./continuation_scripts/initial_periodic_orbit/');
  rmpath('./continuation_scripts/isochrons/');
  rmpath('./continuation_scripts/phase_reset/');

  rmpath('./functions/');
  rmpath('./functions/symbolic/');
  
  rmpath('./plotting_scripts/');
  rmpath('./plotting_scripts/initial_periodic_orbit/');
  rmpath('./plotting_scripts/isochrons/');
  rmpath('./plotting_scripts/phase_reset/');
end