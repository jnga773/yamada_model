function write_bifurcation_data(run_names_in)
  % WRITE_BIFURCATION_DATA: Writes the gamma and A data for each of the
  % different bifucations calculated to .txt files in the ./data/ directory.
  
  %-------------------%
  %     Read Data     %
  %-------------------%
  % Run names
  H_run = run_names_in.hopf_bifurcations;
  S_run = run_names_in.saddle_nodes;
  T_run = run_names_in.transcritical;

  % L_run = run_names_in.approx_homo.continue_homoclinics;
  L_run = run_names_in.lins_method.continue_homoclinics;

  D_run = run_names_in.limit_cycle.follow_limit_cycle;

  % Read COCO data matrices
  bd_H = coco_bd_read(H_run);
  bd_S = coco_bd_read(S_run);
  bd_T = coco_bd_read(T_run);
  bd_L = coco_bd_read(L_run);
  bd_D = coco_bd_read(D_run);
  
  % Hopf bifurcation line (H)
  A_run8 = coco_bd_col(bd_H, 'A');
  gamma_run8 = coco_bd_col(bd_H, 'gamma');

  % Find the minimum to split into H line and NSA line
  [~, idx] = min(A_run8);
  A_NSA = A_run8(1:idx); gamma_NSA = gamma_run8(1:idx);
  A_H = A_run8(idx+1:end); gamma_H = gamma_run8(idx+1:end);

  % Saddle-Node bifurcation line (A_S)
  A_SN = coco_bd_col(bd_S, 'A');
  gamma_SN = coco_bd_col(bd_S, 'gamma');

  % Transcritical bifurcation line (A_T)
  A_T = coco_bd_col(bd_T, 'A');
  gamma_T = coco_bd_col(bd_T, 'gamma');

  % Approximate homoclinic line
  A_homoclinic = coco_bd_col(bd_L, 'A');
  gamma_homoclinic = coco_bd_col(bd_L, 'gamma');

  % Approximate double limit cycle line
  A_double_limit = coco_bd_col(bd_D, 'A');
  gamma_double_limit = coco_bd_col(bd_D, 'gamma');

  %--------------------------------------------%
  %     Write Bifurcation Data to txt File     %
  %--------------------------------------------%
  % % Hopf bifurcation
  % fileID = fopen('./data/H.txt', 'w');
  % fprintf(fileID, '%.12E        %.12E \n', [A_H; gamma_H]);

  % % Neutral Saddle
  % fileID = fopen('./data/SQ.txt', 'w');
  % fprintf(fileID, '%.12E        %.12E \n', [A_NSA; gamma_NSA]);

  % % Saddle-Node
  % fileID = fopen('./data/S.txt', 'w');
  % fprintf(fileID, '%.12E        %.12E \n', [A_SN; gamma_SN]);

  % % Transcritical
  % fileID = fopen('./data/T.txt', 'w');
  % fprintf(fileID, '%.12E        %.12E \n', [A_T; gamma_T]);

  % % Approximate homoclinic
  % fileID = fopen('./data/L.txt', 'w');
  % fprintf(fileID, '%.12E        %.12E \n', [A_homoclinic; gamma_homoclinic]);

  % % Approximate double limit
  % fileID = fopen('./data/D.txt', 'w');
  % fprintf(fileID, '%.12E        %.12E \n', [A_double_limit; gamma_double_limit]);

  % Save to data matrix
  save('./data/Bifurcation_Data.mat', ...
       'A_H', 'gamma_H', ...
       'A_NSA', 'gamma_NSA', ...
       'A_SN', 'gamma_SN', ...
       'A_T', 'gamma_T', ...
       'A_homoclinic', 'gamma_homoclinic', ...
       'A_double_limit', 'gamma_double_limit');
end