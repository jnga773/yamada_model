function data_out = read_data_lins(run_in, label_in, read_lingap)
  % data_out = read_data_lins(run_in, label_in)
  %
  % Reads the chart data from the previous solution [label_in] of [run_in],
  % grabs the unstable and stable eigenvectors and eigenvalues, and outputs
  % them in arrays.
  %
  % Parameters
  % ----------
  % run_in : str
  %     The string identifier for the previous COCO run that we will
  %     read information from.
  % label_in : int
  %     The integer solution label from the previous run.
  % read_lingap : boolean
  %     Logical check to read lingap data or not
  %
  % Returns
  % -------
  % data_out: data structure
  %     Data structure containing the previous run's values of eps1, esp2,
  %     theta, and the unstable and stable eigenvectors and eigenvalues.
  %         - normal : Normal vector at the intersection point.
  %         - pt0 : Intersection point for hyperplane \Sigma.
  %         - epsilon : Epsilon vector [eps1; eps2; theta].
  %         - vu : Unstable eigenvector.
  %         - vs1 : Stable eigenvector 1.
  %         - vs2 : Stable eigenvector 2.
  %         - lu : Unstable eigenvalue.
  %         - ls1 : Stable eigenvalue 1.
  %         - ls2 : Stable eigenvalue 2.
  %         - lingap : Updated value of the Lin gap (distance).
  %
  % See Also
  % --------
  % coco_read_solution, coll_read_solution

  %-------------------%
  %     Arguments     %
  %-------------------%
  arguments
    run_in char
    label_in double
    read_lingap logical = true
  end

  %----------------------------------%
  %     Read Data: Epsilon Stuff     %
  %----------------------------------%
  % Extract stored deviations of stable and unstable manifolds from
  % stationary equilibrium point
  [datae, charte] = coco_read_solution('bcs_initial', run_in, label_in);
  epsilon = charte.x(datae.epsilon_idx);

  % Normal vector
  normal = datae.normal;
  % Intersection point
  pt0    = datae.pt0;

  %--------------------------------%
  %     Read Data: Eigen Stuff     %
  %--------------------------------%
  % Unstable eigenvector and eigenvalue
  [datau, chartu] = coco_read_solution('bcs_eig_unstable', run_in, label_in);
  vu = chartu.x(datau.vu_idx);
  lu = chartu.x(datau.lu_idx);

  % Stable eigenvector and eigenvalue 1
  [datas1, charts1] = coco_read_solution('bcs_eig_stable1', run_in, label_in);
  vs1 = charts1.x(datas1.vs1_idx);
  ls1 = charts1.x(datas1.ls1_idx);

  % Stable eigenvector and eigenvalue 2
  [datas2, charts2] = coco_read_solution('bcs_eig_stable2', run_in, label_in);
  vs2 = charts2.x(datas2.vs2_idx);
  ls2 = charts2.x(datas2.ls2_idx);

  %----------------------------------%
  %     Read Data: Epsilon Stuff     %
  %----------------------------------%
  if read_lingap
    % Read Linsgap data structure
    data_lingap = coco_read_solution('lins_data', run_in, label_in);

    % Extract stored lingap value from previous run
    [datal, chartl] = coco_read_solution('bcs_lingap', run_in, label_in);
    lingap = chartl.x(datal.lingap_idx);
  end
  
  %----------------%
  %     Output     %
  %----------------%
  % Normal vector
  data_out.normal  = normal;
  % Intersection point
  data_out.pt0     = pt0;
  % Epsilon things
  data_out.epsilon = epsilon;

  % Eigenvectors
  data_out.vu      = vu;
  data_out.vs1     = vs1;
  data_out.vs2     = vs2;

  % Eigenvalues
  data_out.lu      = lu;
  data_out.ls1     = ls1;
  data_out.ls2     = ls2;

  % Lingap
  if read_lingap
    data_out.vgap    = data_lingap.vgap;
    data_out.lingap0 = data_lingap.lingap0;
    data_out.lingap  = lingap;
  end

end