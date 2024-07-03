%-------------------------------------------------------------------------%
%%                               READ DATA                               %%
%-------------------------------------------------------------------------%
run_name = 'run10_phase_reset_PTC_multi';

% Read data
[theta_old_read, theta_new_read, A_perturb_read] = read_PTC_scan_data(run_name);

% Fix data
[theta_old_fix, theta_new_fix] = fix_data(theta_old_read, theta_new_read);

% Pad arrays
[theta_old_gt1, theta_new_gt1, A_perturb_gt1] = pad_arrays(theta_old_fix.gt1, theta_new_fix.gt1, A_perturb_read);
[theta_old_lt1, theta_new_lt1, A_perturb_lt1] = pad_arrays(theta_old_fix.lt1, theta_new_fix.lt1, A_perturb_read);

% Save to .mat file
save('./PTC_scan_test.mat');

%-------------------------------------------------------------------------%
%%                               FUNCTIONS                               %%
%-------------------------------------------------------------------------%
function [theta_old_out, theta_new_out, A_perturb_out] = read_PTC_scan_data(run_in)
  % Reads the data innint

  % Folder name
  data_dir = sprintf('./data/%s/', run_in);
  % List all directories
  dirs = dir(data_dir);
  % Remove ./ and ../
  dirs = dirs(~ismember({dirs.name}, {'.', '..', '.DS_Store'}));
  % Sub folder names
  dir_sub = {dirs.name};

  % Empty arrays for theta_old and theta_new
  theta_old_gt1 = {};
  theta_old_lt1 = {};
  theta_new_gt1 = {};
  theta_new_lt1 = {};
  A_perturb_out = {};

  for i = 1 : length(data_dir)
    % Sub directory name
    dir_read = {run_in, dir_sub{i}};

    % Bifurcation data
    bd_read = coco_bd_read(dir_read);

    %-------------------%
    %     Read Data     %
    %-------------------%
    % theta_old
    theta_old = coco_bd_col(bd_read, 'theta_old');
    % theta_new
    theta_new = coco_bd_col(bd_read, 'theta_new');
    % A_perturb
    A_perturb_read = coco_bd_val(bd_read, 1, 'A_perturb');

    %--------------------%
    %     Split Data     %
    %--------------------%
    % Create mask for theta_old > 1
    mask_gt1 = theta_old > 1.0;
    % Create mask for theta_old <= 1
    mask_lt1 = theta_old <= 1.0;

    % Split up theta_old data
    old_gt1_read = theta_old(mask_gt1);
    old_lt1_read = theta_old(mask_lt1);
    % Split up theta_new data
    new_gt1_read = theta_new(mask_gt1);
    new_lt1_read = theta_new(mask_lt1);

    %---------------------%
    %     Append Data     %
    %---------------------%
    % theta_old
    theta_old_gt1{i} = old_gt1_read;
    theta_old_lt1{i} = old_lt1_read;
    % theta_new
    theta_new_gt1{i} = new_gt1_read;
    theta_new_lt1{i} = new_lt1_read;
    % A_perturb
    A_perturb_out{i} = A_perturb_read;

  end

  %----------------%
  %     Output     %
  %----------------%
  % Make data structure of theta_old and theta_new
  theta_old_out.gt1 = theta_old_gt1;
  theta_old_out.lt1 = theta_old_lt1;
  theta_new_out.gt1 = theta_new_gt1;
  theta_new_out.lt1 = theta_new_lt1;

end

function [theta_old_out, theta_new_out] = fix_data(theta_old_in, theta_new_in)
  % Fixes the data and moves things around

  % Empty arrays
  theta_old_gt1 = {};
  theta_old_lt1 = {};
  theta_new_gt1 = {};
  theta_new_lt1 = {};

  % Fix up data for i in range(len(theta_old))
  for i = 1 : length(theta_old_in.gt1)
    % theta_old values
    old_gt1 = theta_old_in.gt1{i} - 1.0;
    old_lt1 = theta_old_in.lt1{i};
    % theta_new values
    new_gt1 = theta_new_in.gt1{i};
    new_lt1 = theta_new_in.lt1{i};

    % Check if theta_new goes below 0 at all
    if min(new_gt1) < 0.0
      new_gt1 = new_gt1 + 1.0;
    end
    if min(new_lt1) < 0.0
      new_lt1 = new_lt1 + 1.0;
    end

    %---------------------%
    %     Append Data     %
    %---------------------%
    % theta_old
    theta_old_gt1{i} = old_gt1;
    theta_old_lt1{i} = old_lt1;
    % theta_new
    theta_new_gt1{i} = new_gt1;
    theta_new_lt1{i} = new_lt1;

  end
  
  %----------------%
  %     Output     %
  %----------------%
  % Make data structure of theta_old and theta_new
  theta_old_out.gt1 = theta_old_gt1;
  theta_old_out.lt1 = theta_old_lt1;
  theta_new_out.gt1 = theta_new_gt1;
  theta_new_out.lt1 = theta_new_lt1;

end

function [theta_old_out, theta_new_out, A_perturb_out] = pad_arrays(theta_old_in, theta_new_in, A_perturb_in)
  % Pads the theta_old and theta_new arrays to match the arrays with the max
  % length. This makes things easier to plot using MATLAB's surf.

  % cycle through and calculate max length
  array_lengths = zeros(length(theta_old_in), 1);
  for i = 1 : length(theta_old_in)
    array_lengths(i) = length(theta_old_in{i});
  end

  % find maximum array length to pad to
  max_len = max(array_lengths);

  % Empty arrays
  theta_old_out = [];
  theta_new_out = [];
  A_perturb_out = [];

  % Cycle through data and pad
  for i = 1 : length(theta_old_in)
    %-------------------%
    %     Read Data     %
    %-------------------%
    A_read   = A_perturb_in{i};
    old_read = theta_old_in{i};
    new_read = theta_new_in{i};

    %------------------%
    %     Pad Data     %
    %------------------%
    % Pad length
    pad_len = max_len - length(old_read);

    % A_perturb values
    A_pad = A_read * ones(1, max_len);

    % Pad arrays
    old_pad = padarray(old_read, [0, pad_len], old_read(1), 'pre');
    new_pad = padarray(new_read, [0, pad_len], new_read(1), 'pre');

    %---------------------%
    %     Append Data     %
    %---------------------%
    theta_old_out = [theta_old_out; old_pad];
    theta_new_out = [theta_new_out; new_pad];
    A_perturb_out = [A_perturb_out; A_pad];

  end

  %----------------%
  %     Output     %
  %----------------%
  theta_old_out = theta_old_out';
  theta_new_out = theta_new_out';
  A_perturb_out = A_perturb_out';

end