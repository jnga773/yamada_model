clc; clear all;

%-------------------%
%     Read Data     %
%-------------------%
load('./data_mat/floquet_solution.mat', 'gamma_read', 'tbp_read', 'T_read');

% Number of periods
k = 25;

%----------------------------------------------%
%     Segment Initial Conditions: Periodic     %
%----------------------------------------------%
t1 = tbp_read;
x1 = gamma_read;

% Otherwise, keep appending periodic solutions
if k > 1
% Cycle through k integers
for j = 1 : k-1
  % Append another period of time data
  t1 = [tbp_read    ; max(tbp_read) + t1(2:end)];
  
  x1 = [gamma_read; x1(2:end, :)];

end
% Normalise time data by integer
t1 = t1 / k;

end

%----------------------------------------------%
%     Segment Initial Conditions: Periodic     %
%----------------------------------------------%
% Segment 4
% If only one period, i.e., k = 1, then the
% input solutions remain unchanged
t2 = tbp_read;
x2 = gamma_read;
t_max  = max(t2);

for i = 1 : k-1
  % Append x_PO
  x2 = [x2; gamma_read(2:end, :)];

  % Append t_PO
  t2 = [t2; t_max + tbp_read(2:end)];

  t_max = t2(end);
end

% Normalise by k
t2 = t2 / k;