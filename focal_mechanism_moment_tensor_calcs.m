


% Calculate the moment tensors for a few 2-D focal mechanisms that I want
% to model in Specfem:

% MW - normal fault
strike = 0;
dip = 45;
rake = 90;
mt1 = sdr2mt(strike, dip, rake)

% IF - low-angle thrust
strike = 0;
dip = 21;
rake = -90;
mt2 = sdr2mt(strike, dip, rake)

% Slab crust - perpendicular thrust
strike = 90;
dip = 45;
rake = -90;
mt3 = sdr2mt(strike, dip, rake)

% output is 6 values:
%     [Mrr Mtt Mpp Mrt Mrp Mtp] in Up - South - East coordinate system
%      Mrr Mrt Mrp
%          Mtt Mtp
%              Mpp