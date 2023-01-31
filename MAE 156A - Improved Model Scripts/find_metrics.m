function [tr_exp, wterm_exp] = find_metrics(omega,t)

% Experimental Terminal Velocity
wterm_exp= mean(omega(end-20:end));% averaging the last 100 data point of a 10sec run [RPM]

% Experimental Rise Time
tr_exp = t(find(omega>=0.63*wterm_exp,1)); % [s]

end


















