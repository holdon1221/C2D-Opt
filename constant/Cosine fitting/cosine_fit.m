function cosine_fit

clear all


load E2_circ.mat E2_circ                   % load extracted hormonal circadian rhythm data E2_circ, P4_circ, LH_circ, or FSH_circ   


time_data = E2_circ(:,1);
data      = E2_circ(:,2);

x1 = (time_data)';                                        
y1 = (data/mean(data))';


cosinor_fit = 'a1*cos((pi/12)*(x - c1)) + d1';
startpoints1 = [0.5 4 1];
        
f1 = fit(x1', y1', cosinor_fit, 'Start', startpoints1)

%plot(f1, x1, y1)

% Convert acrophase c1 from hours to days

