function [dU_a,dU_b,dU_theta,dtheta_r,U_a,U_b] = Potentials(x,y,x0,y0,x_aF,x_bF,DeltaG_a,DeltaG_b,theta,r)

%Total bond extensions
x_a = (x0 + x)/cos(theta) - x0;
x_b = (y0 + y)/cos(theta) - y0 - x_a;

%% Forces
dU_a = 6 * DeltaG_a * x_a .* (x_aF - x_a) / x_aF^3;
dU_b = 6 * DeltaG_b * x_b .* (x_bF - x_b) / x_bF^3;
%from wolframalpha
dU_theta = 6 * DeltaG_a .* (x0 + x)*tan(theta).*sec(theta).*((x0+x)*sec(theta)-x0) .* (x_aF - (x0+x).*sec(theta)+ x0) / x_aF^3 +...
           6 * DeltaG_b .* (y0 + y)*tan(theta).*sec(theta).*((y0+y)*sec(theta)-y0) .* (x_bF - (y0+y)*sec(theta)+ y0) / x_bF^3;
dtheta_r = (y + y0)/((y+y0).^2 + r^2);


%% Potentials
U_a = 3/2 * DeltaG_a * (x_a./x_aF - 1/2) - 2 * DeltaG_a * (x_a./x_aF - 1/2).^3;
U_b = 3/2 * DeltaG_b * (x_b./x_bF - 1/2) - 2 * DeltaG_b * (x_b./x_bF - 1/2).^3;
U = U_a + U_b;


% wolfram d/dx (3/2 *G *(((x+o)/cos(t)-o)/F-1/2) - 2 *G *(((x+o)/cos(t)-o)/F-1/2)^3)
% wolfram d/dt (3/2 *G *(((x+o)/cos(t)-o)/F-1/2) - 2 *G *(((x+o)/cos(t)-o)/F-1/2)^3)
end       