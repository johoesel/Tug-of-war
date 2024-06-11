function [Extraction_Probability,bond_coordinates,tvec,Theta,avgForce_on_r] = antigen_extract_dynamics(Tspan,dt,x_aF,x_bF,DeltaG_a,DeltaG_b,gamma_a,gamma_b,gamma_r,KB,T,Forces,Diff,Runs)

N = Tspan/dt; %no. of steps
%Extraction_Probability = zeros(1,length(Forces));
Extraction_Probability = [];
avgForce_on_r = zeros(1,length(Diff));
for D = 1:length(Diff)
    Force_on_r = 0;
for force = 1:length(Forces)
    Extraction_success = 0;
    for run = 1:Runs
        %Sample random Gaussian White Noise
        n1 = randn(1,N+1);
        n2 = randn(1,N+1);
        n3 = randn(1,N+1);
        % n1 = zeros(1,N);
        % n2 = zeros(1,N);

        %Initialise Bond Parameters
        x_a = 0; %total bond extension (noise & angle)
        x_b = 0;
        x = 0; %bond extension due to noise
        y = 0;
        r = 0; %displacement of tether
        theta = 0;
        Theta = 0;
        bond_coordinates = [0, 0, 0];
        %Equilibrium Bond Lengths
        x0 = 0.1e-9; 
        y0 = 0.1e-9;
        
        
        j = 0;
        
        while x_a <= x_aF && x_b <= x_bF
    
            j = j + 1;
            %Finding Derivatives on current time step
            [dU_a,dU_b,dU_theta,dtheta_r] = Potentials(x,y,x0,y0,x_aF,x_bF,DeltaG_a,DeltaG_b,theta,r);
            
            %Solving Langevin Equation (Forward Euler)
            new_x_a = x_a + dt/gamma_a * (-dU_a + dU_b) + sqrt(2*dt*KB*T/gamma_a) * n1(j);
            new_x_b = x_b + x_a - new_x_a + dt/gamma_b * (-dU_b + cos(theta) * Forces(force)) + sqrt(2*dt*KB*T/gamma_b) * n2(j);
            new_r = r + dt/gamma_r * (-dU_theta * dtheta_r * sin(theta)) + sqrt(4*Diff(D)*dt) * n3(j);
            
            

            theta = sin(new_r/(new_x_a + new_x_b + x0 + y0));
            x = cos(theta) * (new_x_a + x0) - x0;
            y = cos(theta) * (new_x_b + y0 + new_x_a) - y0;
            
            %storing coordinates
            bond_coordinates = [bond_coordinates; new_x_a, new_x_b, new_r];
            Theta = [Theta, theta];
            
            %End simulation if taking too long
            if j > N - 1
                disp('particle did not escape the well')
                break
            end

            x_a = new_x_a;
            x_b = new_x_b;
            r = new_r;
        
        end
        %Find the final force on r for each run
        Force_on_r = [Force_on_r, -dU_theta * dtheta_r * sin(theta)];
        %Find Extraction Probabilities
        if x_a >= x_aF
            Extraction_success = Extraction_success + 1;
        end
        tvec = (0:j)*dt;
    end
    %Average maximal restoring forces dependent on diffusion coefficient
    avgForce_on_r(D) =  mean(Force_on_r);
    %disp(force)
    if length(Forces)>1
        disp(force)
    else
        disp(D)
    end
    %Extraction_Probability(force) = Extraction_success/Runs;
    Extraction_Probability = [Extraction_Probability Extraction_success/Runs];
end
end
end