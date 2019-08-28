classdef FixedBetaGaussianHuman < SimHuman
    %SIMFIXEDBETAGAUSSIANHUMAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        xcurr
        mu
        sigma
        uRange
        clippedNorm
        beta
    end
    
    methods
        function obj = FixedBetaGaussianHuman(x0, v, mu, sigma, uRange, beta)
            %SIMFIXEDBETAGAUSSIANHUMAN Construct an instance of this class
            %   Detailed explanation goes here
            obj.xcurr = x0;
            obj.v = v;
            obj.mu = mu;
            obj.sigma = sigma;
            obj.uRange = uRange;
            obj.beta = beta;
            
            if obj.beta == 0
                % Simulated a "rational" human. 
                pd = makedist('Normal','mu',obj.mu,'sigma',obj.sigma);

                % Truncate normal distribution to have finite support that is 
                % within the valid control bounds.
                obj.clippedNorm = truncate(pd,uRange(1),uRange(2));
            elseif beta == 1
                % Simulating an "irrational", or truly random human.
                rng(0,'twister');
            else
                error("This human only works with beta = 0 or 1 values!");
            end
            
        end
        
        %% Simulates an action of the human according to the model.
        function [xnext, u] = simulate(obj, x, t, dt)
            if obj.beta == 0
                % Generate random action from the normal distribution. 
                u = random(obj.clippedNorm,1,1);
            else
                % Generate uniformly random action.
                u = unifrnd(obj.uRange(1), obj.uRange(2));
            end
           
            % Figure out where the human would move in one timestep.
            xdot = dynamics(obj, obj.xcurr, u);
            xnext = obj.xcurr + dt*xdot;
            obj.xcurr = xnext;
        end
    end
end

