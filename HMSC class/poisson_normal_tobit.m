classdef poisson_normal_tobit
    properties
        %data input and output
        Ey
        y
        sigma2
        %latent variable eta ~ N(Ey, sigma2)
        eta
        %internal variables for Hamiltonian Monte Carlo
        eps % step size in Hamiltonian
        accept_rate % acceptance rate
        %other variables
        n
        m
        pos_set
        zero_set
     end
    
    methods
        function obj = poisson_normal_tobit()
            % constructor
            obj.eps = 0.0001;
        end
        
        function result = truncated_normal_rnd(obj, mu, sigma)
            % a very crude truncated normal generator, only generate positive value
            a = cdf('Normal',0, mu, sigma);
            b = 1;
            
            U = rand(length(mu),1);
            prob = a + U .* (b-a);
            
            result = norminv(prob,mu,sigma);
            
        end
        
        function obj = input_Ey_sigma2(obj, y, Ey, sigma2)
            obj.n = numel(y);
            obj.y = reshape(y, [obj.n,1]);
            obj.eta = ones([obj.n,1]);
            obj.pos_set = obj.y>0;
            obj.zero_set = obj.y==0;
            obj.Ey = reshape(Ey,[obj.n,1]);            
            obj.sigma2 = reshape(sigma2, [obj.n,1]);
        end
        
        %gradient of the - log-likelihood
        function result = grad_U(obj, q)
            result = 1.0 - obj.y(obj.pos_set)./q + (q - obj.Ey(obj.pos_set)) ./ obj.sigma2(obj.pos_set);
        end
        
        %the - log-likelihood
        function result = U(obj, q)
            diff = (q - obj.Ey(obj.pos_set));
            result = q - obj.y(obj.pos_set).*log(q) +  ...
                diff.* diff ./ 2 ./ obj.sigma2(obj.pos_set);
        end
        
        %HMC to propose positive random variable for those with Y>0
        %It uses Hamiltonian dynamics to run freely for L steps, then check if all
        %the proposals are >0, if not, keep runnning a few steps to make
        %sure they are all >0.
        %At the second step, it uses Metropolis to accept the proposal
        %See Neal 2012 and Pakman 2012 for details
        function [result, obj] = HMC_for_positive(obj, epsilon, L, current_q)
            
            q = current_q;
            p = normrnd(0,1 , [length(q),1]);  % independent standard normal variates
            current_p = p;
            
            % Make a half step for momentum at the beginning

          
            
            p = p - epsilon .* grad_U(obj, q) ./ 2;
            
            % Alternate full steps for position and momentum
         
             
            for i= 1:L
                %      Make a full step for the position
                q = q + epsilon .* p;
                
   
           
                %      Make a full step for the momentum, except at end of trajectory
                if (i~=L)
                    p = p - epsilon .* grad_U(obj, q);
                end            
                
                   
            end
            
          
            %    Make a half step for momentum at the end.
            
            p = p - epsilon .* grad_U(obj, q) ./ 2;
            
            
            
            
            %   run a few extra steps to make sure all q are >0
            violated_set = q<= 0 | isinf(q) |isnan(q);

%             loop_count = 0;
            while any(violated_set)
                epsilon = epsilon .* 2;
                
                grad_u_q = grad_U(obj, q);
                p(violated_set) = p(violated_set) - epsilon .* grad_u_q(violated_set)./2;
                q(violated_set) = q(violated_set) + epsilon .* p(violated_set);
                grad_u_q = grad_U(obj, q);
                
                p(violated_set) = p(violated_set) - epsilon .* grad_u_q(violated_set)./2;
                
                violated_set = q<= 0 | isinf(q) |isnan(q);

% Mechanism to prevent infinite looping, do not uncomment unless debugging                 
%                 loop_count = loop_count+1;
%                 if loop_count>100
% %                     obj.y(violated_set)
% %                     p(violated_set)
% %                     grad_U(obj, p(violated_set))
%                     break;
%                 end
            end
             
%             if any(violated_set)
%                 q(violated_set)=current_q(violated_set);
%             end

            
            %   # Negate momentum at end of trajectory to make the proposal symmetric
            
            p = -p;
            
            %   # Evaluate potential and kinetic energies at start and end of trajectory
            
            current_U = U(obj, current_q);
            current_K = (current_p.*current_p) ./ 2;
            proposed_U = U(obj, q);
            proposed_K = (p.*p) ./ 2;
            
            %   # Accept or reject the state at end of trajectory, returning either
            %   # the position at the end of the trajectory or the initial position
            accept_set = rand([length(q),1]) < exp(current_U-proposed_U+current_K-proposed_K);
            %
            result = current_q;
            result(accept_set) = q(accept_set);
            obj.accept_rate = sum(accept_set) ./ length(q);
            
        end
        
        %A two step generator for those Y==0
        %First, a Bernoulli is picked to choose eta>0 or eta<=0
        %according to the P(eta>0)/C vs P(eta<=0)/C, (C is the normalizer)
        %Then, if eta<=0, sample from truncated eta ~ N(Ey, sigma2) 
        %if eta>0, sample eta from the truncated posterior with
        % N(Ey, sigma2) Poi( Ey|Y=0)

        function result = Post_eta_for_zero(obj)
            
            sigma = sqrt(obj.sigma2);
            
            p_left = 1- cdf('normal', obj.Ey(obj.zero_set)./sigma(obj.zero_set) );
            p_right = cdf('normal', obj.Ey(obj.zero_set)./sigma(obj.zero_set) - sigma(obj.zero_set) );
            
            p_adjusted = p_left ./(p_left+ p_right);

            positive = rand(sum(obj.zero_set),1) > p_adjusted;
                   
            result = zeros(sum(obj.zero_set),1);
            
            if sum(~positive)>0
                temp = obj.Ey(obj.zero_set);
                sigma_temp = sigma(obj.zero_set);
                result(~positive) = normrnd(temp(~positive),sigma_temp(~positive));
            end
            
%             positive = result >0;
            
            if sum(positive)>0
                temp = obj.Ey(obj.zero_set);
                pos_mean = temp(positive) - obj.sigma2(positive);
                result(positive) = truncated_normal_rnd(obj, pos_mean, sigma(positive));
                
                new_pos = result(positive);
                sigma_pos = sigma(positive);
                violated= isinf(new_pos) | isnan(new_pos);
                while any(violated)
                       new_pos(violated) = truncated_normal_rnd(obj, pos_mean(violated),sigma_pos(violated) );
                       violated= isinf(new_pos) | isnan(new_pos);
%                        pos_mean(violated)
%                        sigma_pos(violated)
                end
                
                result(positive) = new_pos;

            end
            
        end
        
        
        
        function obj = run(obj, steps, burn_in)
            
            
            obj.eps = 0.0001;
           
            for i = 1:steps
                %use Hamiltonina MC to update the eta corresponding to
                % Y>0
                [temp, obj] =  HMC_for_positive(obj, obj.eps, 100, obj.eta(obj.pos_set));
                obj.eta(obj.pos_set) = temp;
                
                %tune the eps to have around 0.65 acceptance rate
                if(burn_in)
                    obj.eps = obj.eps * exp(obj.accept_rate - 0.65);
                end
%                 obj.accept_rate

                %use two-step normal to sample the eta corresponding to
                %Y==0
%                 temp = Post_eta_for_zero(obj);
%                 obj.eta(obj.zero_set) = temp;
                               
            end
        
        end        
        
    end
    
    
end