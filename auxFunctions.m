classdef    auxFunctions 
    methods     ( Static = true )
        % Function to calculate matrices
        function output = calculateMatrices(struct)
            
            % Initialize J H and D
            J = eye(3);
            D = NaN(3,length(struct.beta1)-1);
            H = NaN(3,1);
            % First regression
            H(1)    = struct.beta1(1);
            D(1,:)  = struct.beta1(2:end);
            % Second regression
            J(2,1)  = -struct.beta2(2);
            H(2)    = struct.beta2(1)+struct.beta2(1)*H(1);
            D(2,:)  = struct.beta2(3:end) + struct.beta2(2)*(D(1,:)');
            % Third Regression
            J(3,1:2) = struct.beta3(2:3);
            H(3)     = struct.beta3(1) - J(3,1:2)*H(1:2);
            D(3,:)   = struct.beta3(4:end)'+sum(D(1:2,:).*struct.beta3(2:3));
            % F matrix
            Delta    = (1/(struct.n))*diag([struct.sigma1,struct.sigma2,struct.sigma3]);
            F        = J*sqrt(Delta);
            % Store results
            output.H = H;
            output.D = D;
            output.J = J;
            output.F = F;
            
        end
        % Function to Update up to date Gamma
        function Gamma = calculateGamma(Rt,Gamma0)
            % Initial Prior
            Gamma = Gamma0;
            % Loop to update Gamma
            for i = 1:size(Rt,2)
                Gamma = (Rt(i,:)')*Rt(i,:) + Gamma;
            end
            
        end
        % Function to calculate A
        function output = calculateVAR(struct)
            output    = struct;
            output.A0 = zeros(14,1);
            output.A  = zeros(14);
            output.B  = zeros(14,3);
            % Fill in the ones
            for j = 1:13
                output.A(j+1,j) = 1;
            end
            % Fill the remaining
            idx = [1,5,10]';
            for j = 1:3
                output.A0(idx(j))  = struct.H(j);
                output.A(idx(j),:) = struct.D(j,:);
                output.B(idx(j),:) = struct.F(j,:);
            end
            
        end
        % Function to calculate the IRF
        function irf = calculateIRF(struct,varargin)
            % Define shock
            W = [1,0,0]*(struct.F+struct.D*((eye(14)-struct.A)\struct.B));
            % Normalize shock
            W = W'/norm(W);
            % Check number of arguments
            if(nargin==1)
                % Long term impulse response
                irf =[1,0,0]*(struct.F+struct.D*((eye(14)-struct.A)\struct.B))*W;
            else
                % Impulse response at time tau
                tau = varargin{1};
                % Define the IRF
                irf = struct.F;
                if(tau>1)
                    for j=1:(tau-1)
                        irf = irf + struct.D*(struct.A^(j-1))*struct.B;
                    end
                end
                % Multiply by the shock
                irf = [1,0,0]*irf*W;
                
            end
            
            
        end
        
    end
end