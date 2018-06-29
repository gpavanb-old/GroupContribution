classdef GroupContribution
    properties
        % DEFINE RELATIVE PATHS
        solverInputDir = 'data/Solver_Input/'
        compDescDir = 'data/Compound_Descriptions/'
        initDataDir = 'data/Init_Data/'

        numGroups
        num_comp
        
        propMat
        funcMat
        initMat
        
        MW
        L
                
        PcVec
        TcVec
        TbVec
        VcVec
        omegaVec
        kappaVec
        epsVec
        SigmaVec
        a
        
        x_l_init
        
    end
    
    methods
        function obj = GroupContribution(funcName,varargin)
            relDir = '../include/GroupContribution/';
            obj.propMat = xlsread(strcat(relDir,obj.solverInputDir,'gani_prop_table.xlsx'));
            obj.funcMat = xlsread(strcat(relDir,obj.compDescDir,funcName,'.xlsx'));
            
            % EQUALIZE THEIR SIZES BY PADDING FUNCMAT
            obj.funcMat(size(obj.propMat,2)) = 0; 
            
            obj.numGroups = size(obj.propMat,2);
            obj.num_comp = size(obj.funcMat,1);
            
            if length(varargin)==1
                initName = varargin{1};
                [obj.initMat,~] = xlsread(strcat(relDir,obj.initDataDir,initName,'.xlsx'));
                 obj.x_l_init = (obj.initMat(:,1))';
            else
                obj.x_l_init = ones(1,obj.num_comp);
            end
                
            obj.x_l_init = obj.x_l_init/sum(obj.x_l_init);
                 
            obj.MW = (1e-3 * obj.funcMat*obj.propMat(10,:)')';

            % CALCULATE CRITICAL PROPERTIES
            TcCoVec = obj.propMat(1,:);
            PcCoVec = obj.propMat(2,:);
            VcCoVec = obj.propMat(3,:);
            omegaCoVec = obj.propMat(4,:);
            TbCoVec = obj.propMat(13,:);
            
            % CALCULATE CRITICAL VOLUME VECTOR
            obj.VcVec = 1e-3*(-0.00435 + ( (obj.funcMat*VcCoVec')' ));
            
            % CALCULATE CRITICAL TEMPERATURE VECTOR
            obj.TcVec = 181.128*log( (obj.funcMat*TcCoVec')' );
            
            % CALCULATE ATMOSPHERIC BOILING VECTOR
            obj.TbVec = 204.359*log( (obj.funcMat*TbCoVec')');
            
            % CALCULATE CRITICAL PRESSURE VECTOR
            obj.PcVec = 1.3705 + ( (obj.funcMat*PcCoVec')' + 0.10022).^(-2);
            obj.PcVec = obj.PcVec * 101325;
            
            % CALCULATE ACENTRIC FACTOR VECTOR
            obj.omegaVec = 0.4085 * (log( (obj.funcMat*omegaCoVec')' + 1.1507).^(1/0.5050) );
            
            obj.kappaVec = zeros(1,length(obj.omegaVec));
            for i = 1:length(obj.omegaVec)
                if (obj.omegaVec(i) <= 0.49)
                    obj.kappaVec(i) = 0.37464 + 1.54226*obj.omegaVec(i) - 0.26992*(obj.omegaVec(i)^2);
                else
                    obj.kappaVec(i) = 0.379642 + 1.48503*obj.omegaVec(i) - ...
                        0.164423*(obj.omegaVec(i)^2) + 0.016666*(obj.omegaVec(i)^3);
                end
            end
            
            % ENERGY INTERACTION PARAMETER (PSI)
            obj.a = csvread(strcat(relDir,obj.solverInputDir,'a.csv'));
            
            % CALCULATE LENNARD-JONES PARAMETERS
            % USING TEE, GOTOH, STEWART
            obj.epsVec = (0.7915 + 0.1693 * obj.omegaVec) .* obj.TcVec;
            obj.SigmaVec = (1e-10) * (2.3551 - 0.0874 * obj.omegaVec).*((101325*obj.TcVec./(obj.PcVec)).^(1/3));
           
            % COMPUTE STANDARD HEAT OF VAPORIZATION
            obj.L = 6.829 + ((obj.funcMat*obj.propMat(5,:)'))';
            obj.L = 1e3 * obj.L ./ obj.MW;
            
        end
        
        function r = c_l(obj,T)
            % COMPUTE c_lVec (TODO : Change to non-standard sp. ht. capacity)
            % UNITS : (J/Kg-K)
            theta = (T-298)/700;
            if (T > 700)
                disp('Drop temperature exceeded 700 K');
            end
            r = ((obj.funcMat*obj.propMat(7,:)')' - 19.7779) + ((obj.funcMat*obj.propMat(8,:)')' + 22.5981)*theta + ...
                ((obj.funcMat*obj.propMat(9,:)')' - 10.7983)*(theta^2);
            r = r./obj.MW;
        end
        
        function r = D(obj,p,Tin)
            % READ DIFFUSION PROPERTIES (FUNCTIONAL GROUPS)
            % Currently based on Wilke-Lee
            % TODO : Change to functional group
            eps_air = 78.6;
            Sigma_air = 3.711e-10;
            
            evaVec = obj.epsVec*eps_air;
            Tstar = Tin./evaVec;
            omegaD = @(T) 1.06036./(T.^0.15610) + 0.193./exp(0.47635*T) + 1.03587./exp(1.52996*T) + 1.76474./exp(3.89411*T);
            MW_air = 28.97e-3;
            MWa = 2e3*obj.MW*MW_air./(obj.MW + MW_air);
            r = 101325*(3.03 - 0.98./(MWa.^0.5))*(1e-27)*(Tin^1.5)./(p*(MWa.^0.5).*((obj.SigmaVec + Sigma_air).^2).*omegaD(Tstar));
        end
        
        function r = specVol(obj,T)
            % RACKETT EQUATION
            % COMPUTE LIQUID SPECIFIC VOLUME (m3/mol) AND LIQUID DENSITY (Kg/m3)
            phiVec = zeros(1,length(obj.TcVec));
            for i = 1:length(obj.TcVec)
                if (T > obj.TcVec(i))
                    phiVec(i) = - ((1 - (298./obj.TcVec(i)) ).^(2/7));
                else
                    phiVec(i) = ((1 - (T./obj.TcVec(i)) ).^(2/7)) - ((1 - (298./obj.TcVec(i)) ).^(2/7));
                end
            end
            ZcVec = (0.29056 - 0.08775*obj.omegaVec);
            r = (1e-3)*(-0.00435 + ((obj.funcMat*obj.propMat(6,:)'))');
            r = r.*(ZcVec.^phiVec);
        end
        
        % COMPUTE PSatVec
        % LEE-KESLER SATURATION VAPOR PRESSURE
        function r = Psat(obj,T)
            f0 = @(Tr) 5.92714 - 6.09648./Tr - 1.28862*log(Tr) + 0.169347*(Tr.^6);
            f1 = @(Tr) 15.2518 - 15.6875./Tr - 13.4721*log(Tr) + 0.43577*(Tr.^6);
            r = obj.PcVec.*exp( f0(T./obj.TcVec) + obj.omegaVec.*f1(T./obj.TcVec) );
        end
        
        function r = Tb(obj,p)
            % SWITCH TO BOILING VECTOR FOR ATMOSPHERIC PRESSURE
            if (p >= 1e5 && p <= 1.1e5)
                disp('Using boiling point correlation');
                r = obj.TbVec;
            else
                
                f0 = @(Tr) 5.92714 - 6.09648./Tr - 1.28862*log(Tr) + 0.169347*(Tr.^6);
                f1 = @(Tr) 15.2518 - 15.6875./Tr - 13.4721*log(Tr) + 0.43577*(Tr.^6);
                r = zeros(1,length(obj.TcVec));
                T = sym('T');
                for i = 1:length(obj.TcVec)
                    r(i) = double(vpasolve(obj.PcVec(i).*exp(f0(T./obj.TcVec(i)) + obj.omegaVec(i).*...
                        f1(T./obj.TcVec(i)) ) == p,T,[250 800]));
                end
            end
        end

        function [ ret ] = activity(obj,xVec,T)
                        
            % R TABLE
            R = obj.propMat(11,:);
            assert(length(R) == obj.numGroups);
            
            % Q TABLE
            Q = obj.propMat(12,:);
            assert(length(Q) == obj.numGroups);
            
            % ENERGY INTERACTION PARAMETER
            Psi = exp(-obj.a/T);
            
            % COORDINATION NUMBER
            z = 10;
            
            % r VECTOR
            r = (obj.funcMat*R')';
            
            % q VECTOR
            q = (obj.funcMat*Q')';
            
            % L VECTOR
            LVec = (z/2)*(r - q) - (r - 1);
            
            % THETA VECTOR
            thetaVec = xVec.*q/dot(xVec,q);
            
            % PHI VECTOR
            phiVec = xVec.*r/dot(xVec,r);
            
            % COMBINATORIAL ACTIVITY
            gammaCVec = exp(log(phiVec./xVec) + (z/2)*q.*log(thetaVec./phiVec) + LVec ...
                - (phiVec./xVec)*dot(xVec,LVec));
            
            % SET ACTIVITY OF COMPLETELY VAPORIZED SPECIES TO 1
            gammaCVec(isnan(gammaCVec)) = 1;
            
            %------------------------------------------
            % RESIDUAL ACTIVITY ROUTINE
            %------------------------------------------
            
            % CALCULATE GROUP MOLE FRACTION
            XVec = xVec*obj.funcMat;
            XVec = XVec/sum(XVec);
            
            % CALCULATE Theta VECTOR
            ThetaVec = Q.*XVec/dot(Q,XVec);
            
            % CALCULATE FUNCTIONAL GROUP AND ISOLATED GAMMAS
            
            % AVOID DIVISION-BY-ZERO IN LAST TERM
            divRes = ThetaVec./(ThetaVec*Psi);
            divRes(isnan(divRes)) = 0;
            
            GammaVec = Q.*(1 - log(ThetaVec*Psi) - divRes*Psi' );
            GammaIVec = zeros(obj.num_comp,obj.numGroups);
            
            % REFER TO IJHMT PAPER FOR DETAILS
            for i = 1:obj.num_comp
                XIVec = obj.funcMat(i,:)/sum(obj.funcMat(i,:));
                ThetaIVec = Q.*XIVec/dot(Q,XIVec);
                
                divRes = ThetaIVec./(ThetaIVec*Psi);
                divRes(isnan(divRes)) = 0;
                
                GammaIVec(i,:) = Q.*(1 - log(ThetaIVec*Psi) - divRes*Psi' );
                
            end
            
            % RESIDUAL ACTIVITY
            gammaRVec = zeros(1,obj.num_comp);
            for i = 1:obj.num_comp
                gammaRVec(i) = exp( dot( obj.funcMat(i,:) , GammaVec - GammaIVec(i,:)  ) );
            end
            
            % TOTAL ACTIVITY
            ret = gammaCVec.*gammaRVec;
            
        end
        
    end
    
end
