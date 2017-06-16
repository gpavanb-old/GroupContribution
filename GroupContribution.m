classdef GroupContribution
    properties
        % DEFINE RELATIVE PATHS
        solverInputDir = 'data/Solver_Input/'
        compDescDir = 'data/Compound_Descriptions/'
        initDataDir = 'data/Init_Data/'

        propMat
        funcMat
        initMat
        
        MW
        L
                
        PcVec
        TcVec
        VcVec
        omegaVec
        kappaVec
        epsVec
        SigmaVec
        a
        
        x_l_init
        
    end
    
    methods
        function obj = GroupContribution(funcName,initName)
            relDir = 'include/GroupContribution/data';
            obj.propMat = xlsread(strcat(relDir,obj.solverInputDir,'gani_prop_table.xlsx'));
            obj.funcMat = xlsread(strcat(relDir,obj.compDescDir,funcName,'.xlsx'));
            
            [obj.initMat,~] = xlsread(strcat(relDir,obj.initDataDir,initName,'.xlsx'));
            obj.x_l_init = (obj.initMat(:,1))';
            obj.x_l_init = obj.x_l_init/sum(obj.x_l_init);
            
            obj.MW = (1e-3 * obj.funcMat*obj.propMat(10,:)')';
            
            % CALCULATE CRITICAL PROPERTIES
            TcCoVec = obj.propMat(1,:);
            PcCoVec = obj.propMat(2,:);
            VcCoVec = obj.propMat(3,:);
            omegaCoVec = obj.propMat(4,:);
            
            % CALCULATE CRITICAL VOLUME VECTOR
            obj.VcVec = 1e-3*(-0.00435 + ( (obj.funcMat*VcCoVec')' ));
            
            % CALCULATE CRITICAL TEMPERATURE VECTOR
            obj.TcVec = 181.128*log( (obj.funcMat*TcCoVec')' );
            
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
            a = csvread(strcat(relDir,obj.solverInputDir,'a.csv'));
            
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
        
    end
    
end
