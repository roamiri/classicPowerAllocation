classdef FemtoStation
   properties
      X
      Y
      P
      dBS
      dMUE
      dFUE
      FUEX
      FUEY
      M  % distance with MUE
      B  % distance with BS
      dM1 = 15; dM2 = 50; dM3 = 125; 
      dB1 = 50; dB2 = 150; dB3 = 400;
      state = zeros(1,2)
      powerProfile = []
      C_FUE
      C_profile = []
      lambda = 10.0 %lagrange var for Pmax
      gamma = 1.10  %lagrange var for Pmin
      mu = 100     %lagrange var for required rate
      gmf
      gf
      If
      R_FUE % FUE required rate
   end
   methods
      function obj = FemtoStation(xPos, yPos, BS, MUE, dFUE)
        obj.X = xPos;
        obj.Y = yPos;
        obj.dBS = sqrt((xPos-BS.X)^2 + (yPos-BS.Y)^2);
        obj.dMUE = nearest_MUE(xPos, yPos, MUE);% sqrt((xPos-MUE.X)^2 + (yPos-MUE.Y)^2); %distance to nearest MUE
        obj.dFUE = dFUE;
        obj.FUEX = xPos;
        obj.FUEY = yPos+dFUE;
      end
      
      function [lambda, gamma, mu] = getLagrangeVars(obj)
        ss = size(obj.mu,2);
        mu = obj.mu(ss);
        lambda = obj.lambda(ss);
        gamma = obj.gamma(ss);
      end
      
      function obj = updateLagrangeVars(obj,lambda, gamma, mu)
          obj.lambda = [obj.lambda lambda];
          obj.gamma = [obj.gamma gamma];
          obj.mu = [obj.mu mu];
      end
      
      function obj = setGMF(obj, gmf)
          obj.gmf = gmf;
      end
      function obj = setGF(obj, gf)
          obj.gf = gf;
      end
      function obj = setInterf(obj, If)
          obj.If = If;
      end
      function obj = setR_FUE(obj,R)
          obj.R_FUE = R;
      end
      function obj = setPower(obj,power)
%           obj.P = 10^((power-30)/10);
            obj.P = power;
            obj.powerProfile = [obj.powerProfile power];
      end
      
      function obj = setCapacity(obj,c)
        obj.C_FUE = c;
        obj.C_profile = [obj.C_profile c];
      end
      function obj = getDistanceStatus(obj)
          if(obj.dMUE <= obj.dM1 )
              obj.state(1) = 0;
          elseif(obj.dMUE <= obj.dM2 )
              obj.state(1) = 1;
          elseif(obj.dMUE <= obj.dM3 )
              obj.state(1) = 2;
          else
              obj.state(1) = 3;
          end
          
          if(obj.dBS <= obj.dB1 )
              obj.state(2) = 0;
          elseif(obj.dBS <= obj.dB2 )
              obj.state(2) = 1;
          elseif(obj.dBS <= obj.dB3 )
              obj.state(2) = 2;
          else
              obj.state(2) = 3;
          end
      end
   end
end