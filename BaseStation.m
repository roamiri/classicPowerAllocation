classdef BaseStation
   properties
      X
      Y
      P
      eta = 1e12
      gBS_MUE
   end
   methods
      function obj = BaseStation(xPos, yPos, transmitPower)
      if nargin > 0
         if isnumeric(xPos)
            obj.X = xPos;
         else
            error(' X Value must be numeric')
         end
         if isnumeric(yPos)
            obj.Y = yPos;
         else
            error(' Y Value must be numeric')
         end
         if isnumeric(transmitPower)
            %obj.P = 10^((transmitPower-30)/10);
            obj.P = transmitPower;
         else
            error(' P Value must be numeric')
         end
      end
      end
      function variable = getLagrangeVar(obj)
          variable = obj.eta(size(obj.eta,2));
      end
      
      function obj = updateLagrangeVar(obj,eta)
          obj.eta = [obj.eta eta];
      end
      function obj = update_gBS_MUE(obj,value)
          obj.gBS_MUE = value;
      end
   end
end