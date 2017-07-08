classdef UE
   properties
      X
      Y
      SINR
      C
      C_profile = []
   end
   methods
      function obj = UE(xPos, yPos)
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
      end
      end
      
      function obj = setCapacity(obj,c)
        obj.C = c;
        obj.C_profile = [obj.C_profile c];
      end
   end
end