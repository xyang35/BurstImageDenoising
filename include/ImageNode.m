classdef ImageNode
   properties
      range
      pts
      inds
   end
   methods
       function obj = ImageNode(r, p, i)
           if nargin > 0
               obj.range = r;
               obj.pts = p;
               obj.inds = i;
           end
       end
   end
end