classdef test_class
    properties
        a 
    end
    methods
        function obj = test_class(a, b)
            syms f z
            f = diff(z^2*a*b,z);
            obj.a =f;
        end
    end
end