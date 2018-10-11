classdef FData
    % This class includes the FFT Data blocks
    properties
        Data % is a n x m x ep
        Condition
    end
    properties (Dependent)
        Epochs
    end
    
    methods 
        function obj = FData(Data,Condition)
            if exist('Data','var')
                obj.Data = Data;
            else
                obj.Data = [];
            end
            if exist('Condition','var')
                obj.Condition = Condition;
            else
                obj.Condition =[];
            end
        end
        
        function value = get.Epochs(obj)
            value = size(obj.Data,3);
        end
    end
end