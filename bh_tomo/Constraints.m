classdef Constraints < matlab.mixin.Copyable
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        slowness
        attenuation
        ind_reservoir
    end
    
    methods
        function obj = Constraints()
            obj.slowness.data = [];
            obj.slowness.data_xi = [];
            obj.slowness.data_theta = [];
            obj.attenuation.data = [];
            obj.ind_reservoir = [];
        end
    end
    
end

