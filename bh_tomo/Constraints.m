classdef Constraints < matlab.mixin.Copyable
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        slowness
        attenuation
        ind_reservoir
    end
    
    methods
        function obj = Constraints(varargin)
            if nargin == 0
                obj.slowness.data = [];
                obj.slowness.data_xi = [];
                obj.slowness.data_theta = [];
                obj.attenuation.data = [];
                obj.attenuation.data_xi = [];
                obj.attenuation.data_theta = [];
                obj.ind_reservoir = [];
            else
                if ~isstruct(varargin{1})
                    error('Invalid input')
                end
                s = varargin{1};
                obj.slowness.data = s.slowness;
                obj.attenuation.data = s.attenuation;
            end
        end
    end
    
end

