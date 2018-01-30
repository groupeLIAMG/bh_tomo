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
                if isstruct(s.slowness)
                    obj.slowness.data = s.slowness.data;
                    obj.slowness.data_xi = [];
                    obj.slowness.data_theta = [];
                else
                    obj.slowness.data = s.slowness;
                end
                if isstruct(s.attenuation)
                    obj.attenuation.data = s.attenuation.data;
                    obj.attenuation.data_xi = [];
                    obj.attenuation.data_theta = [];
                else
                    obj.attenuation.data = s.attenuation;
                end
                if isfield(s, 'ind_reservoir')
                    obj.ind_reservoir = s.ind_reservoir;
                else
                    obj.ind_reservoir = [];
                end
            end
        end
    end
    
end

