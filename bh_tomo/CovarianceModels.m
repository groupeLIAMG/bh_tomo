classdef CovarianceModels
    %COVARIANCEMODELS Admissible covariance models
    
    % Models taken from Covardm, which are:
    %'h==0                                     '; nugget
    %'exp(-h)                                  '; exponential
    %'exp(-(h).^2)                             '; gaussian
    %'1-(1.5*min(h,1)/1-.5*(min(h,1)/1).^3)    '; spherical
    %'1-h                                      '; linear
    %'1-3*min(h,1).^2+2*min(h,1).^3            '; modele cubique
    %'(h.^2).*log(max(h,eps))                  '; spline plaque mince
    %'(h.^2+1).^(-0.5)                         '; modele gravimetrique (Cauchy avec b=0.5)
    %'(h.^2+1).^(-1.5)                         '; modele magnetique (Cauchy avec b=1.5) 
    %'sin(max(eps,h*2*pi))./max(eps,h*2*pi)    '; effet de trou sinusoidal
    %'cos(h*2*pi)                              '; effet de trou cosinusoidal
    
    enumeration
        nugget, exponential, gaussian, spherical, linear, cubic, ...
        thin_plate, gravimetric, magnetic, hole_effect_sine, ...
        hole_effect_cosine
    end
end
