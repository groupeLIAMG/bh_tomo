classdef CovarianceModels < int8
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
        Cubic (1)
        Spherical (2)
        Gaussian (3)
        Exponential (4)
        Linear (5)
        Thin_Plate (6)
        Gravimetric (7)
        Magnetic (8)
        Hole_Effect_Sine (9)
        Hole_Effect_Cosine (10)
        Nugget (11)
    end
    methods (Static=true)
        function c = buildCov(type,r,a,s)
            switch(type)
                case CovarianceModels.Cubic
                    c = CovarianceCubic(r,a,s);
                case CovarianceModels.Spherical
                    c = CovarianceSpherical(r,a,s);
                case CovarianceModels.Gaussian
                    c = CovarianceGaussian(r,a,s);
                case CovarianceModels.Exponential
                    c = CovarianceExponential(r,a,s);
                case CovarianceModels.Linear
                    c = CovarianceLinear(r,a,s);
                case CovarianceModels.Thin_Plate
                    c = CovarianceThinPlate(r,a,s);
                case CovarianceModels.Gravimetric
                    c = CovarianceGravimetric(r,a,s);
                case CovarianceModels.Magnetic
                    c = CovarianceMagnetic(r,a,s);
                case CovarianceModels.Hole_Effect_Sine
                    c = CovarianceHoleEffectSine(r,a,s);
                case CovarianceModels.Hole_Effect_Cosine
                    c = CovarianceHoleEffectCosine(r,a,s);
                case CovarianceModels.Nugget
                    c = CovarianceNugget(s);
            end
                    
        end
        function c = getDefault2D()
            c = CovarianceSpherical([4 4],0,1);
        end
        function c = getDefault3D()
            c = CovarianceSpherical([4 4 4],[0 0 0],1);
        end
    end
end
