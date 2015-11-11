function [covara,dif]=ajuster_Cov(options,covar,covar_code,L,...
	gridx,gridz,x,iktt,afi,lclas,c0,ax1,ax2)
%  fonction pour ajuster simultanement un modele dans plusieurs directions


% verifier si l'on est avec un cas 2D anisotrope ou 3D anisotrope
%[n,p] = size(covar.model);

id{1} = covar_code.model==0;
p1 = sum(sum(id{1}));
x0libre = covar.model(id{1});
x0libre = x0libre(:);  % on s'assure qu'on a un vecteur colonne

id{2} = covar_code.c==0;
p2 = sum(sum(id{2}));
x0libre = [x0libre; covar.c(id{2})];

id{3} = covar_code.nugget_t==0;
p3 = sum(id{3});
x0libre = [x0libre; covar.nugget_t(id{3})];

id{4} = covar_code.nugget_l==0;
p4 = sum(id{4});
x0libre = [x0libre; covar.nugget_l(id{4})];

if covar.aniso > 0
    id{5} = covar_code.model_xi==0;
    p5 = sum(sum(id{5}));
    tmp = covar.model_xi(id{5});
    x0libre = [x0libre; tmp(:)];
    
    id{6} = covar_code.c_xi==0;
    p6 = sum(sum(id{6}));
    x0libre = [x0libre; covar.c_xi(id{6})];
    
    id{7} = covar_code.nugget_xi==0;
    p7 = sum(sum(id{7}));
    x0libre = [x0libre; covar.nugget_xi(id{7})];
end
if covar.aniso > 1
    id{8} = covar_code.model_th==0;
    p8 = sum(sum(id{8}));
    tmp = covar.model_th(id{8});
    x0libre = [x0libre; tmp(:)];
    
    id{9} = covar_code.c_th==0;  %%YH  .c_th_code  --> .c_th
    p9 = sum(sum(id{9}));
    x0libre = [x0libre; covar.c_th(id{9})];
    
    id{10} = covar_code.nugget_th==0;
    p10 = sum(sum(id{10}));
    x0libre = [x0libre; covar.nugget_th(id{10})];
end


if ~isempty(x0libre)
	[t1,dif,flag]=fminsearch('fitModeliKss',x0libre,options,covar,...
		id,L,gridx,gridz,x,iktt,afi,lclas,c0,ax1,ax2);
end

if p1>0;
	t11 = t1(1:p1);
	covar.model(id{1}) = t11;
end
if p2>0
	t21 = t1(p1+1:p1+p2);
	covar.c(id{2}) = t21;
end
if p3>0
	t31 = t1(p1+p2+1:p1+p2+p3);
	covar.nugget_t(id{3}) = t31;
end
if p4>0
	t41 = t1(p1+p2+p3+1:p1+p2+p3+p4);
	covar.nugget_l(id{4}) = t41;
end
if covar.aniso > 0
    if p5>0
        t51 = t1(p1+p2+p3+p4+1:p1+p2+p3+p4+p5);
        covar.model_xi(id{5}) = t51;
    end
    if p6>0
        t61 = t1(p1+p2+p3+p4+p5+1:p1+p2+p3+p4+p5+p6);
        covar.c_xi(id{6}) = t61;
    end
    if p7>0
        t71 = t1(p1+p2+p3+p4+p5+p6+1:p1+p2+p3+p4+p5+p6+p7);   %%YH  p5+1 --> p5+p6+1
        covar.nugget_xi(id{7}) = t71;
    end
end
if covar.aniso > 1
    if p8>0
        t81 = t1(p1+p2+p3+p4+p5+p6+p7+1:p1+p2+p3+p4+p5+p6+p7+p8);  %%YH  p6+1 --> p6+p7+1
        covar.model_th(id{8}) = t81;
    end
    if p9>0
        t91 = t1(p1+p2+p3+p4+p5+p6+p7+p8+1:p1+p2+p3+p4+p5+p6+p7+p8+p9);  %%YH  p7+1 --> p7+p8+1
        covar.c_th(id{9}) = t91;
    end
    if p10>0
        t101 = t1(p1+p2+p3+p4+p5+p6+p7+p8+p9+1:p1+p2+p3+p4+p5+p6+p7+p8+p9+p10);  %%YH
        covar.nugget_th(id{10}) = t101;
    end
end
covara = covar;
