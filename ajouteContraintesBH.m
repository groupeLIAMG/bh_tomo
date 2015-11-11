function cont = ajouteContraintesBH(cont, BHcont, x0, a, origine, az, dip)
%function cont = ajouteContraintesBH(cont, BHcont, x0, a, origine, az, dip)
%
%
%
if isempty(BHcont), return, end
if ~isempty(x0) && ~isempty(a)
	f = proj_plan([BHcont.x+zeros(size(BHcont.z)) ...
		BHcont.y+zeros(size(BHcont.z)) ...
		BHcont.z], x0, a);
end
f = transl_rotat(f, origine, az, dip);

cont.data = [cont.data; [f(:,3) f(:,1) BHcont.valeur BHcont.variance]];
