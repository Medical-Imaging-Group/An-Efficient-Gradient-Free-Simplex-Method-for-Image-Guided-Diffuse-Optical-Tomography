function [ omega]= ObjectiveFunction_Opt_HardPrior(mu,fwd_mesh, anom,R,RL,qvec)

for i = 1:length(mu)  
    fwd_mesh.mua(R(1:RL(i),i)) =  (mu(i));   
end
fwd_mesh.kappa = (1./(3.*(fwd_mesh.mus+fwd_mesh.mua)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frequency =0;

% modulation frequency
omega = 2*pi*frequency*1e6;
mesh = fwd_mesh;
% Create FEM matricex
if mesh.dimension == 2
  [i,j,s] = gen_matrices_2d(mesh.nodes(:,1:2),...
			    sort(mesh.elements')', ...
			    mesh.bndvtx,...
			    mesh.mua,...
			    mesh.kappa,...
			    mesh.ksi,...
			    mesh.c,...
			    omega);

end
            
%junk = length(find(i==0));
%MASS = sparse(i(1:end-junk),j(1:end-junk),s(1:end-junk));
%%next two lines modified by subha, 6/11/07 to make it more memory efficient
nz_i = nonzeros(i);
MASS = sparse(i(1:length(nz_i)),j(1:length(nz_i)),s(1:length(nz_i)));
clear junk i j s omega nz_i  
% Catch zero frequency (CW) here
if frequency == 0
  MASS = real(MASS);
  qvec = real(qvec);
end
% Calculate field for all sources
[data.phi,mesh.R]=get_field(MASS,mesh,qvec);

% Calculate boundary data
[data.complex]=get_boundary_data(mesh,data.phi);
data.link = mesh.link;

% Map complex data to amplitude and phase
data.amplitude = abs(data.complex);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ref = log(data.amplitude);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega = norm(ref-anom, 2)^2;

end

