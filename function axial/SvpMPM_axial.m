function [S, G, T] = SvpMPM_axial(rp, rv, h)

%1D axisymmetric material point basis functions (standard MPM)
%--------------------------------------------------------------------------
% Author: [Your Name]
% Date:   [Date]
% Description:
% Standard MPM (single point integration) for axisymmetric problems.
% 
% Unlike planar MPM which outputs [Svp, dSvp], axisymmetric requires 
% three functions:
%   S = N_i(rp)                 - standard shape function (scalar mapping)
%   G = dN_i/dr (rp)           - standard gradient (radial strain)
%   T = N_i(rp)/rp             - hoop coupling (axisymmetric specific,
%                                used for hoop strain eps_theta = u_r/r
%                                and hoop stress contribution to Fint)
%
% Truncation handling:
%   When particle is on/near symmetry axis (rp -> 0), T = N/rp becomes
%   singular. This function implements axis truncation by regularizing
%   T to zero when rp < tol, relying on v_r=0 boundary condition at axis.
%   Physically, at r=0, eps_theta = eps_rr, but numerical treatment uses
%   the G value at axis for consistency (if needed, modify else clause).
%--------------------------------------------------------------------------
% [S, G, T] = SVPMPMAXISYM(rp, rv, h)
%--------------------------------------------------------------------------
% Input(s):
% rp    - particle radial position (can be 0)
% rv    - grid node radial position  
% h     - grid spacing (element length)
%--------------------------------------------------------------------------
% Output(s);
% S     - shape function value at particle
% G     - shape function gradient at particle
% T     - hoop coupling term S/rp (regularized at axis)
%--------------------------------------------------------------------------
% See also:
%   SvpGIMPaxisymTrunc (for GIMP version with domain truncation)
%--------------------------------------------------------------------------

delta = rp - rv;
tol   = 1e-12;  % truncation tolerance for symmetry axis

%==========================================================================
% Standard linear MPM shape functions (identical to planar MPM)
%==========================================================================
if -h < delta && delta <= 0                                                 % Particle in "left" element (rv is right node of element)
  S = 1 + delta/h;                                                          % N_i = 1 + (rp-rv)/h
  G = 1/h;                                                                  % dN_i/dr = 1/h
  
elseif 0 < delta && delta <= h                                              % Particle in "right" element (rv is left node of element)
  S = 1 - delta/h;                                                          % N_i = 1 - (rp-rv)/h  
  G = -1/h;                                                                 % dN_i/dr = -1/h
  
else                                                                        % Particle outside node influence domain
  S = 0; 
  G = 0;
end

%==========================================================================
% Axisymmetric hoop term T = S/rp with axis truncation
%==========================================================================
if rp > tol
    % Standard case: particle away from axis
    T = S / rp;
    
else
    % TRUNCATED CASE: Particle on or very near symmetry axis (rp -> 0)
    % 
    % Physical background:
    % At r=0, eps_theta = u_r/r is indeterminate (0/0), but by symmetry 
    % eps_theta = eps_rr = du_r/dr. Thus T should theoretically equal G 
    % for the axis node (rv=0) and 0 for others.
    %
    % Numerical treatment: Set T=0 for all nodes when rp<tol. This is safe
    % because:
    % 1. Boundary condition enforces v_r=0 at r=0, so u_r=0 at axis
    % 2. For off-axis nodes, S=0 when rp=0, so S/rp is technically 0/0 but
    %    physically the influence is zero
    % 3. For axis node (rv=0), if needed, use G instead of T for hoop strain
    
    T = 0;
    
    % Alternative (if you need strict consistency with eps_theta = eps_rr):
    % if abs(rv) < tol && S > 0.5  % essentially the axis node
    %     T = G;  % enforce eps_theta = eps_rr at axis
    % else
    %     T = 0;
    % end
end

end