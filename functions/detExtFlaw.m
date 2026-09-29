function Q_fext = detExtFlaw(nodes,mpData)


nmp  = size(mpData,2);                                                      % number of material points & dimensions 
Q_fext = zeros(nodes,1);                                                    % zero the external force vector
for mp = 1:nmp
   nIN = mpData(mp).nIN;                                                    % nodes associated with MP
   Svp = mpData(mp).Svp;                                                    % basis functions
   fp  = mpData(mp).Flow_P_current*Svp;                                             % material point body & point nodal forces
   Q_fext(nIN) = Q_fext(nIN) + fp';                                          % combine into external force vector
end
end