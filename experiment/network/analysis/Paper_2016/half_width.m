function [hw, pnts ] = half_width(vt, rmp)

% calculate half width for single pulse:
% Trace must have single pulse!
% output: hw half width in vt units
% output: pnts start and end points  in vt units.
% TODO: plotting, to enable one look validation.

N = length(vt);
[m,I] = max(vt);
vtd = vt-rmp;
M = (m-rmp);
possibleCrossings = ( (vtd > M*0.4) & (vtd < M*0.6));

[S,E] = regexp( sprintf('%i',possibleCrossings), '1+', 'match' );

s = cellfun(@(x) size(x,2),S);

[~,Is] = min( abs(vtd(E(1):E(1)+s(1)) - M/2) );
[~,Ie] = min( abs(vtd(E(2):E(2)+s(2)) - M/2) );

pnts(1) = E(1) + Is;
pnts(2) = E(2) + Ie;

hw = pnts(2) - pnts(1);