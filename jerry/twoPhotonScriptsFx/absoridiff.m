function degrees = absoridiff(alpha,beta);
%ABSORIDIFF the angle between two ORIENTATIONS (NOT DIRECTIONS)
% DEGREES = ABSORIDIFF(ALPHA,BETA) where
% DEGREES is the absolute difference in degrees of alpha and beta
% ALPHA and BETA are orientations in degrees with a range of (0, 179)

nordeg = mod(alpha-beta,180);
degrees = min(180-nordeg, nordeg);

end