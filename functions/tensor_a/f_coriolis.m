function f=f_coriolis(coriolis,y_i)
% Compute the Coriolis frequency
%

if nargin > 1
    y_origin = coriolis.y_origin;
    y_i = y_i + y_origin;
else
    if ~ strcmp(coriolis.f_model,'f_plane')
        error('Location information is needed to compute the Coriolis frequency');
    end
end

switch coriolis.f_model
    case 'f_plane'
        f=coriolis.f0;
    case 'beta_plane'
        f= coriolis.beta * y_i + coriolis.f0;
    case 'general'
        f= 2* coriolis.OMEGA * sin( 2*pi/coriolis.earth_radius * y_i );
    otherwise
        error('incorrect coriolis model');
end
end