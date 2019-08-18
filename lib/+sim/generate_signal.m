function [ x,varargin] = generate_signal( source,angle,r,distant,beta,varargin )
%----------------------------------------------------------------------
%                 generate micarray simulation signal 
%function
%
% Usage: (1) x = generate_signal( source,[90,0],0.032,1,0.2 ) 
%
% Inputs:
%   source         source signal
%   angle          incident angle
%   r              radius of array
%   beta           reverberation parameter
%   varargin
%       varargin{1}  scale of room-impulse-response
%
% Outputs:
%   x            multi-channel signal
%  varargin{1}   room-impulse-response
% dependencies:
%   rir-generator
%
%   Authors: Wang wei
%   $Revision: 0.0 $  $Date:  $
%
% ---------------------------------------------------------------------

if nargin<2
    fprintf('Usage: [x]=generate_signal(source,angle) \n');
    return;
end;
if nargin<3
    r = 0.032;
end;
if nargin<4
    distant = 1.0;
end;
if nargin<5
    beta = 0.2;
end;
if nargin==6
    scale = varargin{1};
else
    scale = 10;
end;


angle = [angle(1) angle(2)]/180*pi;              % source direction [0,180]
[x1,y1,z1]=sph2cart(angle(1),angle(2),distant);  % source position 1
source_pos = [x1,y1,z1];                         % Source position [x y z] (m)

N = 4;                                           % number of sensor

h = sim.RIR_generator_URA( source_pos,beta,r);
h = h*scale;
x = zeros(length(source)+size(h,2)-1,N);
for i=1:N
    x(:,i) = conv(source,h(i,:));
end

if nargout>1
    varargin{1} = h;
end

end

