function [vertices] = rectgen(x,y,z)

% Function to produce the vertices of a rectangular prism, with side
% lengths x, y, and z. x, y, and z are assumed to be aligned with the
% normal coordinate axes, and the shape is assumed to be centered on the
% origin

vertices =	[
            
            -0.5,	-0.5,	-0.5
            0.5,    -0.5,   -0.5
            -0.5,	-0.5,	0.5
            0.5,    -0.5,   0.5
            -0.5,	0.5,	-0.5
            0.5,    0.5,	-0.5
            -0.5,	0.5,	0.5
            0.5,    0.5,	0.5
            
            ];

vertsizes = [8,3];

sizemult = repmat([x,y,z],[vertsizes(1),1]);

vertices =	vertices.*sizemult;

end