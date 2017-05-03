%% Test

close all
clear
clc

% Box side lengths

x = 2;
y = 3;
z = 4;

%Extra space in plot (percentage of max axis size)

extraspace = .5;

% Vector representing graph axes to show

axnew = max([x,y,z]/2)*(1 + extraspace)*reshape([-ones(1,3);ones(1,3)],6,1);

% Number of interpolation steps for each scheme (linear, Euler, SLERP)

interpsize = 1001;

% VARIABLES FOR PLOTS AND THINGS BELOW

% Maximum and minimum box transperencies

maxtrans = .8;
mintrans = .2;

% Transparency transition time (there is an animation where a box becomes
% transparent)

transtime = 1;

% Frames per second for all movies; 120FPS is used because the gamers like
% it.

FPS = 120;

% Establishing a timer for animations; ONLY active if video is NOT being
% saved (saving the video takes longer, so the timer is not useful).
% The timer isn't very good at enforcing the FPS, so generally, the timer
% is not very accurate. But it's better than letting it go at LIGHTNING
% speed.

t = timer('StartDelay',1/FPS,'TimerFcn','0;');

% Number of time steps for the transparency animation

tranits = ceil(transtime*FPS)+1;

% Transparency vector: each point corresponds to changing (diminishing)
% transparency at each time step (of 1/FPS seconds).

transvec = linspace(maxtrans,mintrans,tranits);

% Extra time to view the glory of the fully rotated model

etime = 2;

% Extra time number of total frames

etframes = etime*FPS;

% Variable to allow or disallow saving of video; saving the video takes a
% TON of time, so for demonstration of this code, set it to 0.

savevid = 0;

% Order of vertices to produce 6 faces. If you plot these in sequence, it
% "spirals" through them!

faceorder =	[
            
            1,  3,  4,  2;
            3,  4,  8,  7;
            2,  6,  8,  4;
            1,  2,  6,  5;
            1,  5,  7,  3;
            5,  7,  8,  6;
            
            ];

% Generating vertices for rectangle using a simple function

vertices = rectgen(x,y,z);

% Rectangle face colors, to make distinguishing easier.

col =	[
        
        1,      0,      0;
        1,      0.3,    1;
        0,      0       1;
        0.6,	0.6,	0.6;
        1,      1,      0;
        0,      1,      1;
        
        ];

Eulers =	[
            
            3,  2,	1;
            70,	50,	10;
            
            ];

tic

Eulerinterps = linspaceNDim(zeros(1,3),Eulers(2,:),interpsize)';

linEultime = toc;

% Finding a rotation matrix for the final rotation given; note that
% "spincalc" creates a matrix that is transposed (it gives rotation of
% WORLD, not of object), so we need to transpose it).

eval(sprintf('rotmat = SpinCalc(''EA%i%i%itoDCM'',Eulers(2,:))'';',Eulers(1,:)))
eval(sprintf('rotqt = SpinCalc(''EA%i%i%itoQ'',Eulers(2,:));',Eulers(1,:)))

% For whatever reason, the user that made SpinCalc has the quaternion
% ordered like q = [q2, q3, q4, q1] so I have to reorder it; the logic
% behind that move baffles me... I ALSO have to invert the angle, just like
% with the matrix, again because the reported angle is from the world
% frame, not the object frame.

rotqt = [rotqt(4),-rotqt(1:3)];

tic

rotang = acos(rotqt(1))*2;

slerptime = toc;

tic

linrots = linspaceNDim(eye(3),rotmat,interpsize);

lintime = toc;

linEuls = zeros([3,3,interpsize]);
linEuls(:,:,1) = eye(3);
linEuls(:,:,end) = rotmat;

slerprots = linEuls;

slerpquats =	(repmat([1,0,0,0],[interpsize,1]).*repmat(sin(linspace(1,0,interpsize)'*rotang),[1,4])+...
                repmat(rotqt,[interpsize,1]).*repmat(sin(linspace(0,1,interpsize)'*rotang),[1,4]))/sin(rotang);

% linear Euler calculation follows

% For performance reasons, I will simply loop through an Euler angle
% calcuation for a 321 sequence, but will then follow it with a generic one
% using sprintf (so I can control earlier in the code what Euler sequence I
% use: I can chooce 313, 321, etc.) eval(sprintf()) takes  a non-trivial
% amount of time to run, so it's not fair to time that operation for
% comparison.

tic

for ix = 2:interpsize-1
    
    linEuls(:,:,ix) = SpinCalc('EA321toDCM',Eulerinterps(ix,:))';
    
end

linEultime = linEultime+toc;

for ix = 2:interpsize-1
    
    eval(sprintf('linEuls(:,:,ix) = SpinCalc(''EA%i%i%itoDCM'',Eulerinterps(ix,:))'';',Eulers(1,:)))
    
end

verticesfinal = (rotmat*vertices')';
verticesfinalqt = quatrotate(rotqt,vertices);

verticeslin = repmat(vertices,[1,1,interpsize]);
verticeslin(:,:,1) = vertices;
verticeslin(:,:,end) = verticesfinal;
verticeslinEul = verticeslin;
verticesslerp = verticeslinEul;

for ix = 2:interpsize-1
    
    tic
    
    verticeslin(:,:,ix) = (linrots(:,:,ix)*vertices')';
    
    lintime = lintime+toc;
    
    tic
    
    verticeslinEul(:,:,ix) = (linEuls(:,:,ix)*vertices')';
    
    linEultime = linEultime+toc;
    
    tic
    
    verticesslerp(:,:,ix) = quatrotate(slerpquats(ix,:),vertices);
    
    slerptime = slerptime + toc;
    
end

% Print statements showing time to calculate matrices/quaternions AND apply
% them to the 8 vertices

fprintf('\nLinear interpolation time was %.4f seconds\n',lintime)
fprintf('\nEuler Angle Linear interpolation time was %.4f seconds\n',linEultime)
fprintf('\nSLERP time was %.4f seconds\n\n',slerptime)

%% Plotting, initial and final boxes

figme = figure(1);

clf

pause(10^-8)
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);

% This maximizes the window for us automatically

subplot(1,2,2)

hold on

% Building final prism

patch('Vertices',verticesfinal,'Faces',faceorder,'FaceColor','flat','FaceVertexCData',col,'FaceAlpha',maxtrans,'EdgeColor','k','LineWidth',2);

% Miscellaneous plot commands to make everything look nicer (we're using
% an isometric 3D view angle)

axis equal
grid on

axis(axnew)

% Isometric view angle
view(-45,35.264)

xlabel('Global X axis')
ylabel('Global Y axis')
zlabel('Global Z axis')

title('Final Prism')

subplot(1,2,1)

hold on

% Building initial prism

patch('Vertices',vertices,'Faces',faceorder,'FaceColor','flat','FaceVertexCData',col,'FaceAlpha',maxtrans,'EdgeColor','k','LineWidth',2);

% Miscellaneous plot commands to make everything look nicer (we're using
% an isometric 3D view angle)

axis equal
grid on

axis(axnew)

% Isometric view angle
view(-45,35.264)

xlabel('Global X axis')
ylabel('Global Y axis')
zlabel('Global Z axis')

title('Initial Prism')

figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')
set(findall(figureHandle,'type','axes'),'fontSize',14,'fontWeight','bold')

% Larger text

set(figureHandle,'PaperPositionMode','auto')
print('Prism Rotation Preview','-dpng','-r0')

if ~savevid
    
    pause(etime)
    
end

%% Plotting, linear interpoolation (awful, includes distortion)

figure(2)

clf

pause(10^-8)
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);

% This maximizes the window for us automatically

subplot(1,2,2)

hold on

% Building final prism

patch('Vertices',verticesfinal,'Faces',faceorder,'FaceColor','flat','FaceVertexCData',col,'FaceAlpha',maxtrans,'EdgeColor','k','LineWidth',2);

% Miscellaneous plot commands to make everything look nicer (we're using
% an isometric 3D view angle)

axis equal
grid on

axis(axnew)

% Isometric view angle
view(-45,35.264)

xlabel('Global X axis')
ylabel('Global Y axis')
zlabel('Global Z axis')

title('Final Prism, Linear Interpolation')

subplot(1,2,1)

hold on

% Building initial prism

initpris = patch('Vertices',vertices,'Faces',faceorder,'FaceColor','flat','FaceVertexCData',col,'FaceAlpha',maxtrans,'EdgeColor','k','LineWidth',2);

% Miscellaneous plot commands to make everything look nicer (we're using
% an isometric 3D view angle)

axis equal
grid on

axis(axnew)

% Isometric view angle
view(-45,35.264)

xlabel('Global X axis')
ylabel('Global Y axis')
zlabel('Global Z axis')

title('Initial Prism, Linear Interpolation')

figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')
set(findall(figureHandle,'type','axes'),'fontSize',14,'fontWeight','bold')

% Larger text

if savevid
    
    vid = VideoWriter('rotationlin.avi');
    vid.FrameRate = 120;
    open(vid);

end

for ix = 1:tranits
    
    if ~savevid
        
        start(t)
        
    end
    
    set(initpris,'FaceAlpha',transvec(ix))
    
    if savevid
        
        writeVideo(vid,getframe(gcf));
        
    else
        
        wait(t)
        
    end
    
end

for ix = 1:interpsize
    
    if ~savevid
        
        start(t)
        
    end
    
    rotme = patch('Vertices',verticeslin(:,:,ix),'Faces',faceorder,'FaceColor','flat','FaceVertexCData',col,'FaceAlpha',maxtrans,'EdgeColor','k','LineWidth',2);
    
    if savevid
        
        writeVideo(vid,getframe(gcf));
        
    else
        
        wait(t)
        
    end
    
    delete(rotme)
    
end

patch('Vertices',verticeslin(:,:,end),'Faces',faceorder,'FaceColor','flat','FaceVertexCData',col,'FaceAlpha',maxtrans,'EdgeColor','k','LineWidth',2);

if savevid
        
    finalframe = getframe(gcf);
    
    for ix = 1:etframes
        
        writeVideo(vid,finalframe);
        
    end
    
    close(vid)
    
else
    
    pause(etime)
    
end


%% Plotting, linear interpoolation of Euler angles (no distortion, but jumpy)

figure(3)

clf

pause(10^-8)
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);

% This maximizes the window for us automatically

subplot(1,2,2)

hold on

% Building final prism

patch('Vertices',verticesfinal,'Faces',faceorder,'FaceColor','flat','FaceVertexCData',col,'FaceAlpha',maxtrans,'EdgeColor','k','LineWidth',2);

% Miscellaneous plot commands to make everything look nicer (we're using
% an isometric 3D view angle)

axis equal
grid on

axis(axnew)

% Isometric view angle
view(-45,35.264)

xlabel('Global X axis')
ylabel('Global Y axis')
zlabel('Global Z axis')

title('Final Prism, Linear Euler Angle Interpolation')

subplot(1,2,1)

hold on

% Building initial prism

initpris = patch('Vertices',vertices,'Faces',faceorder,'FaceColor','flat','FaceVertexCData',col,'FaceAlpha',maxtrans,'EdgeColor','k','LineWidth',2);

% Miscellaneous plot commands to make everything look nicer (we're using
% an isometric 3D view angle)

axis equal
grid on

axis(axnew)

% Isometric view angle
view(-45,35.264)

xlabel('Global X axis')
ylabel('Global Y axis')
zlabel('Global Z axis')

title('Initial Prism, Linear Euler Angle Interpolation')

figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')
set(findall(figureHandle,'type','axes'),'fontSize',14,'fontWeight','bold')

% Larger text

if savevid
    
    vid = VideoWriter('rotationlin.avi');
    vid.FrameRate = 120;
    open(vid);

end

for ix = 1:tranits
    
    if ~savevid
        
        start(t)
        
    end
    
    set(initpris,'FaceAlpha',transvec(ix))
    
    if savevid
        
        writeVideo(vid,getframe(gcf));
        
    else
        
        wait(t)
        
    end
    
end

for ix = 1:interpsize
    
    if ~savevid
        
        start(t)
        
    end
    
    rotme = patch('Vertices',verticeslinEul(:,:,ix),'Faces',faceorder,'FaceColor','flat','FaceVertexCData',col,'FaceAlpha',maxtrans,'EdgeColor','k','LineWidth',2);
    
    if savevid
        
        writeVideo(vid,getframe(gcf));
        
    else
        
        wait(t)
        
    end
    
    delete(rotme)
    
end

patch('Vertices',verticeslinEul(:,:,end),'Faces',faceorder,'FaceColor','flat','FaceVertexCData',col,'FaceAlpha',maxtrans,'EdgeColor','k','LineWidth',2);

if savevid
        
    finalframe = getframe(gcf);
    
    for ix = 1:etframes
        
        writeVideo(vid,finalframe);
        
    end
    
    close(vid)
    
    else
    
    pause(etime)
    
end

%% Plotting, SLERP (delicious, spicy algorithm)

figure(4)

clf

pause(10^-8)
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);

% This maximizes the window for us automatically

subplot(1,2,2)

hold on

% Building final prism

patch('Vertices',verticesfinal,'Faces',faceorder,'FaceColor','flat','FaceVertexCData',col,'FaceAlpha',maxtrans,'EdgeColor','k','LineWidth',2);

% Miscellaneous plot commands to make everything look nicer (we're using
% an isometric 3D view angle)

axis equal
grid on

axis(axnew)

% Isometric view angle
view(-45,35.264)

xlabel('Global X axis')
ylabel('Global Y axis')
zlabel('Global Z axis')

title('Final Prism, SLERP')

subplot(1,2,1)

hold on

% Building initial prism

initpris = patch('Vertices',vertices,'Faces',faceorder,'FaceColor','flat','FaceVertexCData',col,'FaceAlpha',maxtrans,'EdgeColor','k','LineWidth',2);

% Miscellaneous plot commands to make everything look nicer (we're using
% an isometric 3D view angle)

axis equal
grid on

axis(axnew)

% Isometric view angle
view(-45,35.264)

xlabel('Global X axis')
ylabel('Global Y axis')
zlabel('Global Z axis')

title('Initial Prism, SLERP')

figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')
set(findall(figureHandle,'type','axes'),'fontSize',14,'fontWeight','bold')

% Larger text

if savevid
    
    vid = VideoWriter('rotationlin.avi');
    vid.FrameRate = 120;
    open(vid);

end

for ix = 1:tranits
    
    if ~savevid
        
        start(t)
        
    end
    
    set(initpris,'FaceAlpha',transvec(ix))
    
    if savevid
        
        writeVideo(vid,getframe(gcf));
        
    else
        
        wait(t)
        
    end
    
end

for ix = 1:interpsize
    
    if ~savevid
        
        start(t)
        
    end
    
    rotme = patch('Vertices',verticesslerp(:,:,ix),'Faces',faceorder,'FaceColor','flat','FaceVertexCData',col,'FaceAlpha',maxtrans,'EdgeColor','k','LineWidth',2);
    
    if savevid
        
        writeVideo(vid,getframe(gcf));
        
    else
        
        wait(t)
        
    end
    
    delete(rotme)
    
end

patch('Vertices',verticesslerp(:,:,end),'Faces',faceorder,'FaceColor','flat','FaceVertexCData',col,'FaceAlpha',maxtrans,'EdgeColor','k','LineWidth',2);

if savevid
        
    finalframe = getframe(gcf);
    
    for ix = 1:etframes
        
        writeVideo(vid,finalframe);
        
    end
    
    close(vid)
    
    else
    
    pause(etime)
    
end

%% Closer comparison of Euler and SLERP

figure(5)

clf

pause(10^-8)
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);

% This maximizes the window for us automatically

subplot(1,2,1)

hold on

% Building initial prism

patch('Vertices',vertices,'Faces',faceorder,'FaceColor','flat','FaceVertexCData',col,'FaceAlpha',mintrans,'EdgeColor','k','LineWidth',2);

% Miscellaneous plot commands to make everything look nicer (we're using
% an isometric 3D view angle)

axis equal
grid on

axis(axnew)

% Isometric view angle
view(-45,35.264)

xlabel('Global X axis')
ylabel('Global Y axis')
zlabel('Global Z axis')

title('Initial Prism, Linear Euler Angle Interpolation')

subplot(1,2,2)

hold on

% Building initial prism

patch('Vertices',vertices,'Faces',faceorder,'FaceColor','flat','FaceVertexCData',col,'FaceAlpha',mintrans,'EdgeColor','k','LineWidth',2);

% Miscellaneous plot commands to make everything look nicer (we're using
% an isometric 3D view angle)

axis equal
grid on

axis(axnew)

% Isometric view angle
view(-45,35.264)

xlabel('Global X axis')
ylabel('Global Y axis')
zlabel('Global Z axis')

title('Initial Prism, SLERP Interpolation')

figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')
set(findall(figureHandle,'type','axes'),'fontSize',14,'fontWeight','bold')

% Larger text

if savevid
    
    vid = VideoWriter('rotationscompare.avi');
    vid.FrameRate = 120;
    open(vid);
    
end

for ix = 1:interpsize
    
    if ~savevid
        
        start(t)
        
    end
    
    subplot(1,2,1)
    rotme1 = patch('Vertices',verticeslinEul(:,:,ix),'Faces',faceorder,'FaceColor','flat','FaceVertexCData',col,'FaceAlpha',maxtrans,'EdgeColor','k','LineWidth',2);
    subplot(1,2,2)
    rotme2 = patch('Vertices',verticesslerp(:,:,ix),'Faces',faceorder,'FaceColor','flat','FaceVertexCData',col,'FaceAlpha',maxtrans,'EdgeColor','k','LineWidth',2);
    
    if savevid
        
        writeVideo(vid,getframe(gcf));
        
    else
        
        wait(t)
        
    end
    delete(rotme1)
    delete(rotme2)
    
end

subplot(1,2,1)
patch('Vertices',verticeslinEul(:,:,end),'Faces',faceorder,'FaceColor','flat','FaceVertexCData',col,'FaceAlpha',maxtrans,'EdgeColor','k','LineWidth',2);
subplot(1,2,2)
patch('Vertices',verticesslerp(:,:,end),'Faces',faceorder,'FaceColor','flat','FaceVertexCData',col,'FaceAlpha',maxtrans,'EdgeColor','k','LineWidth',2);

if savevid
    
    finalframe = getframe(gcf);
    
    for ix = 1:etframes
    
        writeVideo(vid,finalframe);
    
    end
    
    close(vid)
    
else
    
    pause(etime)
    
end