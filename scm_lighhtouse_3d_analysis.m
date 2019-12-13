close all 

POSTCAL = true;
FINAL = true;
FILTER = true;

lh_data_file = 'SCM_lighthouse_data_log_20190806T133924';
frame_rate = 119.997597;
optitrak_t = scumlighthouserftake1optitrack.VarName1 ./frame_rate;
scum_quat = quaternion(scumlighthouse006.RigidBody3,scumlighthouse006.RigidBody,scumlighthouse006.RigidBody1,scumlighthouse006.RigidBody2);
scum_position = [scumlighthouse006.RigidBody4,scumlighthouse006.RigidBody5,scumlighthouse006.RigidBody6];

lh_quat = quaternion(mean(scumlighthouse006.RigidBody11,'omitnan'),mean(scumlighthouse006.RigidBody8,'omitnan'),mean(scumlighthouse006.RigidBody9,'omitnan'),mean(scumlighthouse006.RigidBody10,'omitnan'));
lh_position = [scumlighthouse006.RigidBody12,scumlighthouse006.RigidBody13,scumlighthouse006.RigidBody14];

figure
plot3(lh_position(:,1), lh_position(:,2), lh_position(:,3))
hold on

plot3(scum_position(:,1), scum_position(:,2), scum_position(:,3))

theta = 3.14159/2
Ry=[cos(theta), 0, sin(theta);
    0, 1, 0;
    -sin(theta),0,cos(theta)]

Rx=[1, 0, 0;
    0, cos(theta), -sin(theta);
    0,sin(theta),cos(theta)]

scum_rotated = scum_position*Ry'*Rx';
lh_rotated = lh_position*Ry'*Rx';

plot3(scum_rotated(:,1), scum_rotated(:,2), scum_rotated(:,3))
plot3(lh_rotated(:,1), lh_rotated(:,2), lh_rotated(:,3))
figure

plot( -scum_position(1:5000,1),scum_position(1:5000,2))

lh_rot = rotmat(lh_quat,'frame')

%calculate plane from 3 points on the back of lighthouse in room frame
m1 = Rx*Ry*[mean(scumlighthouse006.RigidBodyMarker16,'omitnan'),mean(scumlighthouse006.RigidBodyMarker17,'omitnan'),mean(scumlighthouse006.RigidBodyMarker18,'omitnan')]'
m2 =  Rx*Ry*[mean(scumlighthouse006.RigidBodyMarker20,'omitnan'),mean(scumlighthouse006.RigidBodyMarker21,'omitnan'),mean(scumlighthouse006.RigidBodyMarker22,'omitnan')]'
m3 =  Rx*Ry*[mean(scumlighthouse006.RigidBodyMarker24,'omitnan'),mean(scumlighthouse006.RigidBodyMarker25,'omitnan'),mean(scumlighthouse006.RigidBodyMarker26,'omitnan')]'
data = [m1';m2';m3']

%plane equation
vec1 = m1 - m2
vec2 = m3 - m2

normal = cross(vec1',vec2')
figure 
scatter3(data(:,1), data(:,2), data(:,3))
hold on 
axis equal
quiver3(mean(data(:,1)),mean(data(:,2)), mean(data(:,3)),normal(1), normal(2), normal(3))
quiver3(mean(data(:,1)),mean(data(:,2)), mean(data(:,3)),vec1(1), vec1(2), vec1(3))
quiver3(mean(data(:,1)),mean(data(:,2)), mean(data(:,3)),vec2(1), vec2(2), vec2(3))

%calculate rotation of lighthouse, earth to lh
if POSTCAL
    lh_x_correction = -50e-3;
else
    lh_x_correction = 0
end
lh_position = [mean(data(:,1))+lh_x_correction,mean(data(:,2)), mean(data(:,3))]
lh_rotation_earth_to_lh = vrrotvec2mat(vrrotvec([1,0,0],normal))
dot(normal,[-1,0,0])
theta = acos(dot(normal,[-1,0,0])/norm(normal))

if POSTCAL
    theta_offset = deg2rad(-.5);
else
    theta_offset = 0;
    
end
theta = theta+theta_offset;
lh_rotation_earth_to_lh = [cos(theta), 0, sin(theta);
    0, 1, 0;
    -sin(theta),0,cos(theta)]

%rotate by 180 over z
if POSTCAL    
    alignment_offset = -2.5;
else
    alignment_offset = 0;
end
angle = deg2rad(180+alignment_offset);
rz180 = [cos(angle),-sin(angle),0;
        sin(angle),cos(angle),0;
        0,0,1];
size(scum_rotated)
scum_lhframe = (scum_rotated - repmat(lh_position,size(scum_rotated(:,1))) )*lh_rotation_earth_to_lh'*rz180'; % rotate scum earth to lh frame

%calc elevation and azimuth
elevation = atan2(scum_lhframe(:,3),scum_lhframe(:,1));
azimuth = atan2(scum_lhframe(:,2),scum_lhframe(:,1));

%plot scum trajectory
figure
plot3(scum_rotated(start_idx:end_idx,1), scum_rotated(start_idx:end_idx,2), scum_rotated(start_idx:end_idx,3),'LineWidth',2)
hold on
%plot3(scum_lhframe(:,1), scum_lhframe(:,2), scum_lhframe(:,3))
scatter3(lh_position(1),lh_position(2), lh_position(3))
axis equal
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Experimental Setup with Receiver Trajectory')
set(gca,'FontSize',16)
legend('Receiver Trajectory','Lighthouse Location')

%% data plots

%plot azimuth
figure
%plot(azimuth(1:5000),elevation(1:5000))

%filter the data
clear azimuth_degrees_filtered
clear elevation_degrees_filtered
clear t_filtered
azimuth_degrees_filtered(1) = azimuth_degrees(1) ;
elevation_degrees_filtered(1) = elevation_degrees(1);
t_filtered(1) = t(1)
length(azimuth_degrees(1,:))
filtered_idx = 1;
for i = 2:length(azimuth_degrees)
    
    %max diff filter
    if abs(azimuth_degrees(i) - azimuth_degrees(i-1)) < 5 && abs(elevation_degrees(i) - elevation_degrees(i-1)) < 5
        filtered_idx = filtered_idx + 1;
        azimuth_degrees_filtered(filtered_idx) = azimuth_degrees(i)
        elevation_degrees_filtered(filtered_idx) = elevation_degrees(i)
        t_filtered(filtered_idx) = t(i)
    end
end

%plot azimuths
figure
subplot(2,1,1)
plot(t_filtered,azimuth_degrees_filtered)
set(gca,'Fontsize',16)
grid on
xlabel('Time [s]')
ylabel('Azimuth [deg]')

subplot(2,1,2)
plot(optitrak_t,azimuth)
set(gca,'Fontsize',16)
grid on
xlabel('Time [s]')
ylabel('Azimuth [deg]')

%offset (optitrak - lighthouse) from plotting
offset = 20.2921 - 10.0346;
time_end = 31.5
%plot azimuths again with offset
figure
subplot(2,1,1)
plot(t_filtered,azimuth_degrees_filtered,'LineWidth',2)
set(gca,'Fontsize',16)
grid on
xlabel('Time [s]')
ylabel('Azimuth [deg]')

hold on
plot(optitrak_t-offset,rad2deg(azimuth) +90,'--','LineWidth',2)
set(gca,'Fontsize',16)
grid on
xlabel('Time [s]')
ylabel('Azimuth [deg]')
axis([0 31 60 140]) 
legend('Lighthouse Tracking', 'Ground Truth')
if POSTCAL
    if FINAL
        title('Lighthouse Tracking of Azimuth and Elevation')
    else
        title('2 deg Yaw Correction, 0.5 Degree pitch correction, 5 cm centroid correction')
    end
    
else 
    title('No Corrections Performed')
end

%plot elevations with offset
subplot(2,1,2)
plot(t_filtered,elevation_degrees_filtered,'LineWidth',2)
set(gca,'Fontsize',16)
grid on
xlabel('Time [s]')
ylabel('Elevation [deg]')

hold on
plot(optitrak_t-offset,rad2deg(elevation) +90, '--','LineWidth',2)
set(gca,'Fontsize',16)
grid on
xlabel('Time [s]')
ylabel('Elevation [deg]')
axis([0 31 60 140]) 
legend('Lighthouse Tracking', 'Ground Truth')

figure
% Plot trajectory

tracking_start = 112;
scatter(azimuth_degrees_filtered(tracking_start:end),elevation_degrees_filtered(tracking_start:end),'.')
xlabel('Azimuth [deg]');
ylabel('Elevation [deg]');
set(gca,'Fontsize',16)
% axis([0 180 0 180])
grid on
axis square
axis([0 180 0 180])
% ax = gca;
% ax.LineWidth = 2
xticks([0:30:180])
yticks([0:30:180])

hold on 
start_idx = round(offset*frame_rate)
end_idx = round((time_end+offset)*frame_rate)
synchronized_az = azimuth(start_idx:end_idx);
synchronized_el = elevation(start_idx:end_idx);

scatter(rad2deg(synchronized_az)+90,rad2deg(synchronized_el)+90,'.');
legend('Lighthouse Tracking','Ground Truth')

if POSTCAL
    if FINAL
        title('Lighthouse Tracking of Trajectory')
    else
        title('2 deg Yaw Correction, 0.5 Degree pitch correction, 5 cm centroid correction')
    end
    
else
    title('No Corrections Performed')
end

%plot error
figure
truth = [rad2deg(synchronized_az)+90,rad2deg(synchronized_el)+90];
tracking = [azimuth_degrees_filtered(tracking_start:end)',elevation_degrees_filtered(tracking_start:end)'];
size(tracking)
size(truth)
truth_interp = [interp1(optitrak_t(start_idx:end_idx)-offset,rad2deg(synchronized_az)+90,t_filtered(tracking_start:end))',interp1(optitrak_t(start_idx:end_idx)-offset,rad2deg(synchronized_el)+90,t_filtered(tracking_start:end))'];
plot(t_filtered(tracking_start:end),tracking-truth_interp, 'LineWidth',2)
set(gca,'FontSize',16)
xlabel("Time (s)")
ylabel("Tracking Error (degrees)")
legend("Azimuth Error", "Elevation Error")


mean(rms(tracking-truth_interp,2),'omitnan');
sqrt(mean((tracking(:,1)-truth_interp(:,1)).^2,'omitnan'))
sqrt(mean((tracking(:,2)-truth_interp(:,2)).^2,'omitnan'))
distance = vecnorm(tracking-truth_interp,2,2);
sqrt(mean(distance.^2,'omitnan'))

el_err = tracking(:,2)-truth_interp(:,2);
filtered_idx = find(abs(el_err)<5);
size(el_err);
bad_idx = find(abs(el_err)>5);
sqrt(mean(el_err(filtered_idx).^2))


%plot y z 
tracking-truth_interp
diff = scum_rotated(start_idx:end_idx,:)-repmat(lh_position,size(scum_rotated(start_idx:end_idx,1)))
r = vecnorm(diff(:,1),2,2);
r_interp = interp1(optitrak_t(start_idx:end_idx)-offset,r,t_filtered(tracking_start:end));

z = tan(deg2rad(tracking(:,2)-98)).*r_interp' + lh_position(3);
y = tan(deg2rad(tracking(:,1)-89)).*r_interp' + lh_position(2);

figure
plot(y,z)

hold on 
plot(-scum_rotated(start_idx:end_idx,2),scum_rotated(start_idx:end_idx,3))