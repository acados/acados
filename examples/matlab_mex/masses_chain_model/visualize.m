% script to visualize the chain of masses

%drawnow update
drawnow('expose')

figure(1), set(gcf, 'Color','white');
clf

subplot(3,3,[1:6]);
tol = 0.00;
p = patch([-0.2, 1.2, 1.2, -0.2], [wall-tol, wall-tol, wall-tol, wall-tol], [-4, -4, 1, 1], 'g');
hold on;



tmp_pos = reshape(cur_pos, 3, length(cur_pos)/3);
cur_pos = zeros(3,1);
for ii=1:nfm
	cur_pos = [cur_pos, tmp_pos(:,1+2*(ii-1))];
end

plot3(cur_pos(1,:), cur_pos(2,:), cur_pos(3,:), '-ob', ...
'MarkerSize', 7.5, 'MarkerFaceColor', 'b', 'linewidth', 0.2);

view([-135 45*3/4]);

xlim([-0.2 1.2]);
ylim([-0.2 1.2]);
zlim([-4 1]);

grid on;

set(gca, 'Box', 'on');

xlabel( 'x [m]' );
ylabel( 'y [m]' );
zlabel( 'z [m]' );
