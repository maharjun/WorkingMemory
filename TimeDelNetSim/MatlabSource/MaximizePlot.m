function [] = MaximizePlot()
	%MAXIMIZEPLOT Summary of this function goes here
	%   Detailed explanation goes here
	pause(0.00001);
	frame_h = get(handle(gcf),'JavaFrame');
	set(frame_h,'Maximized',1);
end

