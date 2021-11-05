clc, clear all, close all,

fprintf('\n testing bash job..')

initial_time    = 0
step            = 2 
final_time    	= 30

for n  = 4:1:16 
	example_butter(n,initial_time,final_time,step);
end


