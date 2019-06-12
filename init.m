% clear all;
clear;
close all;
clc;
dbstop if error;

current_tests=[ pwd '/current_tests' ];
fct = genpath([ pwd '/functions' ]);
mains = [ pwd '/mains'];
addpath(pwd)
addpath(fct)
% addpath(stk_folder)
addpath current_tests mains stk-2.3.4;
% stk_init;
clear current_tests mains fct
home;
% try
%     pdeinit
% end