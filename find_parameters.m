function [v_max, k_m] = find_parameters(v_0, init_S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 132 
% Program Description 
% Given the v_0 data with its corresponding initial substrate concentration
% values [S] (uM), find the parameters V_max and K_m for an enzyme by using
% the Eadie-Hofstee method: v_0 = V_max - K_m * v_0/[S]
%
% Function Call
% [v_max, k_m] = find_parameters(v_0, init_S)
%
% Input Arguments
% v_0 - DATA SET CONTAINING V_0 VALUES FOR AN ENZYME
% init_S - CORRESPONDING SUBSTRATE CONCENTRATION FOR THE V_0 DATA SET
%
% Output Arguments
% v_max - V_MAX PARAMETER FOR THE ENZYME
% k_m - K_M PARAMETER FOR THE ENZYME
%
% Assignment Information
%   Assignment:     M2, Algorithm Development
%   Team member:    Nathan Yao, nhyao@purdue.edu [repeat for each person]
%   Team ID:        004-19
%   Academic Integrity:
%     [X] We worked with one or more peers but our collaboration
%        maintained academic integrity.
%     Peers I worked with: Justin Schwartz, schwarjl@purdue.edu
%                          Kaushik Karthikeyan, karthik5@purdue.edu
%                          Patrick Mullen, mullenp@purdue.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ____________________
%% INITIALIZATION

fraction = v_0 ./ init_S; % obtain the v_0/[S] fraction term (independent variable)

%% ____________________
%% CALCULATIONS

% Eadie-Hofstee method: v_0 = V_max - K_m * v_0/[S]
coeffs = polyfit(fraction, v_0, 1);
k_m = coeffs(1) * -1; % slope will be negative by default.
v_max = coeffs(2);

%% ____________________
%% FORMATTED TEXT/FIGURE DISPLAYS


%% ____________________
%% RESULTS


%% ____________________
%% ACADEMIC INTEGRITY STATEMENT
% We have not used source code obtained from any other unauthorized
% source, either modified or unmodified. Neither have we provided
% access to my code to another. The program we are submitting
% is our own original work.



