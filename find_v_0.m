function [v_0] = find_v_0(time, test)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 132 
% Program Description 
% Given the time, test data, and test duplicate data for a specific enzyme
% working at initial [S] (uM), find its initial reaction velocity v_0.
%
% Function Call
% [v_0] = find_v_0(time, test)
%
% Input Arguments
% time - DATA SET OF TIME ELAPSED (SEC)
% test - DATA SET OF TEST'S PRODUCT CONCENTRATION [P] (uM)
%
% Output Arguments
% v_0 - INITIAL REACTION VELOCITY OF SPECIFIED ENZYME AND INITIAL [S]
%
% Assignment Information
%   Assignment:     M2, Algorithm Development
%   Team member:    Justin Schwartz, schwarjl@purdue.edu [repeat for each person]
%   Team ID:        004-19
%   Academic Integrity:
%     [X] We worked with one or more peers but our collaboration
%        maintained academic integrity.
%     Peers I worked with: Nathan Yao, nhyao@purdue.edu
%                          Patrick Mullen, mullenp@purdue.edu
%                          Kaushik Karthikeyan, karthik5@purdue.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ____________________
%% INITIALIZATION

test = rmmissing(test); % remove empty rows of blank cells from test
test_len = length(test); % see how many rows there are in the test

while_index = 0; % initiate while loop index at 0;
r_squared = 0; % initialize to zero before while loop algorithm
r_squared_threshold = 0.995; % choose 0.995 which indicates a near-perfect fit

v_0_eval = 1; % to find v_0, evaluate at time = 1 sec, when the experiment just begins

%% ____________________
%% CALCULATIONS

% "Trial and error" loop that removes data points until a quadratic polyfit
% function accurately fits the remaining data. Take the derivative of the
% quadratic function and evaluate it at time = 1 s to find v_0

window = 8; % take a moving average to reduce noise in data
test_noise_reduction = movmean(test(window:end), window);
test_noise_reduction = [test(1:(window - 1)); test_noise_reduction]; % avoid taking a moving average of the front

while (r_squared < r_squared_threshold)
    
    % consider fewer data points for each iteration
    limiter = test_len - while_index;
    test = test_noise_reduction(1:limiter);
    time = time(1:limiter);
    
    % implement a quadratic polyfit function
    coeffs = polyfit(time, test, 2);
    quad_model = coeffs(1) * time.^2 + coeffs(2) * time + 0; % force the intercept to be zero
    
    SSE = sum((test - quad_model).^2); % use least squares regression to test quality of fit
    SST = sum((test - mean(test)).^2); 
    r_squared = 1 - SSE / SST;
    
    while_index = while_index + 1;
end

% take the derivative of the quadratic model and evaluate it at time = 1 s.
% Technique involves differentiation power rule
v_0 = 2 * coeffs(1) * v_0_eval + coeffs(2);

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

