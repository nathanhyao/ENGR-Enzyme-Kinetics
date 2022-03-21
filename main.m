function main()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 132 
% Program Description 
% Given the generation of enzyme desired and an initial [S], pass the
% appropriate data into a function that will calculate v_0, v_max, and k_m.
%
% Function Call
% main()
%
% Input Arguments
% NONE
%
% Output Arguments
% NONE
%
% Assignment Information
%   Assignment:     M2 and M3, Algorithm Development
%   Team member:    Patrick Mullen, mullenp@purdue.edu [repeat for each person]
%   Team ID:        004-19
%   Academic Integrity:
%     [X] We worked with one or more peers but our collaboration
%        maintained academic integrity.
%     Peers I worked with: Justin Schwartz, schwarjl@purdue.edu
%                          Nathan Yao, nhyao@purdue.edu
%                          Kaushik Karthikeyan, karthik5@purdue.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ____________________
%% INITIALIZATION

% import the raw experiment data
raw_data = readmatrix("Data_nextGen_KEtesting_allresults.csv"); % imports the excel data
time = raw_data(3:end, 1); % separate a column for time values
init_S_data = sort(raw_data(1, 2:21)); % get a list of [S] values to associate with v_0
init_S_data_unique = raw_data(1, 2:11); % get the unique [S] 

dist_duplicate = 10; % add this many columns to go from test to duplicate
dist_enzyme = 20; % add this many columns to reach the next enzyme's data

%% ____________________
%% CALCULATIONS

% algorithm to determine how many columns to locate data for each initial [S]

init_S_search = 1;
enzyme_search = 1;
v_0_data = zeros(1, 20); % preallocate memory for v_0_data for speed

v_0_insert = 1; % insert first v_0 data at index 1

while (enzyme_search <= 5) % find 10 v_0 values (of different initial [S]) for a given enzyme
    % prepare test and test duplicate data for UDF
    test = raw_data(3:end, 2 + dist_enzyme * (enzyme_search - 1) + (init_S_search - 1));
    duplicate = raw_data(3:end, 2 + dist_enzyme * (enzyme_search - 1) + (init_S_search - 1) + dist_duplicate);
    
    % pass data into UDF to get v_0 for given [S]
    v_0_test = find_v_0(time, test);
    v_0_duplicate = find_v_0(time, duplicate);
    
    % add results for v_0 into a data set
    v_0_data(v_0_insert) = v_0_test;
    v_0_data(v_0_insert + 1) = v_0_duplicate;
    v_0_insert = v_0_insert + 2; % update the insertion index (+2 to make room for test and its duplicate)
        
%% ____________________
%% FORMATTED TEXT/FIGURE DISPLAYS

    % create figures to display the v_0 vs [P] plot for each initial [S].
    % Plot test and its duplicate separately.
    
    % create the figure
    figure(enzyme_search + 5)
    subplot(3, 4, init_S_search)
    
    % readjust data sets for time, test, and duplicate for plot formatting
    test = rmmissing(test);
    time_plotting = time(1:length(test));
    plot(time_plotting, test, "b.")
    
    % choose an appropriate title to display based on initial [S]
    if (init_S_search == 1)
        title("Test 1: [S] = 3.75 \muM")
    elseif (init_S_search == 2)
        title("Test 2: [S] = 7.5 \muM")
    elseif (init_S_search == 3)
        title("Test 3: [S] = 15 \muM")
    elseif (init_S_search == 4)
        title("Test 4: [S] = 30 \muM")
    elseif (init_S_search == 5)
        title("Test 5: [S] = 65 \muM")
    elseif (init_S_search == 6)
        title("Test 6: [S] = 125 \muM")
    elseif (init_S_search == 7)
        title("Test 7: [S] = 250 \muM")
    elseif (init_S_search == 8)
        title("Test 8: [S] = 500 \muM")
    elseif (init_S_search == 9)
        title("Test 9: [S] = 1000 \muM")
    else
        title("Test 10: [S] = 2000 \muM")
    end
    
    xlabel("Time (s)")
    ylabel("[P] (\muM)")
    grid on
    
    % create a v_0 line to plot on the same axes. Only use results from the
    % original test and not the duplicate for simplicity
    hold on
    time_plotting = time(1:floor(length(time) / 5));
    v_0_line = v_0_test * time_plotting;
    plot(time_plotting, v_0_line, "r-")
    legend("Test", "v_0 Line", "Location", "Southeast")
    
    % choose an appropriate subplot title
    if (enzyme_search == 1)
        sgtitle("Enzyme NextGen-A Time vs [P] (\muM)");
    elseif (enzyme_search == 2)
        sgtitle("Enzyme NextGen-B Time vs [P] (\muM)");
    elseif (enzyme_search == 3)
        sgtitle("Enzyme NextGen-C Time vs [P] (\muM)");
    elseif (enzyme_search == 4)
        sgtitle("Enzyme NextGen-D Time vs [P] (\muM)");
    else
        sgtitle("Enzyme NextGen-E Time vs [P] (\muM)");
    end
    
    init_S_search = init_S_search + 1; % move onto the next initial [S]

    if (init_S_search > 10) % move onto the next enzyme if all initial [S] have been exhausted
        
        % find the parameters v_max and k_m
        [v_max, k_m] = find_parameters(v_0_data, init_S_data);
        
        fprintf("Results for Enzyme NextGen-"); % inform which enzyme it is
        if (enzyme_search == 1)
            fprintf("A\n");
        elseif (enzyme_search == 2)
            fprintf("B\n");
        elseif (enzyme_search == 3)
            fprintf("C\n");
        elseif (enzyme_search == 4)
            fprintf("D\n");
        else
            fprintf("E\n");
        end
        
        j = 1; % initialize index variable for displaying v_0 values
        for s = init_S_data_unique
            fprintf("Initial [S] = %8.2f uM yields approx. v_0 = %.3f uM/s\n", s, (v_0_data(j) + v_0_data(j + 1)) / 2);
            j = j + 2;
        end
        
        % display values for parameters V_max and K_m
        fprintf("V_max = %.3f (uM/s)\n", v_max);
        fprintf("K_m = %.3f (uM)\n", k_m);
        
        % display SSE error for model vs raw data
        menten_model = (v_max * init_S_data) ./ (k_m + init_S_data);
        SSE = sum((v_0_data - menten_model).^2);
        fprintf("SSE = %.5f\n", SSE);
        
        % create a Michaelis Menten plot
        figure(enzyme_search)
        plot(init_S_data, v_0_data, "bo")
        
        % choose an appropriate title
        if (enzyme_search == 1)
            title("Enzyme NextGen-A Michaelis-Menten Plot")
        elseif (enzyme_search == 2)
            title("Enzyme NextGen-B Michaelis-Menten Plot")
        elseif (enzyme_search == 3)
            title("Enzyme NextGen-C Michaelis-Menten Plot")
        elseif (enzyme_search == 4)
            title("Enzyme NextGen-D Michaelis-Menten Plot")
        else
            title("Enzyme NextGen-E Michaelis-Menten Plot")
        end
        
        % additional plot formatting
        xlabel("Initial [S] (\muM)")
        ylabel("Initial Reaction Velocity (\muM/s)")
        grid on
        hold on % plot the model on the same axes
        plot(init_S_data, menten_model)
        legend("Raw Data", "Model", "Location", "Southeast")
        
        % reset index and figure number variables for next enzyme
        enzyme_search = enzyme_search + 1;
        init_S_search = 1;
        v_0_insert = 1; % insert first v_0 data at index 1
    end
end

%% ____________________
%% RESULTS


%% ____________________
%% ACADEMIC INTEGRITY STATEMENT
% We have not used source code obtained from any other unauthorized
% source, either modified or unmodified. Neither have we provided
% access to my code to another. The program we are submitting
% is our own original work.
