function heu_val = heuristic_function(y_cand_m, goal_state)
    if y_cand_m(2) < goal_state(2) %off the table
        heu_val = 1000;
    else
        heu_val = 20 * sqrt( (goal_state(1)-y_cand_m(1))^2 + (goal_state(2)-y_cand_m(2))^2);
    end
end