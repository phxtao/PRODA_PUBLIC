function [evolved_ranked_candidate_cost_value, evolved_ranked_candidate_para_value, evolved_ranked_candidate_steady_state, ...
    evolved_ranked_candidate_mod_soc, evolved_ranked_candidate_soc_stock, evolved_ranked_candidate_soc_layer, ...
	evolved_ranked_candidate_r2, max_ss_cost_record_update] = ...
    fun_sce_cce(icomplex, npara, complex_num, point_num, parents_num, max_ss_cost_record, ... % cce variablee
    ranked_candidate_cost_value, ranked_candidate_index, ranked_candidate_para_value, ranked_candidate_steady_state, ...  % cce variablee
    ranked_candidate_mod_soc, ranked_candidate_soc_stock, ranked_candidate_soc_layer, ranked_candidate_r2, ... % cce variablee
    nbedrock, sand_vector, npp_mean, ... % mic model variable
    input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ... % mic model variable
    altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin, ... % mic model variable
    zsoi, soil_decom_num, wosis_layer_depth, wosis_layer_obs, layer_weight) % other variable

% step 4.1: assign weights
triangular_weight = 2*(point_num + 1 - (1:point_num))/(point_num*(point_num+1));

% step 4.2: select parents
subcomplex_point_loc = nan(parents_num, 1);

while true
    selected_samples = randsample(point_num, point_num*100, true, triangular_weight);
    if length(unique(selected_samples)) >= parents_num
        break
    end
end

for iparents = 1:parents_num
    subcomplex_point_loc(iparents) = selected_samples(1);
    selected_samples = selected_samples(selected_samples ~= selected_samples(1));
end

subcomplex_point_cost = ranked_candidate_cost_value(icomplex, subcomplex_point_loc);
subcomplex_point_index = ranked_candidate_index(icomplex, subcomplex_point_loc);
subcomplex_point_ss = ranked_candidate_steady_state(icomplex, subcomplex_point_loc);
subcomplex_point_r2 = ranked_candidate_r2(icomplex, subcomplex_point_loc);
% corresponding para values
subcomplex_para_value = reshape(ranked_candidate_para_value(icomplex, subcomplex_point_loc, :), [parents_num, npara]);
% corresponding mod soc
subcomplex_mod_soc = reshape(ranked_candidate_mod_soc(icomplex, subcomplex_point_loc, :), [parents_num, length(wosis_layer_obs)]);
% corresponding soc stock and layer
subcomplex_soc_stock = reshape(ranked_candidate_soc_stock(icomplex, subcomplex_point_loc, :), [parents_num, 5]);
subcomplex_soc_layer = reshape(ranked_candidate_soc_layer(icomplex, subcomplex_point_loc, :), [parents_num, soil_decom_num]);
% step 4.3: generatinbg offspring
[ranked_subcomplex_point_cost, order_index] = sort(subcomplex_point_cost, 'ascend');
ranked_subcomplex_point_index = subcomplex_point_index(order_index);
ranked_subcomplex_point_ss = subcomplex_point_ss(order_index);
ranked_subcomplex_point_r2 = subcomplex_point_r2(order_index);
ranked_subcomplex_para_value = subcomplex_para_value(order_index, :);
ranked_subcomplex_mod_soc = subcomplex_mod_soc(order_index, :);
ranked_subcomplex_soc_stock = subcomplex_soc_stock(order_index, :);
ranked_subcomplex_soc_layer = subcomplex_soc_layer(order_index, :);
ranked_subcomplex_point_loc = subcomplex_point_loc(order_index);

evolved_ranked_subcomplex_point_cost = ranked_subcomplex_point_cost;
evolved_ranked_subcomplex_point_ss = ranked_subcomplex_point_ss;
evolved_ranked_subcomplex_point_r2 = ranked_subcomplex_point_r2;
evolved_ranked_subcomplex_para_value = ranked_subcomplex_para_value;
evolved_ranked_subcomplex_mod_soc = ranked_subcomplex_mod_soc;
evolved_ranked_subcomplex_soc_stock = ranked_subcomplex_soc_stock;
evolved_ranked_subcomplex_soc_layer = ranked_subcomplex_soc_layer;

% step 4.3.1: computing centriod of parameter values in the subcomplex
para_centriod = nan(npara, 1);
for ipara = 1 : npara
    para_centriod(ipara) = (1/(parents_num-1)*sum(ranked_subcomplex_para_value(1:(parents_num-1), ipara)));
end
% step 4.3.2: new point (reflecttion step)
para_new = 2*para_centriod - ranked_subcomplex_para_value(parents_num, :)';
% step 4.3.3: if new point in para space
if min(para_new) > 0 && max(para_new) < 1
    % do nothing
else
    % (mutation step) propose new para according to the smallest hyper space occupied by the subcomplex
    para_new = nan(npara, 1);
    para_new_origin = rand(complex_num, point_num, npara);
    for ipara = 1 : npara
        sub_max = max(ranked_subcomplex_para_value(:, ipara));
        sub_min = min(ranked_subcomplex_para_value(:, ipara));
        para_new(ipara) = para_new_origin(ipara)*(sub_max - sub_min) + sub_min;
    end
end
% calculate cost function of new para
[~, soc_stock_summary, soc_mod, ss_index] = ...
    matrix_fun_seasonal(para_new, wosis_layer_obs, nbedrock, sand_vector, npp_mean, ...
    input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
    altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);

if isnan(sum(soc_mod(:, 5))) == 0
	optimize_profile_soc = interp1(zsoi(1:soil_decom_num, 1), soc_mod(:, 5), wosis_layer_depth, 'pchip');
else
	optimize_profile_soc = -10^7*ones(length(wosis_layer_depth), 1);
end

cost_new = cost_fun(layer_weight, optimize_profile_soc, wosis_layer_obs);

if ss_index == 0
    cost_new = cost_new + max_ss_cost_record(icomplex);
end

% step 4.3.4: explore new better para
if cost_new < ranked_subcomplex_point_cost(parents_num)
    % do nothing
else
    % propose another para (contraction step)
    para_new = (para_centriod + ranked_subcomplex_para_value(parents_num, :)')/2;
    
    % calculate cost function of new para
    [~, soc_stock_summary, soc_mod, ss_index] = ...
        matrix_fun_seasonal(para_new, wosis_layer_obs, nbedrock, sand_vector, npp_mean, ...
        input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
        altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
    
 
	if isnan(sum(soc_mod(:, 5))) == 0
		optimize_profile_soc = interp1(zsoi(1:soil_decom_num, 1), soc_mod(:, 5), wosis_layer_depth, 'pchip');
	else
		optimize_profile_soc = -10^7*ones(length(wosis_layer_depth), 1);
	end
    
    cost_new = cost_fun(layer_weight, optimize_profile_soc, wosis_layer_obs);
    if ss_index == 0
        cost_new = cost_new + max_ss_cost_record(icomplex);
    end
end

% step 4.3.5: continue explore new better para
if cost_new < ranked_subcomplex_point_cost(parents_num)
    evolved_ranked_subcomplex_point_cost(parents_num) = cost_new;
    evolved_ranked_subcomplex_para_value(parents_num, :) = para_new;
    evolved_ranked_subcomplex_mod_soc(parents_num, :) = optimize_profile_soc;
    evolved_ranked_subcomplex_soc_stock(parents_num, :) = soc_stock_summary;
	evolved_ranked_subcomplex_soc_layer(parents_num, :) = soc_mod(:, 5);
	evolved_ranked_subcomplex_point_ss(parents_num) = ss_index;
    mod_r2 = fun_r2(wosis_layer_obs, optimize_profile_soc);
    evolved_ranked_subcomplex_point_r2(parents_num) = mod_r2;
    
    if ss_index == 1
        max_ss_cost_record(icomplex) = max(max_ss_cost_record(icomplex), cost_new);
    end
else
    % (mutation step) propose new para according to the smallest hyper space occupied by the subcomplex
    para_new = nan(npara, 1);
    para_new_origin = rand(complex_num, point_num, npara);
    for ipara = 1 : npara
        sub_max = max(ranked_subcomplex_para_value(:, ipara));
        sub_min = min(ranked_subcomplex_para_value(:, ipara));
        para_new(ipara) = para_new_origin(ipara)*(sub_max - sub_min) + sub_min;
    end
    
    % calculate cost function of new para
    [~, soc_stock_summary, soc_mod, ss_index] = ...
        matrix_fun_seasonal(para_new, wosis_layer_obs, nbedrock, sand_vector, npp_mean, ...
        input_vector_cwd, input_vector_litter1, input_vector_litter2, input_vector_litter3, ...
        altmax_current_profile, altmax_lastyear_profile, soil_temp_profile, soil_water_profile, xio, xin);
    
	 
	if isnan(sum(soc_mod(:, 5))) == 0
		optimize_profile_soc = interp1(zsoi(1:soil_decom_num, 1), soc_mod(:, 5), wosis_layer_depth, 'pchip');
	else
		optimize_profile_soc = -10^7*ones(length(wosis_layer_depth), 1);
	end
 
    cost_new = cost_fun(layer_weight, optimize_profile_soc, wosis_layer_obs);
    if ss_index == 0
        cost_new = cost_new + max_ss_cost_record(icomplex);
    else
        max_ss_cost_record(icomplex) = max(max_ss_cost_record(icomplex), cost_new);
    end
    
    % replace with new para no matter better or not
    evolved_ranked_subcomplex_point_cost(parents_num) = cost_new;
    evolved_ranked_subcomplex_para_value(parents_num, :) = para_new;
    evolved_ranked_subcomplex_mod_soc(parents_num, :) = optimize_profile_soc;
    evolved_ranked_subcomplex_soc_stock(parents_num, :) = soc_stock_summary;
	evolved_ranked_subcomplex_soc_layer(parents_num, :) = soc_mod(:, 5);
	evolved_ranked_subcomplex_point_ss(parents_num) = ss_index;
    mod_r2 = fun_r2(wosis_layer_obs, optimize_profile_soc);
    evolved_ranked_subcomplex_point_r2(parents_num) = mod_r2;
end

% step 4.4: assign back costs and para values
ranked_subcomplex_point_cost = evolved_ranked_subcomplex_point_cost;
ranked_subcomplex_para_value = evolved_ranked_subcomplex_para_value;
ranked_subcomplex_mod_soc = evolved_ranked_subcomplex_mod_soc;
ranked_subcomplex_soc_stock = evolved_ranked_subcomplex_soc_stock;
ranked_subcomplex_soc_layer = evolved_ranked_subcomplex_soc_layer;
ranked_subcomplex_point_ss = evolved_ranked_subcomplex_point_ss;
ranked_subcomplex_point_r2 = evolved_ranked_subcomplex_point_r2;

subcomplex_point_cost(order_index) = ranked_subcomplex_point_cost;
subcomplex_para_value(order_index, :) = ranked_subcomplex_para_value;
subcomplex_mod_soc(order_index, :) = ranked_subcomplex_mod_soc;
subcomplex_soc_stock(order_index, :) = ranked_subcomplex_soc_stock;
subcomplex_soc_layer(order_index, :) = ranked_subcomplex_soc_layer;
subcomplex_point_ss(order_index) = ranked_subcomplex_point_ss;
subcomplex_point_r2(order_index) = ranked_subcomplex_point_r2;

ranked_candidate_cost_value(icomplex, subcomplex_point_loc) = subcomplex_point_cost;
ranked_candidate_para_value(icomplex, subcomplex_point_loc, :) = subcomplex_para_value;
ranked_candidate_mod_soc(icomplex, subcomplex_point_loc, :) = subcomplex_mod_soc;
ranked_candidate_soc_stock(icomplex, subcomplex_point_loc, :) = subcomplex_soc_stock;
ranked_candidate_soc_layer(icomplex, subcomplex_point_loc, :) = subcomplex_soc_layer;
ranked_candidate_steady_state(icomplex, subcomplex_point_loc) = subcomplex_point_ss;
ranked_candidate_r2(icomplex, subcomplex_point_loc) = subcomplex_point_r2;

evolved_ranked_candidate_cost_value = ranked_candidate_cost_value;
evolved_ranked_candidate_para_value = ranked_candidate_para_value;
evolved_ranked_candidate_mod_soc = ranked_candidate_mod_soc;
evolved_ranked_candidate_soc_stock = ranked_candidate_soc_stock;
evolved_ranked_candidate_soc_layer = ranked_candidate_soc_layer;
evolved_ranked_candidate_steady_state = ranked_candidate_steady_state;
evolved_ranked_candidate_r2 = ranked_candidate_r2;

max_ss_cost_record_update = max_ss_cost_record;
end
