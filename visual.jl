# basic plot settings
Plots.resetfontsizes()
Plots.scalefontsizes(1.25)

# """
#     Regression plot of Ignition Delay Times (IDTs)
#     It will plot the correlations of training data and prediciton data in validation set
# """
function regression_plot(; max=10, name="")
    @printf("...regression plot ...\n")
    data = conds[n_train.+(1:minimum([max,n_exp])), end-2:end]
    plt = plot(xlabel="Training IDT [s]", ylabel="Prediction IDT [s]", legend=:topleft)
    Plots.scatter!(plt, data[:, 1], data[:, 3], xscale=:log10, yscale=:log10, label="optimized")
    Plots.scatter!(plt, data[:, 1], data[:, 2], xscale=:log10, yscale=:log10, label="initial")
    plot!(plt, [minimum(data), maximum(data)], [minimum(data), maximum(data)], label="y=x")
    png(plt, string(fig_path, "/regression$name"))
end

# """
#     Regression plot of Laminar Flame Speeds (Sus)
#     It will plot the correlations of training data and prediciton data in validation set
# """
function regression_f_plot(; max=10, name="")
    @printf("...regression flame plot ...\n")
    plt = plot(xlabel="Training Su [m/s]", ylabel="Prediction Su [m/s]", legend=:topleft)
    data = conds_f[1:minimum([max,n_exp_f]), end-2:end]
    Plots.scatter!(plt, data[:,1], data[:,3], xscale=:log10, yscale=:log10, label="optimized")
    Plots.scatter!(plt, data[:,1], data[:,2], xscale=:log10, yscale=:log10, label="initial")
    plot!(plt, [minimum(data), maximum(data)], [minimum(data), maximum(data)], label="y=x")
    png(plt, string(fig_path, "/regression$(name)_f"))
end

# """
#     validation_plot, will plot figures in $expr_name/figs/condtions
#     Each figure contains profiles of species and temperature, including profiles of
#         master_mech with p .= zeros  (detailed mech)
#         target_mech with p .= zeros  (skeletal mech)
#         target_mech with p           (optimized mech)
#     at different conditions.
#     By default, 10 conditions in validation set will be plotted.
# """
function validation_plot(;max=10, name="")
    @printf("...validation plot ...\n")
    max = minimum([max, n_exp - n_train]);

    # get profiles data of master mech
    if isfile("$ckpt_path/valid_master.bson")
        @load "$ckpt_path/valid_master.bson" vm_profiles vm_max;
        @printf("valid_master.bson loaded with max=%d\n", vm_max);
    else
        vm_profiles = []
        vm_max = 0
    end
    if vm_max < max
        global gas = CreateSolution(master_mech);
        global ns = gas.n_species;
        global nr = gas.n_reactions;
        for i_exp = n_train .+ ((vm_max+1):max)
            ts, pred = get_Tcurve_idt(conds[i_exp,:], zeros(nr*npr); dT=dT, tfinal=10.0)
            push!(vm_profiles, [ts, pred])
            @printf("vm %d idt %.2e\n", i_exp, ts[end])
        end
        vm_max = max
        @save "$ckpt_path/valid_master.bson" vm_profiles vm_max
    end

    # get profiles data of target mech with initial p
    if isfile("$ckpt_path/valid_target.bson")
        @load "$ckpt_path/valid_target.bson" vt_profiles vt_max
        @printf("valid_target.bson loaded with max=%d\n", vt_max);
    else
        vt_profiles = []
        vt_max = 0
    end
    if vt_max < max
        global gas = CreateSolution(target_mech);
        global ns = gas.n_species;
        global nr = gas.n_reactions;
        for i_exp = n_train .+ ((vt_max+1):max)
            ts, pred = get_Tcurve_idt(conds[i_exp,:], zeros(nr*npr); dT=dT, tfinal=10.0)
            push!(vt_profiles, [ts, pred])
            @printf("vt %d idt %.2e\n", i_exp, ts[end])
        end
        vt_max = max
        @save "$ckpt_path/valid_target.bson" vt_profiles vt_max
    end

    # get profiles data of target mech with current p
    global gas = CreateSolution(target_mech);
    global ns = gas.n_species;
    global nr = gas.n_reactions;
    ct_profiles = []
    for i_exp = n_train .+ (1:max)
        ts, pred = get_Tcurve_idt(conds[i_exp,:], p; dT=dT, tfinal=10.0)
        push!(ct_profiles, [ts, pred])
        @printf("ct %d idt %.2e\n", i_exp, ts[end])
    end

    # compare each profile and save
    for i_exp = n_train .+ (1:max)
        l_plt = []
        phi, Pj, T0 = conds[i_exp, 1:3];

        # prepare figure handles
        plt = plot(xlabel="Time [s]", ylabel="Temperature [K]", xscale=:log10,
                   title=@sprintf("%.1f K, %.1f atm, phi=%.1f", T0, Pj, phi))
        push!(l_plt, plt)
        for (s,v) in fuel
            plt = plot(xlabel="Time [s]", ylabel="Y $s", xscale=:log10)
            push!(l_plt, plt)
        end

        # plot master profile
        label = "Ref";
        line = (2, :dot);
        gas = CreateSolution(master_mech);
        ts, pred = vm_profiles[i_exp-n_train];
        ts .+= 1.e-6
        plot!(l_plt[1], ts, pred[end, :], l=line, label=label)
        for (i,(s,v)) in enumerate(fuel)
            plot!(l_plt[i+1], ts, pred[species_index(gas, "$s"), :], l=line, label=label)
        end

        # plot target profile with initial p
        label = "Sk";
        line = (2, :dash);
        gas = CreateSolution(target_mech);
        ts, pred = vt_profiles[i_exp-n_train];
        ts .+= 1.e-6
        plot!(l_plt[1], ts, pred[end, :], l=line, label=label)
        for (i,(s,v)) in enumerate(fuel)
            plot!(l_plt[i+1], ts, pred[species_index(gas, "$s"), :], l=line, label=label)
        end

        # plot target profile with current p
        label = "Train";
        line = (2, :solid);
        gas = CreateSolution(target_mech);
        ts, pred = ct_profiles[i_exp-n_train];
        ts .+= 1.e-6
        plot!(l_plt[1], ts, pred[end, :], l=line, label=label)
        for (i,(s,v)) in enumerate(fuel)
            plot!(l_plt[i+1], ts, pred[species_index(gas, "$s"), :], l=line, label=label)
        end

        pltsum = plot(l_plt..., framestyle=:box, legend=false)
        png(pltsum, string(fig_path, "/conditions/valid$(name)_$i_exp"))
    end
end

# """
#     Plot utils for profiles comparison
# """
function compare_profiles(figs, ts, pred, sp_names, line, color, name)
    plot!(figs[1], ts, pred[end,:], msw=0, l=line, c=color, lw=2, label=name)
    for (i,sp) in enumerate(sp_names)
        plot!(figs[1+i], ts, pred[species_index(gas, sp),:],
                xlim=[ts[end]*0.0, ts[end]], msw=0, l=line, c=color, lw=2, xscale=:identity)
    end
end

# """
#     Update parameters and save to mech
# """
function updateyaml(mech, p; name="_op")
    yaml = YAML.load_file(mech)
    n_species = length(yaml["phases"][1]["species"])
    n_reactions = length(yaml["reactions"])
    species_names = yaml["phases"][1]["species"]
    elements = yaml["phases"][1]["elements"]
    n_elements = length(elements)

    if length(p) == n_reactions
        _p = zeros(n_reactions, 3)
        _p[:,1] = p
    else
        _p = reshape(p, n_reactions, 3)
    end

    for i = 1:n_reactions
        p_vec = _p[i, :]
        fA = exp(p_vec[1])
        fb = p_vec[2]
        fEa = p_vec[3]
        reaction = yaml["reactions"][i]
        if haskey(reaction, "type")
            if (reaction["type"] == "falloff")
                yaml["reactions"][i]["low-P-rate-constant"]["A"] *= fA
                yaml["reactions"][i]["low-P-rate-constant"]["b"] += fb
                yaml["reactions"][i]["low-P-rate-constant"]["Ea"] += fEa
                yaml["reactions"][i]["high-P-rate-constant"]["A"] *= fA
                yaml["reactions"][i]["high-P-rate-constant"]["b"] += fb
                yaml["reactions"][i]["high-P-rate-constant"]["Ea"] += fEa
            else
                yaml["reactions"][i]["rate-constant"]["A"] *= fA
                yaml["reactions"][i]["rate-constant"]["b"] += fb
                yaml["reactions"][i]["rate-constant"]["Ea"] += fEa
            end
        else
            yaml["reactions"][i]["rate-constant"]["A"] *= fA
            yaml["reactions"][i]["rate-constant"]["b"] += fb
            yaml["reactions"][i]["rate-constant"]["Ea"] += fEa
        end
    end

    YAML.write_file(mech[1:end - 5] * "$name.yaml", yaml)
end
