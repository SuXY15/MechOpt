include("header.jl")
include("simulator.jl")
include("visual.jl")

phi_arr = minimum(phi_arr):0.1:maximum(phi_arr);

using PyCall
global ct = pyimport("cantera");
global ct_gas = ct.Solution(mech_arr[1]);
global npr=1;

## prepare figure handles
fig = plot(legend=:topleft, legendfontsize=10, foreground_color_legend=nothing,
           xlabel="phi", ylabel="Su [cm/s]");
lines = [:scatter, :dash, :solid];
colors = [:black, :red, :blue];
symbols = [:square, :circle, :diamond];

## get IDTs of each mechanism under each condition
for m in 1:length(mech_arr)
    name = name_arr[m];
    mech = mech_arr[m];

    if is_restart && isfile("$exp_path/Sus_$name.bson")
        println("Loading Sus for $name ...");
        @load "$exp_path/Sus_$name.bson" Sus;
    else
        println("Calculating Sus for $name ...");

        global ct_gas = ct.Solution(mech);
        global gas = CreateSolution(mech);
        global ns = gas.n_species;
        global nr = gas.n_reactions;

        Sus = zeros(length(phi_arr), length(P_arr));
        t_costs = zeros(length(phi_arr), length(P_arr));

        for (k,phi) in enumerate(phi_arr)
            for (j,Pj) in enumerate(P_arr)
                t_cost = @elapsed Su = get_Su([phi, Pj], zeros(nr*npr));
                Sus[k,j] = Su;
                t_costs[k,j] = t_cost;
                @printf("phi %.1f P %2d [atm] Su = %.1e [cm/s] \n",
                         phi,     Pj,         Su*100);
            end
        end
        @save "$exp_path/Sus_$name.bson" Sus t_costs;
    end

    println("Ploting Sus for $(name) ... \n");
    for (j,Pj) in enumerate(P_arr)
        if j==2
            continue
        end
        plot!(fig, phi_arr, Sus[:,j], line=lines[m],
                c=colors[m], label=["$(name)", ""][(j!=1)+1],
                m=symbols[m], mc=:white, msc=colors[m], ma=0.8);
    end
end

# save results
display(fig);
png(fig, "$fig_path/Sus_compare.png");

## Counting for Sus error
@load "$exp_path/Sus_$(name_arr[1]).bson" Sus; Sus_0 = deepcopy(Sus); # true
for m in 2:length(mech_arr)
    name = name_arr[m];
    @load "$exp_path/Sus_$(name).bson" Sus;
    err = abs.(Sus .- Sus_0) ./ Sus_0 .* 100;
    @printf("Error of %s\n
             mean=%.1f%%, min=%.2f%%, max=%.2f%%, std=%.2f%%\n", name,
             mean(err), minimum(err), maximum(err), std(err));
end

## Counting for Sus costs
for m in 1:length(mech_arr)
    name = name_arr[m];
    @load "$exp_path/Sus_$(name).bson" t_costs;
    @printf("cost of %s=%.6f ± %.6f [s]\n", name, mean(t_costs), std(t_costs));
end
