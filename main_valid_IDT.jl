include("header.jl")
include("simulator.jl")
include("callback.jl")
include("visual.jl")

global npr = 1

fm, fs = "Computer Modern", 12
# fm, fs = "sans-serif", 10
fnt = Plots.font(fs, fm)
default(titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt, legendtitlefontsize=fs)

## prepare figure handles
figs = [];
for phi in phi_arr
    push!(figs, plot(legend=:topleft, yscale=:log10, framestyle=:box, size=(480, 360),
                     title="", foreground_color_legend=nothing, fontfamily=fm,
                     xlabel="1000/T [1/K]", ylabel="IDT [s]"));
end
lines = [:scatter, :dash, :solid, :dashdot, :dashdash, :dashdashdot, :solid];
colors = [:black, :blue, :red, :green, :magenta, :purple, :orange];
symbols = [:square, :circle, :diamond, :star, :pentagon, :octagon, :uptriangle];
labels = ["Data", "SK34", "SK34OPa", "SK34OPb"]
# labels = ["C7_sk42", "C7_sk34", "C7_sk34_opt"]

## get IDTs of each mechanism under each condition
for m in 1:length(mech_arr)
    name = name_arr[m];
    mech = mech_arr[m];

    if is_restart && isfile("$exp_path/IDTs_$name.bson")
        println("Loading IDTs for $name ...");
        @load "$exp_path/IDTs_$name.bson" IDTs;
    else
        println("Calculating IDTs for $name ...");

        global gas = CreateSolution(mech);
        global ns = gas.n_species;
        global nr = gas.n_reactions;

        # pre-run one time, to make sure the time cost is accurate
        t_cost = @elapsed idt = get_idt(
            [phi_arr[1], P_arr[1], T0_arr[1]], zeros(nr*npr); dT=dT, dTabort=dTabort, tfinal=100
        );

        IDTs = zeros(length(phi_arr), length(P_arr), length(T0_arr));
        t_costs = zeros(length(phi_arr), length(P_arr), length(T0_arr));
        for (k,phi) in enumerate(phi_arr)
            for (j,Pj) in enumerate(P_arr)
                for i in 1:length(T0_arr)
                    T0 = T0_arr[i];
                    t_cost = @elapsed idt = get_idt(
                        [phi, Pj, T0], zeros(nr*npr); dT=dT, dTabort=dTabort, tfinal=100
                    );
                    IDTs[k,j,i] = idt;
                    t_costs[k,j,i] = t_cost;
                    @printf("phi %.1f P %2d [atm] T %6.1f [K] IDT = %.2e [s] \n",
                             phi,     Pj,           T0,        idt );
                end
            end
        end
        @save "$exp_path/IDTs_$name.bson" IDTs t_costs;
    end

    println("Ploting IDTs for $name ... \n");
    for (k,phi) in enumerate(phi_arr)
        for (j,Pj) in enumerate(P_arr)
            if j==2
                continue
            end
            plot!(figs[k], 1000 ./ T0_arr, IDTs[k,j,:], line=lines[m],
                    c=colors[m], label=["$(labels[m])", ""][(j!=1)+1],
                    m=symbols[m], mc=:white, msc=colors[m], ma=0.8);
        end
    end
end

# START external one with npr=3
global npr = 3
m = 4
name = name_arr[3]
@load "results/nordin_sk34_op_npr=3/p.bson" p;
@load "results/nordin_sk34_op_npr=3/IDTs_$name.bson" IDTs t_costs;
IDTs_3 = deepcopy(IDTs)
for (k,phi) in enumerate(phi_arr)
    for (j,Pj) in enumerate(P_arr)
        if j==2
            continue
        end
        plot!(figs[k], 1000 ./ T0_arr, IDTs[k,j,:], line=lines[m],
                c=colors[m], label=["$(labels[m])", ""][(j!=1)+1],
                m=symbols[m], mc=:white, msc=colors[m], ma=0.8);
    end
end
err = abs.(IDTs .- IDTs_0) ./ IDTs_0 .* 100;
@printf("Error of %s npr=3
         mean=%.1f%%, min=%.2f%%, max=%.2f%%, std=%.2f%%\n", name,
         mean(err), minimum(err), maximum(err), std(err));
@load string("results/nordin_sk34_op_npr=1", "/p.bson") p;
global npr = 1
# FINISH external one with npr=3

annotate!(figs[2], 1.05, 8e-3, text("10 atm", :left, fs), fnt)
annotate!(figs[2], 0.95, 1e-4, text("60 atm", :left, fs), fnt)
# save results
display(figs[2])
for (k,phi) in enumerate(phi_arr)
    Plots.savefig(figs[k], "$fig_path/IDTs_phi=$phi.svg");
end

## Counting for IDTs error
@load "$exp_path/IDTs_$(name_arr[1]).bson" IDTs; IDTs_0 = deepcopy(IDTs); # true
for m in 2:length(mech_arr)
    name = name_arr[m];
    @load "$exp_path/IDTs_$(name).bson" IDTs;
    err = abs.(IDTs .- IDTs_0) ./ IDTs_0 .* 100;
    @printf("Error of %s
             mean=%.1f%%, min=%.2f%%, max=%.2f%%, std=%.2f%%\n", name,
             mean(err), minimum(err), maximum(err), std(err));
end

## Counting for IDTs costs
cost_compare_arr = Vector();
for m in 1:length(mech_arr)
    name = name_arr[m];
    mech = mech_arr[m];
    global gas = CreateSolution(mech);
    @load "$exp_path/IDTs_$(name).bson" t_costs;
    @printf("cost of %s (%ds %dr)
             = %.6f ± %.6f [s]\n", name, gas.n_species, gas.n_reactions,
             mean(t_costs), std(t_costs));
    push!(cost_compare_arr, [mean(t_costs), gas.n_species, gas.n_reactions])
end

## Compare costs against number of species and reactions
if length(mech_arr)>3
    cost_arr = vcat(cost_compare_arr'...);
    h = plot(xlabel="Time Cost [s]", ylabel="Number of Species/Reactions",
             xscale=:log10, yscale=:log10,
             legend=:topleft);
    plot!(h, cost_arr[:,1], cost_arr[:,2], l=:scatter, c=:red, label="Ns");
    plot!(h, cost_arr[:,1], cost_arr[:,3], l=:scatter, c=:blue, label="Nr");

    t_ratio = maximum(cost_arr[:,1]) / minimum(cost_arr[:,1]);
    s_ratio = maximum(cost_arr[:,2]) / minimum(cost_arr[:,2]);
    r_ratio = maximum(cost_arr[:,3]) / minimum(cost_arr[:,3]);

    s_slope = log(t_ratio)/log(s_ratio);
    r_slope = log(t_ratio)/log(r_ratio);

    cost_s = minimum(cost_arr[:,2]) .* (cost_arr[:,1] ./ minimum(cost_arr[:,1])) .^ (1/s_slope);
    cost_r = minimum(cost_arr[:,3]) .* (cost_arr[:,1] ./ minimum(cost_arr[:,1])) .^ (1/r_slope);
    plot!(h, cost_arr[:,1], cost_s, l=:dash, c=:red, label= @sprintf("slope=%.2f", s_slope));

    plot!(h, cost_arr[:,1], cost_r, l=:dash, c=:blue, label= @sprintf("slope=%.2f", r_slope));
    display(h)
end
