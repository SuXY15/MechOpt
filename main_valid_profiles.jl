include("header.jl")
include("simulator.jl")
include("visual.jl")
using LaTeXStrings

Plots.resetfontsizes()

fm, fs = "Computer Modern", 14
# fm, fs = "sans-serif", 10
fnt = Plots.font(fs, fm)
default(titlefont=fnt, guidefont=fnt, tickfont=fnt,
        legendfont=fnt, legendfontsize=fs)

# ==============================================================================
# Validation profiles under single condition
phi = 1.0;
Pj = 1;
T0 = 1400;
tfinal = 10.0;
dTabort = 2000;
sp_names = ["OH"]
# sp_names = ["OH", "CH4", "H"]

# generate figure handles
lines = [:dot, :dash, :solid, :dashdot]
colors = [:black, :blue, :red, :green]
labels = ["Data", "SK34", "SK34OPa", "SK34OPb"]
widths = [5, 3, 3, 4]

figs = [
    plot(legend=false, grid=false, fontfamily=fm,
        foreground_color_legend=nothing, ylabel="T [K]"),
];
for sp in sp_names
    plt = plot(legend=:right, grid=false, fontfamily=fm, ylabel="Y($(sp))",
            foreground_color_legend=nothing, background_color_legend=nothing)
    push!(figs, plt)
end

# calculate and plot
idt = 0.
for (m,mech) in enumerate(mech_arr)
    global gas = CreateSolution(mech);
    global ns = gas.n_species;
    global nr = gas.n_reactions;
    cond = [phi, Pj, T0];
    @time ts, pred = get_Tcurve(cond, zeros(nr*npr); dT=dT, dTabort=dTabort, tfinal=tfinal);
    ts *= 1000
    if m==1
        idt = interpx(ts, pred[end,:], pred[end,1]+dT);
    end
    plot!(figs[1], ts, pred[end,:], msw=0, alpha=0.8, l=lines[m],
            c=colors[m], lw=widths[m], label=labels[m], xscale=:identity)
    for (i,sp) in enumerate(sp_names)
        plot!(figs[1+i], ts, pred[species_index(gas, sp),:], alpha=0.8,
                xlim=[ts[end]*0.0, ts[end]], msw=0, label=labels[m], l=lines[m], c=colors[m], lw=widths[m], xscale=:identity)
    end
    # compare_profiles(figs, ts, pred, sp_names, lines[m], colors[m], name_arr[m]);
end


# START external one with npr=3
m = 4
mech = "mechanism/nordin_34s121r_op3.yaml"
global gas = CreateSolution(mech);
@load string("results/nordin_sk34_op_npr=3", "/p.bson") p;
cond = [phi, Pj, T0];
@time ts, pred = get_Tcurve(cond, zeros(nr*npr); dT=dT, dTabort=dTabort, tfinal=tfinal);
ts *= 1000
if m==1
    idt = interpx(ts, pred[end,:], pred[end,1]+dT);
end
plot!(figs[1], ts, pred[end,:], msw=0, alpha=0.8, l=lines[m],
        c=colors[m], lw=widths[m], label=labels[m], xscale=:identity)
for (i,sp) in enumerate(sp_names)
    plot!(figs[1+i], ts, pred[species_index(gas, sp),:], alpha=0.8,
            xlim=[ts[end]*0.0, ts[end]], msw=0, label=labels[m], l=lines[m], c=colors[m], lw=widths[m], xscale=:identity)
end
@load string("results/nordin_sk34_op_npr=1", "/p.bson") p;
# FINISH external one with npr=3

# save results
for fig in figs
    # plot!(fig, xlim=[idt*0.5, idt*2])
    plot!(fig, xlim=[idt*0.8, idt*1.4])
    xlabel!(fig, "Time [ms]")
    # xticks!(fig, [0.8, 1.0, 1.2, 1.4])
end
# annotate!(figs[1], 1.05*idt, 1700, Plots.text(L"$P=40 \, atm $", :left))
# annotate!(figs[1], 1.05*idt, 1400, Plots.text(L"$T_0=700 \, K$", :left))
# annotate!(figs[1], 1.05*idt, 2100, Plots.text(L"$P=1 \, atm $", :left, fs))
# annotate!(figs[1], 1.05*idt, 1900, Plots.text(L"$T_0=1400 \, K$", :left, fs))
figs_high_T = figs;

# Low-temperature
# ==============================================================================
# Validation profiles under single condition
phi = 1.0;
Pj = 40;
T0 = 700;
tfinal = 10.0;
dTabort = 2000;

figs = [
    plot(legend=false, grid=false, fontfamily=fm,
        foreground_color_legend=nothing, ylabel="T [K]"),
];
for sp in sp_names
    plt = plot(legend=false, grid=false, fontfamily=fm, ylabel="Y($(sp))",
            foreground_color_legend=nothing)
    push!(figs, plt)
end

# calculate and plot
idt = 0.
for (m,mech) in enumerate(mech_arr)
    global gas = CreateSolution(mech);
    global ns = gas.n_species;
    global nr = gas.n_reactions;
    cond = [phi, Pj, T0];
    ts, pred = get_Tcurve(cond, zeros(nr*npr); dT=dT, dTabort=dTabort, tfinal=tfinal);
    ts *= 1000
    if m==1
        idt = interpx(ts, pred[end,:], pred[end,1]+dT);
    end
    plot!(figs[1], ts, pred[end,:], msw=0, alpha=0.8, l=lines[m],
            c=colors[m], lw=widths[m], label=labels[m], xscale=:identity)
    for (i,sp) in enumerate(sp_names)
        plot!(figs[1+i], ts, pred[species_index(gas, sp),:], alpha=0.8,
                xlim=[ts[end]*0.0, ts[end]], msw=0, label=labels[m], l=lines[m], c=colors[m], lw=widths[m], xscale=:identity)
    end
    # compare_profiles(figs, ts, pred, sp_names, lines[m], colors[m], name_arr[m]);
end

# START external one with npr=3
m=4
mech = "mechanism/nordin_34s121r_op3.yaml"
global gas = CreateSolution(mech);
@load string("results/nordin_sk34_op_npr=3", "/p.bson") p;
cond = [phi, Pj, T0];
@time ts, pred = get_Tcurve(cond, zeros(nr*npr); dT=dT, dTabort=dTabort, tfinal=tfinal);
ts *= 1000
if m==1
    idt = interpx(ts, pred[end,:], pred[end,1]+dT);
end
plot!(figs[1], ts, pred[end,:], msw=0, alpha=0.8, l=lines[m],
        c=colors[m], lw=widths[m], label=labels[m], xscale=:identity)
for (i,sp) in enumerate(sp_names)
    plot!(figs[1+i], ts, pred[species_index(gas, sp),:], alpha=0.8,
            xlim=[ts[end]*0.0, ts[end]], msw=0, label=labels[m], l=lines[m], c=colors[m], lw=widths[m], xscale=:identity)
end
@load string("results/nordin_sk34_op_npr=1", "/p.bson") p;
# FINISH external one with npr=3

# save results
for fig in figs
    # plot!(fig, xlim=[idt*0.5, idt*2])
    plot!(fig, xlim=[idt*0.8, idt*1.4])
    xlabel!(fig, "Time [ms]")
    # xticks!(fig, [0.8, 1.0, 1.2, 1.4])
end
# annotate!(figs[1], 1.05*idt, 1800, Plots.text(L"$P=40 \, atm $", :left))
# annotate!(figs[1], 1.05*idt, 1500, Plots.text(L"$T_0=700 \, K$", :left))
figs_low_T = figs;

plot!(figs_high_T[1], xticks=[], xlabel="", title=L"$P=1\, atm, T_0=1400\, K$")
plot!(figs_low_T[1], xticks=[], xlabel="",  title=L"$P=40\, atm, T_0=700\, K$")
plot!(figs_low_T[1], ylabel="")
plot!(figs_low_T[2], ylabel="")
using Measures
h = plot(figs_high_T[1], figs_low_T[1],
         figs_high_T[2], figs_low_T[2],
         size=(800,600), framestyle=:box)
display(h)
png(h, "$fig_path/profiles_compare_both.png");
