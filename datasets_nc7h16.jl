# """
#     master_mech: detailed one, used as ground truth to provide data for training
#     target_mech: skeletal or reduced one, used as initial model to be trained
# """
master_mech = mech_arr[1];
target_mech = mech_arr[2];

## Override settings
# """
#     npr = 1, which means only pre-exponential factor `K` is modified
# """
# global npr = 1;


#    the dudt function for solving ODE Problem, get du/dt from current status u, p, t
#    the symbol `!` means the variables are not changed inside the function
@inbounds function dudt!(du, u, p, t)
    Y = @view(u[1:ns])
    T = u[end]
    mean_MW = 1.0 / dot(Y, 1 ./ gas.MW)
    ρ_mass = P / R / T * mean_MW
    X = Y2X(gas, Y, mean_MW)
    C = Y2C(gas, Y, ρ_mass)
    cp_mole, cp_mass = get_cp(gas, T, X, mean_MW)
    h_mole = get_H(gas, T, Y, X)
    S0 = get_S(gas, T, P, X)

    # # for npr = 3
    # # The original Arrhenius equation:     ω =  A *      T^b  * exp(-Ea/RT)
    # # kp =  exp( p1 + p2 * log(T) - p3 / RT) = exp(p1) * T^p2 * exp(-p3/RT)
    # # so ω .* kp =  A * exp(p1) * T^(b+p2) *  exp(-(Ea+p3)/RT)
    # # which modifies the parameters to `A * exp(p1)`, `b+p2`, `Ea+p3`
    #
    if npr==3
        # print("in dudt! npr=$npr (3)\r\n")
        _p = reshape(p, nr, 3)
        kp = @. @views(exp(_p[:, 1] + _p[:, 2] * log(T) - _p[:, 3] * 4184.0 / R / T))
    end

    # for npr = 1
    # The original Arrhenius equation:     ω =  A *      T^b  * exp(-Ea/RT)
    # kp =  exp(p)
    # so ω .* kp =  A * exp(p) * T^b *  exp(-Ea/RT)
    # which modifies the parameters to `A * exp(p)`, (`b` and `Ea` not changed)
    if npr==1
        # print("in dudt! npr=$npr (1)\r\n")
        kp = exp.(p)
    end

    qdot = wdot_func(gas.reaction, T, C, S0, h_mole; get_qdot=true) .* kp
    wdot = gas.reaction.vk * qdot
    Ydot = wdot / ρ_mass .* gas.MW
    Tdot = -dot(h_mole, wdot) / ρ_mass / cp_mass
    du .= vcat(Ydot, Tdot)
    return du
end

# """
#     User override dudt!, for different npr setting
#     since dudt! is not changed compared from simulator.jl, it is not written here
# """

# """
#     User override make_prob, for different problem settings
# """
function make_prob(cond, p; tfinal=10.0)
    phi, Pj, T0 = cond[1:3];
    global P = Pj * one_atm;
    local fuel2air = 11;
    X0 = zeros(ns);
    X0[species_index(gas, "C7H16")] = phi
    X0[species_index(gas, "O2")] = fuel2air
    X0[species_index(gas, "N2")] = fuel2air * 3.76
    X0 = X0 ./ sum(X0);
    Y0 = X2Y(gas, X0, dot(X0, gas.MW));
    u0 = vcat(Y0, T0);

    prob = ODEProblem(dudt!, u0, (0.0, tfinal), p);
    return prob;
end

## generate IDT datasets
# """
#     If previous calculated conditions is avaliable, the program will load it
#     Otherwise the program will calculate conditions data and save it to file
#
#     For n-heptane, the conditions array `conds` stores data in 6 columns:
#     phi, P/one_atm, T0, true IDT, initial IDT, optimized IDT
# """
if is_restart && isfile("$exp_path/conds.bson")
    @load "$exp_path/conds.bson" conds;
    n_exp = size(conds, 1);
    n_train = Int64(n_exp * 0.8);
    println("Conds data loaded from $exp_path/conds.bson");
else
    Random.seed!(1);
    n_exp = Int64(conf["n_exp"]);
    plan = randomLHC(n_exp, 6) ./ n_exp;
    conds = scaleLHC(plan, [(Float64(conf["phi"]["lb"]),
                             Float64(conf["phi"]["ub"])),
                            (Float64(conf["pressure"]["lb"]),
                             Float64(conf["pressure"]["ub"])),
                            (Float64(conf["temperature"]["lb"]),
                             Float64(conf["temperature"]["ub"])),
                            (0.0, 1.0),
                            (0.0, 1.0),
                            (0.0, 1.0) ] );

    # calculate true IDTs by master mech
    global gas = CreateSolution(master_mech);
    global ns = gas.n_species;
    global nr = gas.n_reactions;
    p = zeros(nr * npr);
    for i_exp in 1:n_exp # this line could not be run under threads, since P is global
        cond = conds[i_exp, :];
        conds[i_exp, end-2] = get_idt(cond, p; dT=dT, tfinal=10.0);
        @printf("(%3d) phi %.1f P %4.1f [atm] T %6.1f [K] IDT0 = %.2e [s] \n",
                i_exp, cond[1], cond[2], cond[3], conds[i_exp, end-2] );
    end

    # calculate initial IDTs by target mech
    global gas = CreateSolution(target_mech);
    global ns = gas.n_species;
    global nr = gas.n_reactions;
    p = zeros(nr * npr);
    for i_exp in 1:n_exp # this line could not be run under threads, since P is global
        cond = conds[i_exp, :];
        conds[i_exp, end-1] = get_idt(cond, p; dT=dT, tfinal=10.0);
        @printf("(%3d) phi %.1f P %4.1f [atm] T0 %.1f[K] IDT0 = %.2e [s] IDT1 = %.2e [s] \n",
                i_exp, cond[1], cond[2], cond[3], conds[i_exp, end-2], conds[i_exp, end-1] );
    end

    @save "$exp_path/conds.bson" conds;
    n_train = Int64(n_exp * 0.8);
    println("Conds data calculated and saved in $exp_path/conds.bson");
end
