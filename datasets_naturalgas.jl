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
global npr = 1;

# """
#     User override dudt, for different npr setting
# """
@inbounds function dudt!(du, u, p, t)
    Y = @view(u[1:ns]);
    T = u[end];
    mean_MW = 1.0 / dot(Y, 1 ./ gas.MW);
    ρ_mass = P / R / T * mean_MW;
    X = Y2X(gas, Y, mean_MW);
    C = Y2C(gas, Y, ρ_mass);
    cp_mole, cp_mass = get_cp(gas, T, X, mean_MW);
    h_mole = get_H(gas, T, Y, X);
    S0 = get_S(gas, T, P, X);

    # """ for npr = 3
    #     The original Arrhenius equation:     ω =  A *      T^b  * exp(-Ea/RT)
    #     kp =  exp( p1 + p2 * log(T) - p3 / RT) = exp(p1) * T^p2 * exp(-p3/RT)
    #     so ω .* kp =  A * exp(p1) * T^(b+p2) *  exp(-(Ea+p3)/RT)
    #     which modifies the parameters to `A * exp(p1)`, `b+p2`, `Ea+p3`
    # """

    # _p = reshape(p, nr, 3)
    # kp = @. @views(exp(_p[:, 1] + _p[:, 2] * log(T) - _p[:, 3] * 4184.0 / R / T))

    # """ for npr = 1
    #     The original Arrhenius equation:     ω =  A *      T^b  * exp(-Ea/RT)
    #     kp =  exp(p)
    #     so ω .* kp =  A * exp(p) * T^b *  exp(-Ea/RT)
    #     which modifies the parameters to `A * exp(p)`, (`b` and `Ea` not changed)
    # """
    kp = exp.(p);

    qdot = wdot_func(gas.reaction, T, C, S0, h_mole; get_qdot=true) .* kp;
    wdot = gas.reaction.vk * qdot;
    Ydot = wdot / ρ_mass .* gas.MW;
    Tdot = -dot(h_mole, wdot) / ρ_mass / cp_mass;
    du .= vcat(Ydot, Tdot);
    return du;
end

# """
#     User override make_prob, for different problem settings
# """
function make_prob(cond, p; tfinal=2000.0)
    phi, Pj, T0, aCH4, aC3H8 = cond[1:5];
    aC2H6 = 1 - aCH4 - aC3H8;
    global P = Pj * one_atm;
    local fuel2air = aCH4 * 2 + aC2H6 * 3.5 + aC3H8 * 5;
    X0 = zeros(ns);
    X0[species_index(gas, "CH4")] = phi * aCH4;
    X0[species_index(gas, "C2H6")] = phi * aC2H6;
    X0[species_index(gas, "C3H8")] = phi * aC3H8;
    X0[species_index(gas, "O2")] = fuel2air;
    X0[species_index(gas, "N2")] = fuel2air * 3.76;
    X0 = X0 ./ sum(X0);
    Y0 = X2Y(gas, X0, dot(X0, gas.MW));
    u0 = vcat(Y0, T0);

    prob = ODEProblem(dudt!, u0, (0.0, tfinal), p);
    return prob;
end

## For 1D Flame simulation
# """
#     Define global variables for flame speed optimization
#     cantera is used for calculating flame speed and flame speed sensitivity
#     For sensitiviy in native Julia, please refer Arrhenius_Flame_1D/src/flame_1d.jl
# """
if "optimize_flame_speed" in keys(conf) && conf["optimize_flame_speed"]
    using PyCall
    global ct = pyimport("cantera");
    global ct_gas = ct.Solution(master_mech);
end
function get_Su(cond, p)
    phi, Pj, aCH4, aC3H8 = cond[1:4];
    aC2H6 = 1 - aCH4 - aC3H8;
    local fuel2air = aCH4 * 2 + aC2H6 * 3.5 + aC3H8 * 5;
    X0 = zeros(ns);
    X0[species_index(gas, "CH4")] = phi * aCH4;
    X0[species_index(gas, "C2H6")] = phi * aC2H6;
    X0[species_index(gas, "C3H8")] = phi * aC3H8;
    X0[species_index(gas, "O2")] = fuel2air;
    X0[species_index(gas, "N2")] = fuel2air * 3.76;
    X0 = X0 ./ sum(X0);

    ct_gas.set_multiplier(1.0);
    for (i,p_i) in enumerate(p)
        ct_gas.set_multiplier(exp(p_i), i-1);
    end
    ct_gas.TPX = 300.0, Pj*ct.one_atm, X0;

    width = 0.05;  # m
    f = ct.FreeFlame(ct_gas, width=width);
    f.set_refine_criteria(ratio=3, slope=0.07, curve=0.14);
    f.solve(loglevel=0, auto=true);
    Su = f.velocity[1];
    sens = f.get_flame_speed_reaction_sensitivities();
    return Su, sens;
end

## generate IDT datasets
# """
#     If previous calculated conditions is avaliable, the program will load it
#     Otherwise the program will calculate conditions data and save it to file
#
#     For natural gas, the conditions array `conds` stores data in 8 columns:
#     phi, P/one_atm, T0, CH4 ratio, C3H8 ratio, true IDT, initial IDT, optimized IDT
# """
if is_restart && isfile("$exp_path/conds.bson")
    @load "$exp_path/conds.bson" conds;
    n_exp = size(conds, 1);
    n_train = Int64(n_exp * 0.8);
    println("Conds data loaded from $exp_path/conds.bson");
else
    Random.seed!(1);
    n_exp = Int64(conf["n_exp"]);
    plan = randomLHC(n_exp, 8) ./ n_exp;
    conds = scaleLHC(plan, [(Float64(conf["phi"]["lb"]),
                             Float64(conf["phi"]["ub"])),
                            (Float64(conf["pressure"]["lb"]),
                             Float64(conf["pressure"]["ub"])),
                            (Float64(conf["temperature"]["lb"]),
                             Float64(conf["temperature"]["ub"])),
                            (Float64(conf["aCH4"]["lb"]),
                             Float64(conf["aCH4"]["ub"])),
                            (Float64(conf["aC3H8"]["lb"]),
                             Float64(conf["aC3H8"]["ub"])),
                            (0.0, 1.0),
                            (0.0, 1.0),
                            (0.0, 1.0)]);
    # make sure the ratio sum does not exceed 1
    @threads for i_exp in 1:n_exp
        conds[i_exp,5] = minimum([conds[i_exp,5], 1-conds[i_exp,4]]);
    end

    # calculate true IDTs by master mech
    global gas = CreateSolution(master_mech);
    global ns = gas.n_species;
    global nr = gas.n_reactions;
    p = zeros(nr * npr);
    for i_exp in 1:n_exp # this line could not be run under threads, since P is global
        cond = conds[i_exp, :];
        conds[i_exp, end-2] = get_idt(cond, p; dT=dT, tfinal=10.0);
        @printf("(%3d) phi %.1f P %4.1f [atm] T %.1f [K] CH4 %.3f C3H8 %.3f IDT0 = %.2e [s] \n",
                i_exp, cond[1], cond[2], cond[3], cond[4], cond[5], conds[i_exp, end-2] );
    end

    # calculate initial IDTs by target mech
    global gas = CreateSolution(target_mech);
    global ns = gas.n_species;
    global nr = gas.n_reactions;
    p = zeros(nr * npr);
    for i_exp in 1:n_exp # this line could not be run under threads, since P is global
        cond = conds[i_exp, :];
        conds[i_exp, end-1] = get_idt(cond, p; dT=dT, tfinal=10.0);
        @printf("(%3d) phi %.1f P %4.1f [atm] T0 %.1f[K] CH4 %.3f C3H8 %.3f IDT0 = %.2e [s] IDT1 = %.2e [s] \n",
                i_exp, cond[1], cond[2], cond[3], cond[4], cond[5], conds[i_exp, end-2], conds[i_exp, end-1] );
    end

    conds[:,end] = conds[:,end-1];
    @save "$exp_path/conds.bson" conds
    n_train = Int64(n_exp * 0.8);
    println("Conds data calculated and saved in $exp_path/conds.bson");
end


## generate Su datasets
# """
#     If previous calculated conditions is avaliable, the program will load it
#     Otherwise the program will calculate conditions data and save it to file
#
#     For natural gas, the conditions array `conds_f` stores data in 7 columns:
#     phi, P/one_atm, CH4 ratio, C3H8 ratio, true Su, initial Su, optimized Su
# """
if "optimize_flame_speed" in keys(conf) && conf["optimize_flame_speed"]
    if is_restart && isfile("$exp_path/conds_f.bson")
        @load "$exp_path/conds_f.bson" conds_f;
        n_exp_f = size(conds_f, 1);
        n_train_f = Int64(n_exp_f * 0.8);
        println("Conds Flame data loaded from $exp_path/conds_f.bson");
    else
        Random.seed!(0x7777777);
        n_exp_f = Int64(conf["n_exp_f"]);
        plan = randomLHC(n_exp_f, 7) ./ n_exp_f;
        conds_f = scaleLHC(plan, [(Float64(conf["phi"]["lb"]),
                                   Float64(conf["phi"]["ub"])),
                                  (Float64(conf["pressure"]["lb"]),
                                   Float64(conf["pressure"]["ub"])),
                                  (Float64(conf["aCH4"]["lb"]),
                                   Float64(conf["aCH4"]["ub"])),
                                  (Float64(conf["aC3H8"]["lb"]),
                                   Float64(conf["aC3H8"]["ub"])),
                                  (0.0, 1.0),
                                  (0.0, 1.0),
                                  (0.0, 1.0)] );
        # make sure the ratio sum does not exceed 1
        @threads for i_exp in 1:n_exp_f
            conds_f[i_exp,4] = minimum([conds_f[i_exp,4], 1-conds_f[i_exp,3]]);
        end

        # calculate true Su by master mech
        global ct_gas = ct.Solution(master_mech);
        global gas = CreateSolution(master_mech);
        global ns = gas.n_species;
        global nr = gas.n_reactions;
        p = zeros(nr * npr);
        for i_exp in 1:n_exp_f # this line could not be run under threads, since P is global
            cond = conds_f[i_exp, :];
            conds_f[i_exp, end-2], _ = get_Su(cond, p);
            @printf("(%3d) phi %.1f P %4.1f [atm] CH4 %.3f C3H8 %.3f Su0 = %.2e [m/s] \n",
                    i_exp, cond[1], cond[2], cond[3], cond[4], conds_f[i_exp, end-2]);
        end

        # calculate initial IDTs by target mech
        global ct_gas = ct.Solution(target_mech);
        global gas = CreateSolution(target_mech);
        global ns = gas.n_species;
        global nr = gas.n_reactions;
        p = zeros(nr * npr);
        for i_exp in 1:n_exp_f # this line could not be run under threads, since P is global
            cond = conds_f[i_exp, :];
            conds_f[i_exp, end-1], _ = get_Su(cond, p);
            @printf("(%3d) phi %.1f P %4.1f [atm] CH4 %.3f C3H8 %.3f Su0 = %.2e [m/s] Su1 = %.2e [m/s] \n",
                    i_exp, cond[1], cond[2], cond[3], cond[4], conds_f[i_exp, end-2], conds_f[i_exp, end-1]);
        end
        conds_f[:,end] = conds_f[:,end-1];

        @save "$exp_path/conds_f.bson" conds_f;
        n_train_f = Int64(n_exp_f * 0.8);
        println("Conds Flame data calculated and saved in $exp_path/conds_f.bson");
    end
end
