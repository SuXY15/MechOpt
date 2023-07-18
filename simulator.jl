# """
#     npr is the ratio of parameter numbers
#         size(p) = nr * npr
#     npr = 1 means only consider optimizing pre-exponential factor K
#     npr = 3 means will optimize K, b, Ea in equation ` ω = A T^b exp(-Ea/RT) `
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
#     Create ODE Problem instance for further integration
#     The defualt make_prob will only use species in conf["fuel"] and conf["oxyd"]
#     You can override this function in `datasets_\$fuelname.jl` for specific usage
#     `cond` represents conditions
# """
function make_prob(cond, p; tfinal=10.0)
    phi, Pj, T0 = cond[1:3];
    global P = Pj * one_atm;

    X0 = zeros(ns);
    for (f,v) in fuel
        X0[species_index(gas, f)] = v * phi
    end
    for (o,v) in oxyd
        X0[species_index(gas, o)] = v
    end
    X0 = X0 ./ sum(X0);

    Y0 = X2Y(gas, X0, dot(X0, gas.MW));
    u0 = vcat(Y0, T0);
    prob = ODEProblem(dudt!, u0, (0.0, tfinal), p);
    return prob
end

# """
#     This function gets the ignition delay time without interpolation
# """
function get_ind_ign(sol; dT=400)
    if maximum(sol[end, :]) - sol[end, 1] > dT
        return findfirst(sol[end, :] .> sol[end, 1] + dT)
    else
        println("Warning, no ignition")
        return length(sol.t)
    end
end

# """
#     Get profiles of IDT integration
# """
function get_Tcurve(cond, p; dT=400, dTabort=800, tfinal=10.0, saveat=[])
    phi, Pj, T0 = cond[1:3];
    prob = make_prob(cond, p; tfinal=tfinal)

    print("Tcurve, npr=$npr\r\n")
    condition(u, t, integrator) = u[end] > T0 + dTabort
    affect!(integrator) = terminate!(integrator)
    _cb = DiscreteCallback(condition, affect!)

    sol = solve(prob, CVODE_BDF(), saveat=saveat,
                reltol=1e-9, abstol=1e-12, callback=_cb);

    ts = sol.t
    pred = clamp.(Array(sol), 0, Inf);
    return ts, pred
end

# """
#     Get profiles of IDT integration, till IDT
# """
function get_Tcurve_idt(cond, p; dT=400, dTabort=800, tfinal=10.0, saveat=[])
    phi, Pj, T0 = cond[1:3];
    prob = make_prob(cond, p; tfinal=tfinal)

    condition(u, t, integrator) = u[end] > T0 + dTabort
    affect!(integrator) = terminate!(integrator)
    _cb = DiscreteCallback(condition, affect!)

    sol = solve(prob, CVODE_BDF(), saveat=saveat,
                reltol=1e-9, abstol=1e-12, callback=_cb);

    ind_ign = get_ind_ign(sol; dT=dT)
    ts = sol.t[1:ind_ign]
    pred = clamp.(Array(sol)[:, 1:ind_ign], 0, Inf);
    return ts, pred
end

# """
#     Interpolation for more accurate ignition times
#     by default, 10 points are used to fit the 1D curve
# """
function interpx(t, T, T_ign; N=10)
    pos = findfirst(T .> T_ign);
    _t = t[pos-N:pos];
    _T = T[pos-N:pos];
    f = Spline1D(_T, _t);
    t_ign = f(T_ign);
    return t_ign;
end

# """
#     Interpolation for more accurate IDTs
#     by default, 10 points are used to fit the 1D curve
# """
function get_idt(cond, p; dT=400, dTabort=800, tfinal=10.0)
    ts, pred = get_Tcurve(cond, p; dT=dT, dTabort=dTabort, tfinal=tfinal);
    idt = interpx(ts, pred[end,:], pred[end,1]+dT);
    return idt;
end

# """
#     Get Flame Speed via cantera
# """
function get_Su(cond, p)
    phi, Pj = cond[1:2];

    X0 = zeros(ns);
    for (f,v) in fuel
        X0[species_index(gas, f)] = v * phi
    end
    for (o,v) in oxyd
        X0[species_index(gas, o)] = v
    end
    X0 = X0 ./ sum(X0);

    ct_gas.set_multiplier(1.0)
    for (i,p_i) in enumerate(p)
        ct_gas.set_multiplier(exp(p_i), i-1)
    end
    ct_gas.TPX = 300.0, Pj*ct.one_atm, X0;

    width = 0.05;  # m
    f = ct.FreeFlame(ct_gas, width=width);
    f.set_refine_criteria(ratio=3, slope=0.07, curve=0.14);
    f.solve(loglevel=0, auto=true);
    Su = f.velocity[1];
    return Su
end
