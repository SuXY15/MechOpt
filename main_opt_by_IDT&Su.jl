include("header.jl")
include("simulator.jl")
include("visual.jl")

global npr = 1

# include("datasets_naturalgas.jl")
include("datasets_nc7h16.jl")

## Preparations
if "optimize_flame_speed" in keys(conf) && conf["optimize_flame_speed"]
    global npr = 1;
    global ct_gas = ct.Solution(target_mech);
    flame_grad_multiplier = 10;
end

global gas = CreateSolution(target_mech);
global ns = gas.n_species;
global nr = gas.n_reactions;
p = zeros(nr*npr);

print("npr=$npr\r\n")

ind_sl = Int64.(conf["ind_sl"]);
include("sensitivity.jl");

## Training settings
n_epoch = 500;
p_weight = 1e-4;
batchsize = 40;
batchsize_val = 10;
batchsize_f = 10;
grad_max = 10^(2);
opt = ADAMW(Float64(conf["lr"]), (0.9, 0.999), Float64(conf["weight_decay"]));

## Training process
include("callback.jl")
epochs = ProgressBar(iter:n_epoch);
l_epoch = ones(batchsize+batchsize_val);
# grad_epoch = ones(nr*npr, batchsize)
grad_norm = ones(batchsize);
l_epoch_f = ones(batchsize_f);
grad_norm_f = ones(batchsize_f);
for epoch in epochs
    global p;
    l_epoch .= 1.0;
    grad_norm .= 1.0;

    # prepare samples for this epoch (size(epoch_samples) = batchsize)
    epoch_samples = sample(1:n_train, batchsize, replace=false);
    epoch_samples_val = n_train .+ sample(1:(n_exp-n_train), batchsize_val, replace=false);
    for (i,i_exp) in enumerate(epoch_samples)
        idt0 = conds[i_exp, end-2];
        ts, pred = get_Tcurve_idt(conds[i_exp,:], p; dT=dT, tfinal=10.0);
        idt = ts[end];
        conds[i_exp, end] = idt;
        l_epoch[i] = (log(idt / idt0))^2 + p_weight*norm(p, 2)^2;

        @printf("(E%3d:IDT%2d/%2d) i%3d nts: %d idt %.2e idt0 %.2e idt/idt0 %.3f\n",
                epoch, i, batchsize, i_exp, length(ts), idt, idt0, idt/idt0);

        grad = 2*log(idt / idt0) .* sensBVP_mthread(ts, pred, p) .+ 2*p_weight .* p;
        grad[ind_sl] .= 0.0;

        grad_norm[i] = norm(grad, 2);
        if grad_norm[i] > grad_max
            @. grad = grad / grad_norm[i] * grad_max;
        end
        if (iter > 1) # for calculating initial losses and grad_norms
            update!(opt, p, grad);
        end
        # grad_epoch[:,i] = deepcopy(grad)
    end

    # # mean grad is smoother to converge, but no other improvement
    # grad_mean = mean(grad_epoch, dims=2)[:]
    # if (iter > 1) # for calculating initial losses and grad_norms
    #     update!(opt, p, grad_mean);
    # end

    loss = mean(l_epoch[1:batchsize]);
    g_norm = mean(grad_norm[1:batchsize]);

    # calculate loss in validation sets
    for (i,i_exp) in enumerate(epoch_samples_val)
        idt0 = conds[i_exp, end-2];
        ts, pred = get_Tcurve_idt(conds[i_exp,:], p; dT=dT, tfinal=10.0);
        idt = ts[end];
        conds[i_exp, end] = idt;
        l_epoch[batchsize+i] = (log(idt / idt0))^2 + p_weight*norm(p, 2)^2;
    end
    loss_val = mean(l_epoch[batchsize+1:end]);

    # training for Su
    if "optimize_flame_speed" in keys(conf) && conf["optimize_flame_speed"]
        epoch_samples = sample(1:n_train_f, batchsize_f, replace=false);
        for (i,i_exp) in enumerate(epoch_samples)
            Su0 = conds_f[i_exp, end-2];
            Su, sens = get_Su(conds_f[i_exp,:], p);
            conds_f[i_exp, end] = Su;
            l_epoch_f[i] = (log(Su / Su0))^2 + p_weight*norm(p, 2)^2;

            @printf("(E%3d:Su %2d/%2d) i%4d Su %.2e Su0 %.2e Su/Su0 %.3f\n",
                    epoch, i, batchsize_f, i_exp, Su, Su0, Su/Su0);

            grad = 2 * flame_grad_multiplier * log(Su / Su0) .* sens .+ 2*p_weight .* p;
            grad[ind_sl] .= 0.0;

            grad_norm_f[i] = norm(grad, 2);
            if grad_norm_f[i] > grad_max
                @. grad = grad / grad_norm_f[i] * grad_max;
            end

            if (iter > 1) # for calculating initial losses and grad_norms
                update!(opt, p, grad);
            end
        end
        loss_f = mean(l_epoch_f);
        g_norm_f = mean(grad_norm_f);
        loss_val_f = mean((log.(conds_f[:,end]./conds_f[:,end-2])).^2) + p_weight*norm(p, 2)^2;

        loss = mean([loss, loss_f]);
        g_norm = mean([g_norm, g_norm_f]);
        loss_val = mean([loss_val, loss_val_f]);
    end

    set_description( epochs,
        string(
            @sprintf("loss %.2e val %.2e pnorm %.2f gnorm %.2f",
                      loss,   loss_val,  norm(p),  g_norm );
        )
    )
    cb(p, loss, loss_val, g_norm; doplot=true);
end

## post-processing
validation_plot(; max=50, name="_optimized");
regression_plot(; max=50, name="_optimized");
if "optimize_flame_speed" in keys(conf) && conf["optimize_flame_speed"]
    regression_f_plot(; max=50, name="_optimized");
end

@save string(exp_path, "/p.bson") p;
# @load string(exp_path, "/p.bson") p;

updateyaml(target_mech, p; name="_op");
