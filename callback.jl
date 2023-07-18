l_min = 1.e-3; # minimum acceptable loss
l_loss = [];
l_loss_val = [];
l_grad = [];
l_pnorm = [];
iter = 1;
n_plot = 1;

# """
#     callback function, save results and plot figures
#     will be called during every epoch
# """
cb = function (p, loss, loss_val, g_norm; doplot=true)
    global l_loss, l_grad, iter
    push!(l_loss, loss);
    push!(l_loss_val, loss_val);
    push!(l_grad, g_norm);
    push!(l_pnorm, norm(p));

    if doplot && (iter % n_plot == 0)
        @save "$ckpt_path/mymodel.bson" p opt l_loss l_loss_val l_grad l_pnorm iter;
        plot_loss();

        @save "$exp_path/conds.bson" conds;
        regression_plot(;max=10);

        if "optimize_flame_speed" in keys(conf) && conf["optimize_flame_speed"]
            @save "$exp_path/conds_f.bson" conds_f;
            regression_f_plot();
        end

        if loss_val < l_min
            @save string(ckpt_path, "/model_$(length(l_loss)).bson") p opt l_loss l_loss_val l_grad l_pnorm iter;
            global l_min = loss_val;
        end
    end
    iter += 1;
    return false;
end

if is_restart && isfile("$ckpt_path/mymodel.bson")
    print("loading $ckpt_path/mymodel.bson\r\n")
    @load "$ckpt_path/mymodel.bson" p opt l_loss l_loss_val l_grad l_pnorm iter;
    iter += 1;
end

# """
#     Plot loss result
# """
function plot_loss()
    plt_loss = plot(xlabel="Epoch", ylabel="Loss");
    plot!(plt_loss, l_loss,     yscale=:log10, label="train");
    plot!(plt_loss, l_loss_val, yscale=:log10, label="val");

    plt_grad = plot(xlabel="Epoch", ylabel="Grad Norm");
    plot!(plt_grad, l_grad, yscale=:log10, label="grad_norm");

    plt_pnorm = plot(xlabel="Epoch", ylabel="p Norm");
    plot!(plt_pnorm, l_pnorm .+ 1.e-8, yscale=:identity, label="p_norm");

    plt_all = plot([plt_loss, plt_grad, plt_pnorm]..., legend=:top);
    png(plt_all, string(fig_path, "/loss_grad"));
end

# """
#     Load the pre-trained model of a specific epoch
# """
function load_model(epoch)
    @load string(ckpt_path, "/model_$(epoch).bson") p opt l_loss l_loss_val l_grad l_pnorm iter;
end
