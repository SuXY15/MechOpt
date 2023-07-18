function plot_dots(d_true, d_init, d_pred; ylabel="η")
    h = plot(xlabel="Ref", ylabel=ylabel, legend=:topleft)
    plot!(d_true, d_true, l=:scatter, label="true");
    plot!(d_true, d_init, l=:scatter, label="init");
    plot!(d_true, d_pred, l=:scatter, label="pred");
    mi, ma = minimum(d_true), maximum(d_true);
    plot!([mi,ma], [mi,ma], l=:dash, c=:black, label="");
    display(h)
    png(h, "plot_$ylabel")
end

function plot_bars(d_true, d_init, d_pred; ylabel="UF")
    colors = palette(:default,16);
    x = 1:N;
    h = plot(xlabel="Index", ylabel=ylabel, legend=:bottomright)
    bar!(x.-0.25, d_true, c=colors[1], alpha=0.5, bar_width=0.2, label="true");
    bar!(x.+0.00, d_init, c=colors[2], alpha=0.5, bar_width=0.2, label="init");
    bar!(x.+0.25, d_pred, c=colors[3], alpha=0.5, bar_width=0.2, label="pred");
    display(h)
    png(h, "plot_bar_$ylabel")
end

function plot_pdf(θs_true, θs_init, θs_pred; i=1, j=2)
    function get_pdf(samples; i=1, j=2)
        mvnorm = fit(MvNormal, samples');
        x = range(minimum(samples[:,i]), maximum(samples[:,i]), length=20);
        y = range(minimum(samples[:,j]), maximum(samples[:,j]), length=20);
        z = [pdf(mvnorm, [xi,yj]) for xi in x, yj in y];
        return x, y, z;
    end
    h1 = plot()
    plot!(θs_true[:,i], θs_true[:,j], l=:scatter, msw=0, ms=3, alpha=0.5, label="true");
    # plot!(get_pdf(θs_true; i=i, j=j)...)

    plot!(θs_init[:,i], θs_init[:,j], l=:scatter, msw=0, ms=3, alpha=0.5, label="init");
    # plot!(get_pdf(θs_init; i=i, j=j)...)

    plot!(θs_pred[:,i], θs_pred[:,j], l=:scatter, msw=0, ms=3, alpha=0.5, label="pred");
    # plot!(get_pdf(θs_pred; i=i, j=j)...)

    h2 = plot(legend=false)
    histogram!(θs_true[:,j], alpha=0.3, normalize=true, orientation=:horizontal);
    histogram!(θs_init[:,j], alpha=0.3, normalize=true, orientation=:horizontal);
    histogram!(θs_pred[:,j], alpha=0.3, normalize=true, orientation=:horizontal);

    h3 = plot(legend=false)
    histogram!(θs_true[:,i], alpha=0.3, normalize=true);
    histogram!(θs_init[:,i], alpha=0.3, normalize=true);
    histogram!(θs_pred[:,i], alpha=0.3, normalize=true);

    h4 = plot(legend=false, axis=false, xticks=false, yticks=false, grid=false);

    h = plot(h1, h2, h3, h4, link=:y, layout=@layout[
        a{0.8w, 0.8h} a{0.2w, 0.8h}
        b{0.8w, 0.2h} b{0.2w, 0.2h}
    ]);
    display(h);
    png(h, "dist_x=$(i)_y=$(j)")
end
