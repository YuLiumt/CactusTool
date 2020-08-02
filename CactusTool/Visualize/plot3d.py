def plot_trajectory(ax, x, y, z, color=None, label=None, plot_end=True, lw=1, nmasks=20):
    """
    trajectory with alpha
    """
    nmask_len = int(float(len(x))/nmasks)
    masks = []
    startdex, enddex = 0, nmask_len
    for k in range(nmasks):
        masks.append(range(startdex, enddex))
        startdex = enddex - 1
        enddex = enddex + nmask_len
        if k+1 == nmasks-1:
            enddex = len(x)
            
    alpha_min = 0.1
    alpha_max = 0.99
    for k, mask in enumerate(masks):
        alpha = alpha_min+k*(alpha_max-alpha_min)/(len(masks)-1)
        ax.plot(x[mask], y[mask], z[mask], color=color, alpha=alpha, lw=1)
    if plot_end:
        ax.scatter(x[-1], y[-1], z[-1], color=color, marker='o', s=10, label=label)