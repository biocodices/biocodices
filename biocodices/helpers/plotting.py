def hide_spines_and_ticks(ax, spines=["top", "right", "left"]):
    if spines == "all":
        spines = ["top", "bottom", "right", "left"]

    for spine in spines:
        ax.spines[spine].set_visible(False)

    ax.set_xticks([])
    ax.set_yticks([])


def marker_from_sexcode(sexcode):
    if sexcode == 1:  # Male
        return 's'
    elif sexcode == 2:  # Female
        return 'o'

    return 'd'  # Unspecified
