def percentage_fmt(ratio):
    return '{0:.1f}%'.format(100 * ratio)


def marker_from_sexcode(sexcode):
    if sexcode == 1:  # Male
        return 's'
    elif sexcode == 2:  # Female
        return 'o'

    return 'd'  # Unspecified


def hide_spines_and_ticks(ax, spines=["top", "right", "left"]):
    if spines == "all":
        spines = ["top", "bottom", "right", "left"]

    for spine in spines:
        ax.spines[spine].set_visible(False)

    ax.set_xticks([])
    ax.set_yticks([])


def seconds_to_hms_string(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)

    hms_string = ''

    if h > 0:
        if h == 1:
            hms_string += '{} hour, '.format(h)
        if h > 1:
            hms_string += '{} hours, '.format(h)
    if m > 0:
        if m == 1:
            hms_string += '{} minute, '.format(m)
        if m > 1:
            hms_string += '{} minutes, '.format(m)

    if s == 1:
        hms_string += '{} second'.format(s)
    else:
        hms_string += '{} seconds'.format(s)

    return hms_string
