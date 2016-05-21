import inflect


def percentage_fmt(ratio):
    return '{0:.1f}%'.format(100 * ratio)


def plural(noun, count):
    p = inflect.engine()
    return '{} {}'.format(count, p.plural(noun, count))


def seconds_to_hms_string(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)

    hms_string = ''
    if h > 0:
        hms_string += plural('hour', h) + ', '
    if m > 0:
        hms_string += plural('minute', m) + ', '
    hms_string += plural('second', s)

    return hms_string

