def params_dict_to_str(params_dict):
    params = ['-{} {}'.format(k, v) for k, v in params_dict.items()]
    return ' '.join(params)
