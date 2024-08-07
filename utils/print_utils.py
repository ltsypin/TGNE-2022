# print_utils


def print_aligned_dict(d):
    key_lens = [len(k) for k in d.keys()]
    max_key_len = max(key_lens)

    for k, v in d.items():
        print(f'{" " * (max_key_len - len(k))}{k}:', v)