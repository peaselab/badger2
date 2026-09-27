import sys
import os
from time import time as time_now
import math
import itertools
import numpy as np
import argparse
from scipy.stats import ttest_ind


def generate_argparse(arguments=None):
    parser = argparse.ArgumentParser(
        description=(
            'Bivariate Analysis of Differential Gene Expression Reactions')
    )
    parser.add_argument(
        '--contrasts-a',
        required=True,
        help=("factor A contrasts with factors in a pair "
              "separated by dashes in order positive-negative "
              "and the factor pairs separate by commas, "
              "(example: \"LABEL1-LABEL2, LABEL2-LABEL\")")
    )
    parser.add_argument(
        '--contrasts-b',
        required=True,
        help=("factor B contrasts with factors in a pair "
              "separated by dashes in order positive-negative "
              "and the factor pairs separate by commas, "
              "(example: \"LABEL1-LABEL2, LABEL2-LABEL\")")
    )
    parser.add_argument(
        '--metadata',
        required=True,
        help=("CSV file of experimental design factors "
              "with one sample per row "
              "and first line has factor names as headers")
    )
    parser.add_argument(
        '--input',
        required=True,
        help='CSV file of normalized read counts'
    )
    parser.add_argument(
        '--transform',
        default="cpm",
        choices=('none', 'cpm', 'deseq'),
        required=False,
        help=("transform the data to counts-per-million (cpm) "
              "or using the DESeq method "
              "(scaled pseudocounts to make all libraries "
              "the average library size)")
    )
    parser.add_argument(
        '--min-avg', default=0, type=float,
        help=("minimum average transformed value per row. "
              "Rows not meeting this value are filetered.")
    )
    parser.add_argument(
        '--pair',
        help=("Optional: For paired samples over factor A, "
              "use this metadata field to pair samples. "
              "For unpaired analysis, do not use this flag.")
    )
    parser.add_argument(
        '--factorname-a',
        default="Time",
        help="metadata header name of the factor A variable"
    )
    parser.add_argument(
        '--factorname-b',
        default="Treatment",
        help="metadata header name of the factor B variable"
    )
    parser.add_argument(
        '--max-perms',
        type=int,
        default=10000,
        help="If the number of permutations is large, "
        "this variable caps the number of permutations. "
        "Generally, this only affects unpaired analysis."
    )
    parser.add_argument(
        '--correction',
        choices=("none", "Bonferroni", "BH"),
        default=("none"),
        help=("Multiple testing correction for P-values. "
              "Default is none.")
    )
    parser.add_argument(
        '--seed',
        type=int,
        help="Seed for the randomization of the subset of"
        "permutations when the value grows large. Leave blank for random.")
    return parser.parse_args(arguments)


def transform_filter_data(
        data_filepath, mode='cpm', minimum_average=0, delim=','):
    with open(data_filepath) as input_file:
        headers = input_file.readline().rstrip().split(delim)
        lib_sizes = dict.fromkeys(headers[1:], 0)
        for line in input_file:
            line = line.rstrip().split(delim)
            for i in range(1, len(line)):
                lib_sizes[headers[i]] += float(line[i])
    mean_lib_size = np.mean(list(lib_sizes.values()))
    with open(data_filepath) as input_file, open(
            data_filepath + "." + mode, 'w') as transform_file:
        transform_file.write(input_file.readline())
        for line in input_file:
            line = line.rstrip().split(delim)
            gene = line[0]

            if mode == 'cpm':
                transform_data = [
                    (float(x)) / lib_sizes[headers[i + 1]] * 1000000
                    if lib_sizes[headers[i + 1]] > 0 else 0
                    for i, x in enumerate(line[1:])
                ]

            elif mode == 'deseq':
                transform_data = [
                    round(float(x) / lib_sizes[headers[i + 1]] * mean_lib_size)
                    if lib_sizes[headers[i + 1]] > 0
                    else 0
                    for i, x in enumerate(line[1:])
                ]

            elif mode == 'none':
                transform_data = [
                    float(x) for x in line[1:]
                ]

            if np.mean(transform_data) >= minimum_average:
                transform_file.write(
                    delim.join([gene] + [str(x) for x in transform_data]) +
                    "\n")
    return ''


def parse_factors(file):
    data_dict = {}
    with open(file) as design_file:
        header = design_file.readline().rstrip().split(",")
        for line in design_file:
            entry = line.rstrip().split(",")
            if len(entry) != len(header):
                raise RuntimeError("Metadata entry not same length as header")
            data_dict[entry[0]] = dict([
                (header[i], entry[i]) for i in range(1, len(header))])
    return data_dict


def parse_comma_dash(time_points):
    return [
        tuple(part.strip() for part in x.split("-"))
        for x in time_points.split(",")
        ]


def zero_div(a, b):
    return a / b if b != 0 else 0


def invert_bool(values):
    return [not value for value in values]


def parse_ids(fa, fb, metadata, args):
    xa = args.factorname_a
    xb = args.factorname_b
    if args.pair:
        ret = []
        pair_sets = []
        for (x, y) in ((0, 0), (1, 0), (0, 1), (1, 1)):
            samples = sorted((
                    k for k, v in metadata.items()
                    if (v[xa], v[xb]) == (fa[x], fb[y])
                ),
                key=lambda k: metadata[k][args.pair]
            )
            ret.append(samples)
            pair_sets.append(
                [
                    metadata[s][args.pair] for s in samples
                    ]
            )
        if any(len(s) != len(ret[0]) for s in ret):
            raise RuntimeError("paired group unequal numbers")

        if (pair_sets[0] != pair_sets[1]) or (pair_sets[2] != pair_sets[3]):
            raise RuntimeError(f"pairs not matching: {pair_sets}")

        return ret
    ret = []
    for (x, y) in ((0, 0), (1, 0), (0, 1), (1, 1)):
        s = [k for k, v in metadata.items()
             if (v[xa], v[xb]) == (fa[x], fb[y])]
        ret.append(s)
        print(s, fa[x], fb[y])
    return ret


def benjamini_hochberg(p_values):
    p_values = np.asarray(p_values, dtype=float)
    n = len(p_values)
    order = np.argsort(p_values)
    ranked_p = p_values[order]
    adjusted = ranked_p * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    result = np.empty(n)
    result[order] = adjusted
    return result


def build_output_path(input_file, a, b):
    return (f"badger_output_{os.path.basename(input_file)}"
            f"_{a[0]}v{a[1]}.{b[0]}v{b[1]}.csv")


def build_transform_path(outfile_path, transform):
    return (f"{outfile_path}.{transform}")


def validate_id_sets(headers, samples):
    for xid in samples:
        if xid not in headers:
            raise RuntimeError(f"{xid} is missing from headers")
    return ''


def build_header_index(headers, samples):
    return np.array([
        headers.index(idx) for idx in samples
        ], dtype=int)


def generate_header(fa, fb):
    return [
        "Gene",
        "NPermut",
        "NSig",
        "P",
        "Padj",
        "CritPVal",
        "MaxExpr",
        "SDD",
        f"MeanD_{fb[0]}-MeanD_{fb[1]}",
        f"MeanD_{fb[1]}",
        f"MeanD_{fb[0]}",
        f"Mean_{fb[1]}_{fa[1]}",
        f"Mean_{fb[1]}_{fa[0]}",
        f"Mean_{fb[0]}_{fa[1]}",
        f"Mean_{fb[0]}_{fa[0]}",
        "MaxPerm",
        "MinPerm",
        "MeanPerm",
    ]


def adaptive_p_format(value):
    if value > 0.001:
        return f"{value:.5f}"
    return f"{value:.3e}"


def format_entry(entry):
    return ",".join([
        entry['gene_id'],
        f"{entry['n_perms']}",
        f"{entry['n_sig']}",
        adaptive_p_format(entry['p']),
        adaptive_p_format(entry.get('p_adj', entry['p'])),
        ("NA" if entry['crit_p'] == "NA" else
         adaptive_p_format(entry['crit_p'])
         ),
        f"{zero_div(entry['dd'], entry['max_expr']):.8f}",
        f"{entry['max_expr']:.8f}",
        f"{round(entry['dd'], 8):.8f}",
        f"{round(entry['means'][1],6)}",
        f"{round(entry['means'][0],6)}",
        f"{round(entry['group_means'][0],6)}",
        f"{round(entry['group_means'][1],6)}",
        f"{round(entry['group_means'][2],6)}",
        f"{round(entry['group_means'][3],6)}",
        f"{round(entry['mean_dd_perms'], 8):.8f}",
        f"{round(entry['max_dd_perms'], 8):.8f}",
        f"{round(entry['min_dd_perms'], 8):.8f}",
    ])


def multiple_test_correction(entries, correction="none"):
    new_entries = []
    if correction == "BH":
        adjusted_p_values = benjamini_hochberg([x['p'] for x in entries])
    for i, entry in enumerate(entries):
        new_entry = entry.copy()
        if correction == 'Bonferroni':
            new_entry['p_adj'] = min(entry['p'] * len(entries), 1.0)
        elif correction == "BH":
            new_entry['p_adj'] = adjusted_p_values[i]
        new_entries.append(new_entry)
    return new_entries


def report_progress(n_block, time_prev):
    time_rep = time_now()
    time_x = (
        (time_rep - time_prev)/n_block
        if n_block > 0 else 0
    )
    print(f"{n_block} Entries processed in "
          f"{time_rep - time_prev:.2f} seconds."
          f" ({time_x:.5f})")
    return time_rep


# ###################################
# UNPAIRED MODE
# ###################################


def parse_ids_unpaired(fa, fb, metadata, args):
    xa = args.factorname_a
    xb = args.factorname_b
    treatment_time1 = tuple(
        k for k, v in metadata.items()
        if v[xa] == fa[0] and v[xb] == fb[0]
    )
    treatment_time0 = tuple(
        k for k, v in metadata.items()
        if v[xa] == fa[1] and v[xb] == fb[0]
    )
    control_time1 = tuple(
        k for k, v in metadata.items()
        if v[xa] == fa[0] and v[xb] == fb[1]
    )
    control_time0 = tuple(
        k for k, v in metadata.items()
        if v[xa] == fa[1] and v[xb] == fb[1]
    )
    samples = (
        treatment_time1 +
        treatment_time0 +
        control_time1 +
        control_time0
    )
    mask_b = np.array(
        [True] * len(treatment_time1) +
        [True] * len(treatment_time0) +
        [False] * len(control_time1) +
        [False] * len(control_time0)
    )
    mask_a = np.array(
        [True] * len(treatment_time1) +
        [False] * len(treatment_time0) +
        [True] * len(control_time1) +
        [False] * len(control_time0)
    )
    return samples, mask_a, mask_b


def make_permutations_unpaired(sample_indices, mask_a, mask_b,
                               max_permutations=100000, seed=None):
    rng = np.random.default_rng(seed)
    groups = [
        np.asarray(sample_indices[mask_b & mask_a] - 1, dtype=int),
        np.asarray(sample_indices[mask_b & ~mask_a] - 1, dtype=int),
        np.asarray(sample_indices[~mask_b & mask_a] - 1, dtype=int),
        np.asarray(sample_indices[~mask_b & ~mask_a] - 1, dtype=int),
    ]
    group_sizes = [len(group) for group in groups]
    block_b1 = np.concatenate([groups[0], groups[1]])
    block_b2 = np.concatenate([groups[2], groups[3]])
    n_a1_b1 = group_sizes[0]
    n_a1_b2 = group_sizes[2]

    original = tuple(tuple(sorted(group))
                     for group in groups)
    seen = set()
    perms = []
    perms.append(np.concatenate([
        np.asarray(original[0], dtype=int),
        np.asarray(original[1], dtype=int),
        np.asarray(original[2], dtype=int),
        np.asarray(original[3], dtype=int),
    ]))
    max_possible = (
        math.comb(len(block_b1), n_a1_b1) *
        math.comb(len(block_b2), n_a1_b2)
    )

    if max_possible < 1:
        raise RuntimeError(
            "Insufficient replicates for permutation test"
        )

    if max_possible <= max_permutations:
        print(f"Using all {max_possible} possible permutations.")

        for a1_b1 in itertools.combinations(
                block_b1, n_a1_b1):

            a1_b1 = tuple(sorted(a1_b1))
            a2_b1 = tuple(sorted(
                set(block_b1) - set(a1_b1)
            ))

            for a1_b2 in itertools.combinations(
                    block_b2, n_a1_b2):

                a1_b2 = tuple(sorted(a1_b2))
                a2_b2 = tuple(sorted(
                    set(block_b2) - set(a1_b2)
                ))

                key = (
                    a1_b1,
                    a2_b1,
                    a1_b2,
                    a2_b2,
                )

                if key == original:
                    continue

                perms.append(np.concatenate([
                    np.asarray(a1_b1, dtype=int),
                    np.asarray(a2_b1, dtype=int),
                    np.asarray(a1_b2, dtype=int),
                    np.asarray(a2_b2, dtype=int),
                ]))

    else:
        print(
            f"{max_possible} possible permutations exceeds cutoff, "
            f"using {max_permutations} random permutations."
        )

        seen = {original}

        while len(perms) < max_permutations + 1:
            shuffled_b1 = rng.permutation(block_b1)
            shuffled_b2 = rng.permutation(block_b2)

            new_0 = tuple(
                sorted(shuffled_b1[:n_a1_b1])
            )
            new_1 = tuple(
                sorted(shuffled_b1[n_a1_b1:])
            )
            new_2 = tuple(
                sorted(shuffled_b2[:n_a1_b2])
            )
            new_3 = tuple(
                sorted(shuffled_b2[n_a1_b2:])
            )
            key = (new_0, new_1, new_2, new_3)
            if key in seen:
                continue
            seen.add(key)

            perms.append(np.concatenate([
                np.asarray(new_0, dtype=int),
                np.asarray(new_1, dtype=int),
                np.asarray(new_2, dtype=int),
                np.asarray(new_3, dtype=int),
            ]))
    perms = np.asarray(perms, dtype=int)
    print(f"Using {len(perms)} unpaired permutations.")
    return perms, group_sizes


def analyze_gene_unpaired(values, perms, mask_a, mask_b):
    vals_t1 = values[perms[:, mask_b & mask_a]]
    vals_t0 = values[perms[:, mask_b & ~mask_a]]
    vals_c1 = values[perms[:, ~mask_b & mask_a]]
    vals_c0 = values[perms[:, ~mask_b & ~mask_a]]
    diffs_t = np.mean(vals_t1, axis=1) - np.mean(vals_t0, axis=1)
    diffs_c = np.mean(vals_c1, axis=1) - np.mean(vals_c0, axis=1)
    dd_values = diffs_t - diffs_c
    if dd_values[0] > 0:
        n_sig = np.count_nonzero(dd_values[1:] > dd_values[0])
    else:
        n_sig = np.count_nonzero(dd_values[1:] < dd_values[0])
    return {
        "means": (
            diffs_t[0].mean(),
            diffs_c[0].mean()
        ),
        "group_means": (
            vals_t1[0].mean(),
            vals_t0[0].mean(),
            vals_c1[0].mean(),
            vals_c0[0].mean(),
        ),
        "n_perms": len(perms) - 1,
        "n_sig": n_sig,
        'p': (n_sig) / (len(perms)),
        "dd": dd_values[0],
        "crit_p": 'NA',
        "max_expr": max(vals_t1[0].max(), vals_t0[0].max(),
                        vals_c1[0].max(), vals_c0[0].max()),
        "max_dd_perms": dd_values[1:].max(),
        "min_dd_perms": dd_values[1:].min(),
        "mean_dd_perms": dd_values[1:].mean(),
    }


def run_contrast_unpaired(a_factor, b_factor, metadata, args):
    print(f"Running ==== A: {a_factor} x B: {b_factor}")
    samples, mask_a, mask_b = parse_ids_unpaired(
        a_factor,
        b_factor,
        metadata,
        args
    )
    print(samples)
    out_file_path = build_output_path(
        args.input,
        a_factor,
        b_factor
    )
    transform_file_path = build_transform_path(
        args.input,
        args.transform
    )
    time_prev = time_now()
    n_line = 0
    entries = []
    n_block = 5000
    with open(transform_file_path) as datafile:

        headers = datafile.readline().rstrip().split(",")

        validate_id_sets(headers, samples)

        sample_indices = build_header_index(
            headers,
            samples
        )

        perms, group_sizes = make_permutations_unpaired(
            sample_indices=sample_indices,
            mask_a=mask_a,
            mask_b=mask_b,
            max_permutations=args.max_perms,
            seed=args.seed
        )

        n_perms = len(perms) - 1
        print("total permutations sampled: ", n_perms)
        if len(perms) == 0:
            raise RuntimeError("No valid permutations available")
        print("permutations established")

        for line in datafile:
            gene_id, _, remainder = line.partition(",")
            values = np.fromstring(remainder, sep=",")
            gene_stats = analyze_gene_unpaired(
                values,
                perms,
                mask_a,
                mask_b
            )

            gene_stats['gene_id'] = gene_id

            entries.append(
                gene_stats
            )

            n_line += 1

            if n_line % n_block == 0:
                time_prev = report_progress(
                    n_block,
                    time_prev
                )

    with open(out_file_path, "w") as outfile:
        outfile.write(
            ",".join(generate_header(a_factor, b_factor)) + "\n"
        )
        if args.correction != 'none':
            entries = multiple_test_correction(
                entries,
                correction=args.correction
            )
            print(f"Multiple Testing Correction: {args.correction}")
        outfile.write(
            '\n'.join(
                format_entry(
                    entry)
                for entry in entries
            )
        )
    return ''


# ###################################
# PAIRED MODE
# ###################################


def parse_ids_paired(fa, fb, metadata, args):
    xa = args.factorname_a
    xb = args.factorname_b
    ret = []
    pair_sets = []
    for (x, y) in ((1, 1), (0, 1), (1, 0), (0, 0)):
        samples = sorted(
            (
                k for k, v in metadata.items()
                if (v[xa], v[xb]) == (fa[x], fb[y])
            ),
            key=lambda k: metadata[k][args.pair]
        )
        ret.append(samples)
        pair_sets.append([metadata[s][args.pair] for s in samples])

    if any(len(s) != len(ret[0]) for s in ret):
        raise RuntimeError("paired group unequal numbers")

    if (pair_sets[0] != pair_sets[1]) or (pair_sets[2] != pair_sets[3]):
        raise RuntimeError(f"pairs not matching: {pair_sets}")

    pair_ids_t = [
        pair_id
        for p1, p0 in zip(pair_sets[0], pair_sets[1])
        for pair_id in (p1, p0)
    ]

    pair_ids_c = [
        pair_id
        for p1, p0 in zip(pair_sets[2], pair_sets[3])
        for pair_id in (p1, p0)
    ]

    pair_ids = pair_ids_t + pair_ids_c

    pair_map = {
        value: i
        for i, value in enumerate(dict.fromkeys(pair_ids))
    }

    mask_paired = np.array(
        [pair_map[value] for value in pair_ids],
        dtype=int
    )

    samples_t = [
        sample
        for pair_samples in zip(ret[0], ret[1])
        for sample in pair_samples
    ]

    samples_c = [
        sample
        for pair_samples in zip(ret[2], ret[3])
        for sample in pair_samples
    ]

    samples = samples_t + samples_c
    mask_ab = np.array(
        [True] * len(samples_t) +
        [False] * len(samples_c),
        dtype=bool
    )
    return samples, mask_ab, mask_paired


def make_permutations_paired(
        sample_indices,
        mask_ab,
        mask_pairs,
        max_permutations=100000):

    original_t = sample_indices[mask_ab]
    original_c = sample_indices[~mask_ab]

    pairs_t = mask_pairs[mask_ab]
    pairs_c = mask_pairs[~mask_ab]
    indices_t = [
        tuple(x for i, x in enumerate(original_t)
              if pairs_t[i] == xt)
        for xt in set(pairs_t)
    ]

    indices_c = [
        tuple(x for i, x in enumerate(original_c)
              if pairs_c[i] == xc)
        for xc in set(pairs_c)
    ]

    fs_indices_t = frozenset(indices_t)
    fs_indices_c = frozenset(indices_c)

    size_t = len(indices_t)
    size_c = len(indices_c)

    indices_all = indices_t + indices_c
    n_perms = 0

    original = [
        item for pair in list(sorted(fs_indices_t)) for item in pair
        ] + [
        item for pair in list(sorted(fs_indices_c)) for item in pair
        ]

    perms = [original]

    if size_t == size_c:
        anchor = indices_all[0]
        remaining = indices_all[1:]
        for others in itertools.combinations(remaining, size_c - 1):
            new_t = frozenset((anchor,) + others)
            new_c = frozenset(indices_all) - new_t

            if new_c == fs_indices_c and new_t == fs_indices_t:
                continue

            perms.append(
                [
                    item for pair in list(sorted(new_t)) for item in pair
                ] + [
                    item for pair in list(sorted(new_c)) for item in pair
                ]
            )

            if len(perms) >= max_permutations + 1:
                print(
                    f"Number of possible permutations exceeds cutoff, "
                    f"using first {max_permutations}."
                )
                break

    else:
        for comb_t in itertools.combinations(indices_all, size_t):
            new_t = frozenset(comb_t)
            new_c = frozenset(indices_all) - new_t

            if new_t == fs_indices_t and new_c == fs_indices_c:
                continue

            if new_c == fs_indices_t and new_t == fs_indices_c:
                continue

            perms.append(
                [
                    item
                    for pair in sorted(new_t)
                    for item in pair
                ] + [
                    item
                    for pair in sorted(new_c)
                    for item in pair
                ]
            )
            n_perms += 1
            if n_perms >= max_permutations:
                print(
                    f"Number of possible permutations exceeds cutoff, "
                    f"using first {max_permutations}."
                )
                break

    print(f"Using {len(perms) - 1} paired permutations.")
    return np.array(perms, dtype=int)


def analyze_gene_paired(values, perms, mask_ab):
    vals_t = values[perms[:, mask_ab] - 1]
    vals_c = values[perms[:, ~mask_ab] - 1]

    diffs_t = (vals_t[:, ::2] - vals_t[:, 1::2])
    diffs_c = (vals_c[:, ::2] - vals_c[:, 1::2])

    dd_values = (
        np.mean(diffs_t, axis=1) -
        np.mean(diffs_c, axis=1)
    )

    _, p_vals = ttest_ind(
        diffs_t,
        diffs_c,
        axis=1
    )

    n_sig = np.count_nonzero(p_vals[1:] < p_vals[0])
    permuted_dd = dd_values[1:]
    return {
        "means": (
            diffs_t[0].mean(),
            diffs_c[0].mean()
        ),
        "group_means": (
            vals_t[0, ::2].mean(),
            vals_t[:, 1::2].mean(),
            vals_c[0, ::2].mean(),
            vals_c[:, 1::2].mean(),
        ),
        "n_perms": len(perms) - 1,
        "n_sig": n_sig,
        'p': (n_sig) / (len(perms)),
        "dd": dd_values[0],
        "crit_p": p_vals[0],
        "max_expr": max(
            vals_t[0].max(),
            vals_c[0].max()
        ),
        "max_dd_perms": permuted_dd.max() if len(permuted_dd) else np.nan,
        "min_dd_perms": permuted_dd.min() if len(permuted_dd) else np.nan,
        "mean_dd_perms": permuted_dd.mean() if len(permuted_dd) else np.nan,
    }


def run_contrast_paired(a_factor, b_factor, metadata, args):
    print(f"Running ==== A: {a_factor} x B: {b_factor}")
    samples, mask_ab, mask_pairs = parse_ids_paired(
        a_factor,
        b_factor,
        metadata,
        args
    )
    out_file_path = build_output_path(
        args.input,
        a_factor,
        b_factor
    )
    transform_file_path = build_transform_path(
        args.input,
        args.transform
    )
    time_prev = time_now()
    n_line = 0
    entries = []
    n_block = 5000
    with open(transform_file_path) as datafile:
        headers = datafile.readline().rstrip().split(",")

        validate_id_sets(headers, samples)
        sample_indices = build_header_index(
            headers,
            samples
        )

        perms = make_permutations_paired(
            sample_indices=sample_indices,
            mask_ab=mask_ab,
            mask_pairs=mask_pairs,
            max_permutations=args.max_perms
        )

        n_perms = len(perms)
        print("total permutations: ", n_perms - 1)
        if len(perms) == 0:
            raise RuntimeError("No valid permutations available")
        print("permutations established")

        for line in datafile:
            gene_id, _, remainder = line.partition(",")
            values = np.fromstring(remainder, sep=",")

            gene_stats = analyze_gene_paired(
                values,
                perms,
                mask_ab,
            )

            gene_stats['gene_id'] = gene_id
            entries.append(
                gene_stats
            )

            n_line += 1

            if n_line % n_block == 0:
                time_prev = report_progress(
                    n_block,
                    time_prev
                )

    with open(out_file_path, "w") as outfile:

        outfile.write(
            ",".join(generate_header(a_factor, b_factor)) + "\n"
        )

        if args.correction != 'none':
            print(f"Multiple testing correction: {args.correction}")
            entries = multiple_test_correction(
                entries,
                correction=args.correction
            )

        outfile.write(
            '\n'.join(
                format_entry(entry) for entry in entries)
        )

    return ''


# ###################
# MAIN
# ####################


def main(arguments=sys.argv[1:]):
    TIME0 = time_now()
    args = generate_argparse(arguments)
    metadata = parse_factors(args.metadata)

    transform_filter_data(
        data_filepath=args.input,
        mode=args.transform,
        minimum_average=args.min_avg
        )

    contrasts_a = parse_comma_dash(args.contrasts_a)
    contrasts_b = parse_comma_dash(args.contrasts_b)

    print(f"Initialization seed: {args.seed}")
    print(f"Contrasts A: {contrasts_a}")
    print(f"Contrasts B: {contrasts_b}")

    for a_factor in contrasts_a:
        for b_factor in contrasts_b:
            if args.pair:
                run_contrast_paired(
                    a_factor,
                    b_factor,
                    metadata,
                    args
                )

            else:
                run_contrast_unpaired(
                    a_factor,
                    b_factor,
                    metadata,
                    args
                )

    print(f"Completed in: {(time_now()-TIME0):.1f} seconds")
    return ''


if __name__ == "__main__":
    main()
