import csv
import dataclasses
import json
import os
import sys
import time

import ray

import pyscreener as ps


def check():
    ps.check_env(sys.argv[1], json.loads(sys.argv[2]))
    exit(0)


def main():
    args = ps.args.gen_args()

    if args.smoke_test:
        ps.check_env(args.screen_type, args.metadata_template)
        exit(0)

    print(
        """\
***************************************************************
*      ____  __  ____________________  ___  ____  ___  _____  *
*     / __ \/ / / / ___/ ___/ ___/ _ \/ _ \/ __ \/ _ \/ ___/  *
*    / /_/ / /_/ (__  ) /__/ /  /  __/  __/ / / /  __/ /      *
*   / .___/\__, /____/\___/_/   \___/\___/_/ /_/\___/_/       *
*  /_/    /____/                                              *
***************************************************************"""
    )
    print("Welcome to Pyscreener!\n")

    params = vars(args)
    print("Pyscreener will be run with the following arguments:")
    for param, value in sorted(params.items()):
        print(f"  {param}: {value}")
    print(flush=True)

    try:
        if "redis_password" in os.environ:
            ray.init(
                address=os.environ["ip_head"],
                _node_ip_address=os.environ["ip_head"].split(":")[0],
                _redis_password=os.environ["redis_password"],
            )
        else:
            ray.init(address="auto")
    except ConnectionError:
        ray.init()
    except PermissionError:
        print("Failed to create a temporary directory for ray")
        raise

    print("Ray cluster online with resources:")
    print(ray.cluster_resources())
    print(flush=True)

    print("Preparing and screening inputs ...", flush=True)
    metadata_template = ps.build_metadata(args.screen_type, args.metadata_template)
    virtual_screen = ps.virtual_screen(
        args.screen_type,
        args.receptors,
        args.center,
        args.size,
        metadata_template,
        args.pdbids,
        args.docked_ligand_file,
        args.buffer,
        args.ncpu,
        args.base_name,
        args.output_dir,
        args.score_mode,
        args.repeat_score_mode,
        args.ensemble_score_mode,
        args.repeats,
        args.k,
        args.verbose,
    )
    supply = ps.LigandSupply(
        args.input_files,
        args.input_filetypes,
        args.smis,
        args.use_3d,
        args.optimize,
        args.title_line,
        args.smiles_col,
        args.name_col,
        args.id_property,
    )
    start = time.time()

    S = virtual_screen(supply.ligands)

    total_time = time.time() - start
    print("Done!")

    avg_time = total_time / len(virtual_screen)
    m, s = divmod(total_time, 60)
    h, m = divmod(int(m), 60)
    print(
        f"Total time to dock {len(virtual_screen)} ligands: {h}h {m}m {s:0.2f}s "
        f"({avg_time:0.2f}s/ligand)"
    )

    if args.hist_mode is not None:
        ps.postprocessing.histogram(
            args.hist_mode, S, virtual_screen.path, "score_distribution.png"
        )

    results = virtual_screen.all_results()
    if not args.no_sort:
        results = sorted(results, key=lambda r: r.score if r.score is not None else float("inf"))
    smis_scores = [(r.smiles, r.score) for r in results]

    scores_filename = virtual_screen.path / "scores.csv"
    with open(scores_filename, "w") as fid:
        writer = csv.writer(fid)
        writer.writerow(["smiles", "score"])
        writer.writerows(smis_scores)

    print(f'Scoring data has been saved to: "{scores_filename}"')

    if args.collect_all:
        print("Collecting all input and output files ...", end=" ", flush=True)
        virtual_screen.collect_files()

        extended_filename = virtual_screen.path / "extended.csv"
        with open(extended_filename, "w") as fid:
            writer = csv.writer(fid)
            writer.writerow(field.name for field in dataclasses.fields(results[0]))
            writer.writerows(dataclasses.astuple(r) for r in results)

        print("Done!")
        print(f'Extended data has been saved to: "{extended_filename}"')

    print("Thanks for using Pyscreener!")


if __name__ == "__main__":
    main()
