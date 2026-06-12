#!/usr/bin/env python3
"""Export all retrosynthesis sequences (routes) from IBM RXN for Chemistry to CSV.

The IBM RXN web UI lets you browse routes one at a time via the "Sequence N"
dropdown and export the *current* route to PDF, but there is no single button to
download a table of *all* sequences. This script uses the official RXN4Chemistry
API to fetch every predicted route at once and flatten them into one CSV.

Usage:
    export RXN_API_KEY="your-api-key"        # from your RXN profile page
    python3 rxn_export_sequences.py "COc1ccc(CCN2CCC(Nc3nc4ccccc4n3Cc3ccc(F)cc3)CC2)cc1"
    python3 rxn_export_sequences.py <SMILES> --out routes.csv --project-id <id>

Install dependency first:
    pip install rxn4chemistry

CSV columns:
    sequence       index of the route (0-based, matches the "Sequence N" UI)
    sequence_id    RXN sequenceId for the route
    confidence     overall confidence of the route (0..1)
    n_steps        number of reaction steps in the route
    step           reaction step number within the route (1-based)
    reaction_smiles  reactants>>product for that step
    starting_materials  '.'-joined leaf SMILES (commercial building blocks)
"""

import argparse
import csv
import os
import sys
import time

try:
    from rxn4chemistry import RXN4ChemistryWrapper
except ImportError:
    sys.exit(
        "Missing dependency. Install it with:\n    pip install rxn4chemistry"
    )


def load_api_key():
    """Resolve the API key without ever hardcoding it in this committed file.

    Priority: RXN_API_KEY env var, then a local untracked '.rxn_key' file next
    to this script (kept out of git via .gitignore).
    """
    key = os.environ.get("RXN_API_KEY")
    if key:
        return key.strip()
    local = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".rxn_key")
    if os.path.exists(local):
        with open(local) as fh:
            return fh.read().strip()
    return None


def iter_reactions(node):
    """Yield (reaction_smiles, ) for each reaction node in a retro tree.

    A non-leaf node is a product whose children are its reactants, so the
    reaction is 'reactant1.reactant2...>>product'.
    """
    children = node.get("children") or []
    if children:
        reactants = ".".join(c["smiles"] for c in children)
        yield f"{reactants}>>{node['smiles']}"
        for child in children:
            yield from iter_reactions(child)


def iter_leaves(node):
    """Yield SMILES of leaf nodes (starting materials / building blocks)."""
    children = node.get("children") or []
    if not children:
        yield node["smiles"]
    for child in children:
        yield from iter_leaves(child)


def flatten_path(index, path):
    """Turn one retrosynthetic path into a list of per-step CSV rows."""
    tree = path.get("retrosynthesis") or path  # API nests the tree under either
    reactions = list(iter_reactions(tree))
    starting_materials = ".".join(sorted(set(iter_leaves(tree))))
    rows = []
    for step, rxn in enumerate(reactions, start=1):
        rows.append(
            {
                "sequence": index,
                "sequence_id": path.get("sequenceId", ""),
                "confidence": path.get("confidence", ""),
                "n_steps": len(reactions),
                "step": step,
                "reaction_smiles": rxn,
                "starting_materials": starting_materials,
            }
        )
    if not rows:  # single-node tree (no disconnection found)
        rows.append(
            {
                "sequence": index,
                "sequence_id": path.get("sequenceId", ""),
                "confidence": path.get("confidence", ""),
                "n_steps": 0,
                "step": 0,
                "reaction_smiles": "",
                "starting_materials": starting_materials,
            }
        )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("smiles", help="SMILES of the target molecule")
    parser.add_argument("--out", default="rxn_sequences.csv", help="output CSV path")
    parser.add_argument(
        "--api-key",
        default=load_api_key(),
        help="RXN API key (or set RXN_API_KEY env var, or a local .rxn_key file)",
    )
    parser.add_argument(
        "--project-id",
        default=os.environ.get("RXN_PROJECT_ID"),
        help="optional RXN project id; a temporary one is created if omitted",
    )
    parser.add_argument(
        "--timeout", type=int, default=300, help="seconds to wait for results"
    )
    args = parser.parse_args()

    if not args.api_key:
        sys.exit("No API key. Pass --api-key or set RXN_API_KEY.")

    rxn = RXN4ChemistryWrapper(api_key=args.api_key)

    if args.project_id:
        rxn.set_project(args.project_id)
    else:
        resp = rxn.create_project(f"retrosynthesis_export_{int(time.time())}")
        if not getattr(rxn, "project_id", None):
            # Surface the real cause (auth failure, blocked host, etc.) instead
            # of the cryptic "Project identifier has to be set first." later on.
            detail = resp.get("response") if isinstance(resp, dict) else resp
            sys.exit(f"Could not create RXN project: {detail}")

    print(f"Submitting retrosynthesis for: {args.smiles}", file=sys.stderr)
    response = rxn.predict_automatic_retrosynthesis(product=args.smiles)
    prediction_id = response.get("prediction_id") if isinstance(response, dict) else None
    if not prediction_id:
        # Surface the API's actual reply (auth error, quota, bad SMILES, etc.)
        # instead of a bare KeyError.
        sys.exit(f"No prediction_id returned. Full API response:\n{response}")

    print("Waiting for results", end="", file=sys.stderr, flush=True)
    deadline = time.time() + args.timeout
    results = None
    while time.time() < deadline:
        results = rxn.get_predict_automatic_retrosynthesis_results(prediction_id)
        status = results.get("status")
        if status == "SUCCESS":
            break
        print(".", end="", file=sys.stderr, flush=True)
        time.sleep(10)
    print("", file=sys.stderr)

    if not results or results.get("status") != "SUCCESS":
        sys.exit(f"Prediction did not finish in {args.timeout}s (status: "
                 f"{results.get('status') if results else 'unknown'}).")

    paths = results.get("retrosynthetic_paths") or []
    if not paths:
        sys.exit("No retrosynthetic paths returned.")

    fieldnames = [
        "sequence",
        "sequence_id",
        "confidence",
        "n_steps",
        "step",
        "reaction_smiles",
        "starting_materials",
    ]
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for index, path in enumerate(paths):
            for row in flatten_path(index, path):
                writer.writerow(row)

    print(f"Wrote {len(paths)} sequence(s) to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
