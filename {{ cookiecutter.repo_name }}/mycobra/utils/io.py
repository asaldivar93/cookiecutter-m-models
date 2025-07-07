from __future__ import annotations

import polars as pl
from cobra import Metabolite, Model, Reaction, Gene
from sympy.polys.fglmtools import _update


def read_excel(file_path: str) -> Model:
    """Read an excel model into cobrapy."""
    reactions = pl.read_excel(file_path, sheet_name="reactions")
    metabolites = pl.read_excel(file_path, sheet_name="metabolites")
    genes = pl.read_excel(file_path, sheet_name="genes")
    model = Model("")

def update_from_excel(model: Model, file_path: str):
    # Get current state of the model
    current_rxns = [rxn.id for rxn in model.reactions]
    current_mets = [met.id for met in model.metabolites]
    current_genes = [gene.id for gene in model.genes]

    # Get the updates
    file_path="data/interim/updates_to_model.xlsx"
    reactions = pl.read_excel(file_path, sheet_name="reactions")
    metabolites = pl.read_excel(file_path, sheet_name="metabolites")
    genes = pl.read_excel(file_path, sheet_name="genes")
    deletions = pl.read_excel(file_path, sheet_name="delitions")
    bounds = pl.read_excel(file_path, sheet_name="bounds")

    # Get databases in use
    rxn_databases = reactions[: , "biocyc":].columns
    met_databases = metabolites[: , "biocyc":].columns
    gene_databases = metabolites[: , "ncbiprotein":].columns

    # Update the model by day
    updates_list = []
    dates = reactions["date"].unique()
    for day in dates:
        this_update = {}
        # Get updates by day
        these_rxns = reactions.filter(pl.col("date")==day)
        these_mets = metabolites.filter(pl.col("date")==day)
        these_genes = genes.filter(pl.col("date")==day)
        these_bounds = bounds.filter(pl.col("date")==day)
        these_dels = deletions.filter(pl.col("date")==day)

        # Delete unused reactions
        rxns_to_del = []
        mets_to_del = []
        gens_to_del = []
        for row in these_dels.iter_rows(named=True)
            match row["type"].lower():
                case "reaction":
                    rxns_to_del.extend([row["id"]])
                case "metabolite":
                    mets_to_del.extend([row["id"]])
                case "gene":
                    gens_to_del.extend([row["id"]])
                case _:
                    print(f"Type error in delitions id:{row['id']}")
        model.remove_reactions(rxns_to_del)
        model.remove_metabolites(mets_to_del)

        # Update metabolites
        mets_to_add = []
        for row in these_mets.iter_rows(named=True):
            met_id = row["met_id"]
            name = row["name"]
            if met_id in current_mets:
                met = model.metabolites.get_by_id(met_id)
            else:
                met = Metabolite(met_id)
                mets_to_add.append(met)

            met.name = row["name"]
            met.formula = row["formula"]
            met.charge = row["charge"]
            met.compartment = row["compartment"]

            for db in met_databases:
                met.annotation[db] = row[db]
        model.add_metabolites(mets_to_add)

        # Update reactions
        rxns_to_add = []
        for row in these_rxns.iter_rows(named=True):
            rxn_id = row["rxn_id"]
            name = row["name"]
            if rxn_id in current_rxns:
                rxn = model.reactions.get_by_id(rxn_id)
            else:
                rxn = Reaction(rxn_id, name)
                rxns_to_add.append(rxn)

            rxn.build_reaction_from_string(row["reaction"])
            rxn.gene_reaction_rule = row["gpr"]
            rxn.subsystem = row["subsystem"]
            rxn.lower_bound = row["lower_bound"]
            rxn.upper_bound = row["upper_bound"]
            for db in rxn_databases:
                rxn.annotation[db] = row[db]
        model.add_reactions(rxns_to_add)

        unbalanced_rxn = [rxn.id for rxn in model.reactions
                                 if rxn.check_mass_balance()]

        # Summary of the updates
        this_update["date"] = day
        this_update["metabolites"] = len(model.metabolites)
        this_update["reactions"] = len(model.reactions)
        this_update["genes"] = len(model.genes)
        this_update["unbalanced"] = len(unbalanced_rxn)
        this_update["growth_rate"] = model.slim_optimize()
        updates_list.append(this_update)

    updates_df = pl.from_dicts(updates_list)
    return model, updates_df







import os
os.getcwd()
os.chdir('/home/asaldivargarci1064/cookiecutter-m-models/{{ cookiecutter.repo_name }}')
