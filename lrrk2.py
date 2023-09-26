import io
import os.path
import re
import csv
import pandas as pd
from uniprotparser.betaparser import UniprotParser, UniprotSequence
from copy import deepcopy
from plotly.subplots import make_subplots
import plotly.graph_objects as go
# a regex pattern that capture the variant data from the variant column in the form of A123B where the number should be extracted as position, A would be extracted as original residue while B would be mutated residue

pattern = re.compile(r"([A-Z]+)(\d+)([A-Z]+)")
custom_domain = {
    "Q5S007": [
        {"start": 1, "end": 705, "domain": "ARM"},
        {"start": 706, "end": 800, "domain": "ANK"},
        {"start": 801, "end": 1335, "domain": "LRR"},
        {"start": 1336, "end": 1511, "domain": "ROC"},
        {"start": 1512, "end": 1879, "domain": "COR"},
        {"start": 1880, "end": 2142, "domain": "KIN"},
        {"start": 2143, "end": 2498, "domain": "WD40"},
    ]
}
color_discrete_map = {
    "pathogenic": "#ff5671",
    "ambiguous": "#ffc955",
    "benign": "#4e804e",
    "MDSGene - probable pathogenic": "#4a54db",
    "MDSGene - possibly pathogenic": "#29eeed",
    "MDSGene - clearly pathogenic mutations": "#bc64ff",
}
def extract_variant(row: pd.Series, variant_column_name: str):
    search_result = pattern.search(row[variant_column_name])
    if search_result:
        row["Position"] = int(search_result.group(2))
        row["Original"] = search_result.group(1)
        row["Mutated"] = search_result.group(3)

    return row

def extract_functional_domain(row: pd.Series, domain_column_name: str):
    domains = []
    if pd.notnull(row[domain_column_name]):
        splitted = row[domain_column_name].split(";")
        current_domain = ""
        current_domain_start = ""
        current_domain_end = ""
        previous_domain_end = None
        for i in splitted:
            i = i.strip()
            if i.startswith("DOMAIN"):
                if current_domain:
                    domains.append(deepcopy({
                        "domain": current_domain,
                        "start": current_domain_start,
                        "end": current_domain_end
                    }))
                current_domain = ""
                stripped_domain = i.replace("DOMAIN", "").strip()
                splitted_start_end = stripped_domain.split("..")
                try:
                    current_domain_start = int(splitted_start_end[0])
                    print(i)
                    if previous_domain_end is None:
                        if current_domain_start > 1:
                            domains.append(deepcopy({
                                "domain": "Other",
                                "start": 1,
                                "end": current_domain_start - 1
                            }))
                        else:
                            pass
                    elif (current_domain_start - previous_domain_end) > 1:
                        domains.append(deepcopy({
                            "domain": "Other",
                            "start": previous_domain_end + 1,
                            "end": current_domain_start - 1
                        }))
                except ValueError:
                    current_domain_start = None
                try:
                    current_domain_end = int(splitted_start_end[1])
                    previous_domain_end = int(splitted_start_end[1])
                except ValueError:
                    current_domain_end = None
            elif i.startswith("/note="):
                current_domain = i.replace("/note=", "").replace("\"", "").strip()
        if current_domain != "":
            domains.append(deepcopy({
                "domain": current_domain,
                "start": current_domain_start,
                "end": current_domain_end
            }))
        if previous_domain_end is None or previous_domain_end < len(row["Sequence"]):
            domains.append(deepcopy({
                "domain": "Other",
                "start": current_domain_end + 1,
                "end": len(row["Sequence"])
            }))

    row["domains"] = domains

    return row

if __name__ == "__main__":
    uniprotIDs = ["Q5S007"]
    column_with_mutation = "SNV Mutation"

    mds = pd.read_excel("data/NEW-For TOAN LRRK2 25Sep2023.xlsx", sheet_name="MDSGene classification")
    mds = mds[pd.notnull(mds["MDSGene"])]

    keep = [
        "SNV Mutation",
        "ClinVar Clinical Significance",
        "ClinVar Variation ID",
        "Allele Count",
        "Allele Number",
        "Allele Frequency",
        "Homozygote Count",
        "Hemizygote Count",
    ]


    gnomAd = pd.read_excel("data/NEW-For TOAN LRRK2 25Sep2023.xlsx", sheet_name="Missense gnomAd")
    gnomAd = gnomAd[keep]
    gnomAd["SNV Mutation"] = gnomAd["SNV Mutation"].str.upper()
    gnomAd = gnomAd.merge(mds, left_on="SNV Mutation", right_on="MDSGene", how="left")
    gnomAd["ClinVar Clinical Significance"] = gnomAd["ClinVar Clinical Significance"].fillna("")
    gnomAd["classification"] = gnomAd["classification"].fillna("")
    gnomAd["classification"] = gnomAd["classification"].str.strip()
    gnomAd = gnomAd.apply(lambda x: extract_variant(x, column_with_mutation), axis=1)
    hover_text = []
    for i, row in gnomAd.iterrows():
        hover_text.append("<br>".join([f"{c}: {row[c]}" for c in gnomAd.columns if c not in ["Original", "Mutated"] and pd.notnull(row[c])]))
    gnomAd["hover_text"] = hover_text
    gnomAd["UniprotID"] = "Q5S007"
    gnomAd["Position"] = gnomAd["Position"].astype(int)
    parser = UniprotParser(columns="accession,id,ft_domain,sequence")
    uniprot_results = []
    for i in parser.parse(uniprotIDs):
        uniprot_results.append(pd.read_csv(io.StringIO(i), sep="\t"))
    if len(uniprot_results) == 1:
        uniprot_results = uniprot_results[0]
    else:
        uniprot_results = pd.concat(uniprot_results)
    uniprot_results = uniprot_results.apply(lambda x: extract_functional_domain(x, "Domain [FT]"), axis=1)

    results = []
    if not os.path.exists("data/alpha.txt"):
        with open("data/AlphaMissense_aa_substitutions.tsv", "rt") as alpha:
            reader = csv.reader(alpha, delimiter="\t")
            for line in reader:
                if line[0] in uniprotIDs:
                    results.append(line)

        df = pd.DataFrame(results, columns=["UniprotID", "Variant", "Score", "Pathogenicity"])
        df["Score"] = df["Score"].astype(float)
        df.to_csv("data/alpha.txt", sep="\t", index=False)
    else:
        df = pd.read_csv("data/alpha.txt", sep="\t")
    df = df.apply(lambda x: extract_variant(x, "Variant"), axis=1)
    df["Position"] = df["Position"].astype(int)
    df = df.merge(gnomAd, on=["UniprotID", "Position", "Original", "Mutated"], how="left")

    temp_dict = {}

    for ui, row in uniprot_results.iterrows():
        alpha_df = df[df["UniprotID"] == row["From"]]

        dropdown_option = {k: {} for k in alpha_df["ClinVar Clinical Significance"].unique()}
        dropdown_option["All (ClinVar)"] = {}
        for v in alpha_df["Pathogenicity"].unique():
            for k in dropdown_option:
                dropdown_option[k][v] = []
        dropdown_mds_option = {k: {} for k in alpha_df["classification"].unique()}
        dropdown_mds_option["All (MDS class)"] = {}
        for v in alpha_df["Pathogenicity"].unique():
            for k in dropdown_mds_option:
                dropdown_mds_option[k][v] = []
        sequence_length = len(row["Sequence"])
        for i, g in alpha_df.groupby("Pathogenicity"):
            if i not in temp_dict:
                temp_dict[i] = {
                    "x": [],
                    "y": [],
                    "text": [],
                    "opacity": [],
                    "color": [],
                    "skip": False
                }
                temp_dict[i + " only in Alphamissense"] = {
                    "x": [],
                    "y": [],
                    "text": [],
                    "opacity": [],
                    "color": [],
                    "skip": True
                }

            for i2, row2 in g.iterrows():

                dropdown_option["All (ClinVar)"][i].append(1)
                dropdown_mds_option["All (MDS class)"][i].append(1)
                name = i
                if pd.isna(row2["hover_text"]):
                    name = i + " only in Alphamissense"
                temp_dict[name]["x"].append(row2["Position"])
                temp_dict[name]["y"].append(row2["Score"])
                if pd.notnull(row2["hover_text"]):
                    for clinical_path in dropdown_option:
                        if row2["ClinVar Clinical Significance"] == clinical_path:
                            dropdown_option[clinical_path][i].append(1)
                        else:
                            dropdown_option[clinical_path][i].append(0)
                    if row2["classification"] != "":
                        if f"MDSGene - {row2['classification']}" not in temp_dict:
                            temp_dict[f"MDSGene - {row2['classification']}"] = {
                                "x": [],
                                "y": [],
                                "text": [],
                                "opacity": [],
                                "color": [],
                                "skip": False
                            }
                        temp_dict[f"MDSGene - {row2['classification']}"]["x"].append(row2["Position"])
                        temp_dict[f"MDSGene - {row2['classification']}"]["y"].append(row2["Score"])
                        temp_dict[f"MDSGene - {row2['classification']}"]["text"].append(row2["hover_text"] + "<br>" + f"Score: {row2['Score']: .2f}<br>" + f"Alphamissense: {row2['Pathogenicity']}")
                        temp_dict[f"MDSGene - {row2['classification']}"]["opacity"].append(1)
                        temp_dict[f"MDSGene - {row2['classification']}"]["color"].append(color_discrete_map[f"MDSGene - {row2['classification']}"])

                    for mds_class in dropdown_mds_option:
                        if row2["classification"] == mds_class:
                            dropdown_mds_option[mds_class][i].append(1)
                        else:
                            dropdown_mds_option[mds_class][i].append(0)
                    temp_dict[name]["text"].append(row2["hover_text"] + "<br>" + f"Score: {row2['Score']: .2f}<br>" + f"Alphamissense: {row2['Pathogenicity']}")
                    temp_dict[name]["opacity"].append(1)
                else:
                    temp_dict[name]["text"].append(f"Position: {row2['Position']}<br>Variant: {row2['Variant']}<br>Score: {row2['Score']: .2f}<br>Alphamissense: {row2['Pathogenicity']}")
                    temp_dict[name]["opacity"].append(0.05)
                temp_dict[name]["color"].append(color_discrete_map[row2["Pathogenicity"]])
        fig = make_subplots(rows=2, cols=1, row_heights=[0.2, 0.8], vertical_spacing=0.02, shared_xaxes=True)
        keys = list(temp_dict.keys())
        keys.sort(reverse=True)

        for k in keys:
            if not temp_dict[k]["skip"]:
                fig.add_trace(go.Scatter(
                    x=temp_dict[k]["x"],
                    y=temp_dict[k]["y"],
                    text=temp_dict[k]["text"],
                    mode="markers",
                    marker=dict(
                        color=temp_dict[k]["color"],
                        opacity=temp_dict[k]["opacity"],
                        line=dict(width=2,
                                  color='DarkSlateGrey'),
                        size=10,
                    ),
                    hoverinfo="text",
                    name=k

                ), row=2, col=1)
            else:
                fig.add_trace(go.Scatter(
                    x=temp_dict[k]["x"],
                    y=temp_dict[k]["y"],
                    text=temp_dict[k]["text"],
                    mode="markers",
                    marker=dict(
                        color=temp_dict[k]["color"],
                        opacity=temp_dict[k]["opacity"],
                        size=10,
                    ),
                    hoverinfo="skip",
                    name=k
                ), row=2, col=1)
        buttons = []
        for k in ["All (ClinVar)"] + list(alpha_df["ClinVar Clinical Significance"].unique()):
            args = [{"marker.opacity": []}, []]
            for j in range(len(fig.data)):
                if fig.data[j]["name"] in dropdown_option[k]:
                    args[0]["marker.opacity"].append(dropdown_option[k][fig.data[j]["name"]])
                    args[1].append(j)

            buttons.append(dict(
                label=k,
                method="restyle",
                args=args
            ))
        buttons2 = []

        for k in ["All (MDS class)"] + list(alpha_df["classification"].unique()):
            args = [{"marker.opacity": []}, []]
            for j in range(len(fig.data)):
                if fig.data[j]["name"] in dropdown_mds_option[k]:
                    args[0]["marker.opacity"].append(dropdown_mds_option[k][fig.data[j]["name"]])
                    args[1].append(j)

            buttons2.append(dict(
                label=k,
                method="restyle",
                args=args
            ))
        for i in range(0, sequence_length+1, 50):
            fig.add_shape(
                type="line",
                xref="x",
                yref="y",
                x0=i,
                y0=0,
                x1=i,
                y1=1,
                line=dict(
                    color="black",
                    width=1,
                ),
                layer="below",
                row=2, col=1
            )
            fig.add_shape(
                type="line",
                xref="x",
                yref="y",
                x0=i,
                y0=0,
                x1=i,
                y1=1,
                line=dict(
                    color="black",
                    width=1,
                ),
                layer="below",
                row=1, col=1
            )
        if row["From"] in custom_domain:
            row["domains"] = custom_domain[row["From"]]
        if len(row["domains"]) > 0:
            for d in row["domains"]:
                fig.add_shape(
                    type="rect",
                    xref="x",
                    yref="y",
                    x0=d["start"],
                    y0=0,
                    x1=d["end"],
                    y1=1,
                    fillcolor="#fdfffb",
                    line=dict(
                        color="black",
                    ),
                    opacity=0.5,
                    layer="below",
                    line_width=3, row=2, col=1
                )
                fig.add_shape(
                    type="rect",
                    xref="x",
                    yref="y",
                    x0=d["start"],
                    y0=0,
                    x1=d["end"],
                    y1=1,
                    fillcolor="#fdfffb",
                    line=dict(
                        color="black",
                    ),
                    opacity=0.5,
                    layer="below",
                    line_width=3, row=1, col=1
                )
                fig.add_annotation(
                    text=f"{d['domain']}<br>{d['start']}-{d['end']}",
                    x=(d["start"]+d["end"])/2,
                    y=0.5,
                    showarrow=False,
                    yshift=10,

                    font=dict(
                        size=16,
                        color="black"
                    ),
                    align="center", row=1, col=1
                )
        else:
            fig.add_annotation(
                text="No domain data",
                x=sequence_length/2,
                y=0.5,
                showarrow=False,
                yshift=10,

                font=dict(
                    size=16,
                    color="black"
                ),
                align="center", row=1, col=1
            )

        fig.update_layout(
            template="ggplot2",
            legend_title_text="Alphamissense Pathogenicity",
            updatemenus=[{
                "active": 0,
                "buttons": buttons,
                "direction": "down",
                "y": 1,
                },
                # {
                #     "active": 0,
                #     "buttons": buttons2,
                #     "direction": "down",
                #     "y": 0.6,
                # }
            ]
        )

        fig.update_xaxes(showgrid=False, zeroline=False, range=[1, sequence_length])
        fig.update_yaxes(showgrid=False, zeroline=False, range=[0, 1])
        fig.update_yaxes(showticklabels=False, row=1, col=1)
        fig.write_html("data/{}.html".format(row["Entry Name"]))
        fig.show()

