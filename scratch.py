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

def process_excel_file(file_path: str, variant_column_name: str, uniprotID: str):
    try:
        df = pd.read_excel(file_path, engine="openpyxl", sheet_name="all", skiprows=1)
    except ValueError:
        df = pd.read_excel(file_path, engine="openpyxl", skiprows=1)
    df["Amino acid change"] = df["Amino acid change"].str.replace("p.", "")
    df = df.apply(lambda x: extract_variant(x, variant_column_name), axis=1)
    df = df[["Position", "Original", "Mutated", "Clinical significance (ClinVar)"]]
    df["UniprotID"] = uniprotID
    return df

if __name__ == "__main__":
    select_uniprots = ["Q5S007", "Q96QK1", "Q9BXM7", "O60260", "P37840", "P04062", "Q99497"]
    lrrk2 = process_excel_file("data/PD Variant Browser_LRRK2 variants in humans_ES.xlsx", "Amino acid change", "Q5S007")
    vps35 = process_excel_file("data/PD Variant Browser_VPS35_ES.xlsx", "Amino acid change", "Q96QK1")
    pink1 = process_excel_file("data/PD Variant Browser PINK1.xlsx", "Amino acid change", "Q9BXM7")
    prkn = process_excel_file("data/PD Variant Browser PRKN.xlsx", "Amino acid change", "O60260")
    snca = process_excel_file("data/PD Variant Browser SNCA .xlsx", "Amino acid change", "P37840")
    gba = process_excel_file("data/PD Variant Browser GBA.xlsx", "Amino acid change", "P04062")
    park7 = process_excel_file("data/PD Variant Browser PARK7 (DJ1).xlsx", "Amino acid change", "Q99497")
    concat_df = pd.concat([lrrk2, vps35, pink1, prkn, snca, gba, park7], ignore_index=True)
    if not os.path.exists("data/alpha.txt"):
        results = []
        with open("data/AlphaMissense_aa_substitutions.tsv", "rt") as alpha:
            reader = csv.reader(alpha, delimiter="\t")
            for line in reader:
                if line[0] in select_uniprots:
                    results.append(line)

        df = pd.DataFrame(results, columns=["UniprotID", "Variant", "Score", "Pathogenicity"])
        df["Score"] = df["Score"].astype(float)

        df = df.apply(lambda x: extract_variant(x, "Variant"), axis=1)
        df.to_csv("data/alpha.txt", sep="\t", index=False)
    else:
        df = pd.read_csv("data/alpha.txt", sep="\t")
    color_discrete_map = {
        "pathogenic": "#ff5671",
        "ambiguous": "#ffc955",
        "benign": "#4e804e"
    }
    results =[]
    parser = UniprotParser(columns="accession,id,ft_domain,sequence")
    for i in parser.parse(select_uniprots):
        results.append(pd.read_csv(io.StringIO(i), sep="\t"))
    if len(results) == 1:
        results = results[0]
    else:
        results = pd.concat(results)

    results = results.apply(lambda x: extract_functional_domain(x, "Domain [FT]"), axis=1)
    for i, row in results.iterrows():
        alpha_df = df[df["UniprotID"] == row["From"]]
        study = concat_df[concat_df["UniprotID"] == row["From"]]
        study = study.merge(alpha_df, on=["Position", "Original", "Mutated"])
        study["Position"] = study["Position"].astype(int)
        fig = make_subplots(rows=2, cols=1, row_heights=[0.2, 0.8], vertical_spacing=0.02, shared_xaxes=True)
        sequence_length = len(row["Sequence"])

        dropdown_option = {k: {} for k in study["Clinical significance (ClinVar)"].unique()}
        dropdown_option["All (PD Browser Variants)"] = {}
        for v in study["Pathogenicity"].unique():
            for k in dropdown_option:
                dropdown_option[k][v] = []

        for j, study_sub in study.groupby(["Pathogenicity"]):
            customtext = []

            for k, r in study_sub.iterrows():
                customtext.append(f"Score: {r['Score']: .2f}<br>Position: {r['Position']}<br>Variant: {r['Variant']}<br>AlphaMissense: {r['Pathogenicity']}<br>PD Variant Browser: {r['Clinical significance (ClinVar)']}")
                for clinical_path in dropdown_option:
                    if r["Clinical significance (ClinVar)"] == clinical_path:
                        dropdown_option[clinical_path][j[0]].append(1)
                    else:
                        dropdown_option[clinical_path][j[0]].append(0.1)
                dropdown_option["All (PD Browser Variants)"][j[0]].append(1)

            fig.add_trace(
                go.Scatter(
                    name=j[0],
                    x=study_sub["Position"],
                    y=study_sub["Score"],
                    mode="markers",
                    text=customtext,
                    customdata=study_sub['Clinical significance (ClinVar)'],
                    marker=dict(
                        color=color_discrete_map[j[0]], size=10, opacity=[1]*study_sub.shape[0]
                    ),
                    hovertemplate="%{text}",
                ), row=2, col=1,
            )

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

        buttons = []
        for k in ["All (PD Browser Variants)"]+list(study["Clinical significance (ClinVar)"].unique()):
            args = [{"marker.opacity": []}, [0, 1, 2]]
            for j, _ in study.groupby(["Pathogenicity"]):
                args[0]["marker.opacity"].append(dropdown_option[k][j[0]])

            buttons.append(dict(
                label=k,
                method="restyle",
                args=args
            ))
        print(buttons)
        fig.update_layout(
            template="ggplot2",
            legend_title_text="Alphamissense Pathogenicity",
            updatemenus=[{
                "active": 0,
                "buttons": buttons
            }]


        )

        fig.update_xaxes(showgrid=False, zeroline=False, range=[1, sequence_length])
        fig.update_yaxes(showgrid=False, zeroline=False, range=[0, 1])
        fig.update_yaxes(showticklabels=False, row=1, col=1)

        fig.write_html(f"data/{row['Entry Name']}.html")



    # lrrk2 = lrrk2.merge(df[df["UniprotID"]==select_uniprots[0]], on=["Position", "Original", "Mutated"])
    # # plot scatter plot with plotly express where x is position and y is score value and color is pathogenicity, color scheme is pastel with pastel red being mapped to pathogenic, pastel orange to ambiguous and pastel green to benign
    #
    #
    # fig = make_subplots(rows=1, cols=1)
    # for i, lrrk2_sub in lrrk2.groupby(["Pathogenicity"]):
    #     fig.add_trace(
    #         go.Scatter(
    #             name=i[0],
    #             x=lrrk2_sub["Position"],
    #             y=lrrk2_sub["Score"],
    #             mode="markers",
    #             text=lrrk2_sub["Variant"],
    #             customdata=lrrk2_sub["Pathogenicity"],
    #             marker=dict(color=[color_discrete_map[i] for i in lrrk2_sub["Pathogenicity"]]),
    #             hovertemplate=
    #             "Score: %{x:.2f}<br>Position: %{y}<br>Variant: %{text}<br>Pathogenicity: %{customdata}",
    #
    #         ), row=1, col=1,
    #     )
    # fig.write_html("data/lrrk2.html")
    # lrrk2.to_csv("data/lrrk2.txt", sep="\t", index=False)
    # vps35 = vps35.merge(df[df["UniprotID"]==select_uniprots[1]], on=["Position", "Original", "Mutated"])
    # fig2 = make_subplots(rows=1, cols=1)
    # for i, vps35_sub in vps35.groupby(["Pathogenicity"]):
    #     fig2.add_trace(
    #         go.Scatter(
    #             name=i[0],
    #             x=vps35_sub["Position"],
    #             y=vps35_sub["Score"],
    #             mode="markers",
    #             text=vps35_sub["Variant"],
    #             customdata=vps35_sub["Pathogenicity"],
    #             marker=dict(color=[color_discrete_map[i] for i in vps35_sub["Pathogenicity"]]),
    #             hovertemplate=
    #             "Score: %{x:.2f}<br>Position: %{y}<br>Variant: %{text}<br>Pathogenicity: %{customdata}",
    #
    #         ), row=1, col=1,
    #     )
    # fig2.write_html("data/lrrk2.html")
    # vps35.to_csv("data/vps35.txt", sep="\t", index=False)

