from copy import deepcopy

import pandas as pd


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