import csv
import matplotlib.pyplot as plt
import statsmodels.api as sm
import seaborn as sns
from statsmodels.stats.contingency_tables import Table
import pandas as pd
from Patient import Patient


if __name__ == "__main__":
    # 1. Load patients from MetaData.csv
    Patient.instantiate_from_csv("MetaData.csv")
    print("Loaded patients:")
    for patient in Patient.all_patients:
        print(patient)


    # 2. Print all patients to verify
    cerad_order = {"Absent": 0, "Sparse": 1, "Moderate": 2, "Frequent": 3}
    Patient.all_patients.sort(
        key=lambda p: cerad_order.get(p.CERAD_score.strip(), 99)  # 99 = "unknown" catch
    )
    # Print sorted list
    print("\nPatients sorted by CERAD score:")
    for patient in Patient.all_patients:
        print(patient, getattr(patient, "CERAD_score", "N/A"))


    #3. Sorted dictionary
    def sort_by_cerad(cls):
        # Ensure dictionary is clean
        cls.cerad_groups.clear()
        for patient in cls.all_patients:
            cerad = patient.CERAD_score.strip()
            if cerad not in cls.cerad_groups:
                cls.cerad_groups[cerad] = []
            cls.cerad_groups[cerad].append(patient)

    def subsort_by_age(cls):
        for cerad, patients in cls.cerad_groups.items():
            # Sort in place by numeric age (invalid -> big number at end)
            patients.sort(
                key=lambda p: int(p.age_at_death) if str(p.age_at_death).isdigit() else 999
            )
            cls.cerad_groups[cerad] = patients    


    #.4 Bar graph (CERAD scores vs number of patients)
    # Filter only patients with dementia
    dementia_patients = [
        p for p in Patient.all_patients 
        if p.cognitive_status and "dementia" in p.cognitive_status.lower() 
        and p.CERAD_score
    ]

    # Count how many dementia patients have each CERAD score
    cerad_order = ["Absent", "Sparse", "Moderate", "Frequent"]
    counts = {cerad: 0 for cerad in cerad_order}

    for p in dementia_patients:
        cerad = p.CERAD_score.strip().title()
        if cerad in counts:
            counts[cerad] += 1

    # Convert counts to percentages
    total_dementia = len(dementia_patients)
    percentages = {cerad: (count / total_dementia) * 100 for cerad, count in counts.items()}

    # Plot
    plt.figure(figsize=(6, 4))
    plt.bar(percentages.keys(), percentages.values(), color="salmon", edgecolor="black")
    plt.ylabel("Percentage of Dementia Patients (%)")
    plt.xlabel("CERAD Score")
    plt.title("CERAD Score Distribution Among Dementia Patients")
    plt.ylim(0, 40)
    plt.show()

    #6. percentage of demntia patients vs apoe genotype
    # Filter only patients with dementia
    dementia_patients = [
        p for p in Patient.all_patients
        if p.cognitive_status and "dementia" in p.cognitive_status.lower()
        and p.apoe
    ]

    # Count how many dementia patients have each APOE genotype
    apoe_order = ["2/3", "3/3", "3/4", "4/4"]
    counts = {apoe: 0 for apoe in apoe_order}

    for p in dementia_patients:
        genotype = p.apoe.strip()
        if genotype in counts:
            counts[genotype] += 1

    # Convert counts to percentages
    total_dementia = len(dementia_patients)
    percentages = {genotype: (count / total_dementia) * 100 for genotype, count in counts.items()}

    # Plot
    plt.figure(figsize=(6, 4))
    plt.bar(percentages.keys(), percentages.values(), color="skyblue", edgecolor="black")
    plt.ylabel("Percentage of Dementia Patients (%)")
    plt.xlabel("APOE Genotype")
    plt.title("APOE Genotype Distribution Among Dementia Patients")
    plt.ylim(0, 60)
    plt.show()

    #5. percentage of demntia patients vs apoe genotype
    # Filter only patients with dementia
    dementia_patients = [
        p for p in Patient.all_patients
        if p.cognitive_status and "dementia" in p.cognitive_status.lower()
        and p.apoe
    ]

    # Count how many dementia patients have each APOE genotype
    apoe_order = ["2/3", "3/3", "3/4", "4/4"]
    counts = {apoe: 0 for apoe in apoe_order}

    for p in dementia_patients:
        genotype = p.apoe.strip()
        if genotype in counts:
            counts[genotype] += 1

    # Convert counts to percentages
    total_dementia = len(dementia_patients)
    percentages = {genotype: (count / total_dementia) * 100 for genotype, count in counts.items()}

    # Plot
    plt.figure(figsize=(6, 4))
    plt.bar(percentages.keys(), percentages.values(), color="skyblue", edgecolor="black")
    plt.ylabel("Percentage of Dementia Patients (%)")
    plt.xlabel("APOE Genotype")
    plt.title("APOE Genotype Distribution Among Dementia Patients")
    plt.ylim(0, 60)
    plt.show()

    #6. APOE Genotype vs CERAD scores in the set of dementia patients
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt

    # Filter only patients with dementia
    dementia_patients = [
        p for p in Patient.all_patients
        if p.cognitive_status and "dementia" in p.cognitive_status.lower()
        and p.apoe and p.CERAD_score
    ]

    # Create DataFrame
    data = []
    for p in dementia_patients:
        data.append({
            "APOE": p.apoe.strip(),
            "CERAD": p.CERAD_score.strip().title()
        })

    df = pd.DataFrame(data)

    # Crosstab: rows = APOE, columns = CERAD, values = counts
    crosstab = pd.crosstab(df["APOE"], df["CERAD"])

    # Ensure all CERAD and APOE categories are included
    cerad_order = ["Absent", "Sparse", "Moderate", "Frequent"]
    apoe_order = ["2/3", "3/3", "3/4", "4/4"]

    for cerad in cerad_order:
        if cerad not in crosstab.columns:
            crosstab[cerad] = 0

    crosstab = crosstab[cerad_order]  # reorder columns
    crosstab = crosstab.reindex(apoe_order)  # reorder rows

    # Convert counts to percentages per APOE genotype
    percent_crosstab = crosstab.div(crosstab.sum(axis=1), axis=0) * 100

    # Plot heatmap
    plt.figure(figsize=(8, 5))
    sns.heatmap(
        percent_crosstab,
        annot=True, fmt=".1f", cmap="YlOrRd",
        cbar_kws={'label': 'Percentage of Dementia Patients (%)'}
    )
    plt.xlabel("CERAD Score")
    plt.ylabel("APOE Genotype")
    plt.title("CERAD Score Distribution Across APOE Genotypes (Dementia Patients)")
    plt.show()
