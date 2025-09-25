import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from Patient import Patient 

if __name__ == "__main__":
    # 1. Load patients from MetaData.csv
    Patient.all_patients.clear()  # avoid duplicates if re-run
    Patient.instantiate_from_csv("MetaData.csv")

    # 2. Print all patients to verify
    print("Loaded patients:")
    for patient in Patient.all_patients:
        print(patient)

    # 3. Count patients by APOE genotype
    genotype_counts = {}
    for patient in Patient.all_patients:
        genotype_counts[patient.apoe] = genotype_counts.get(patient.apoe, 0) + 1

    print("\nAPOE Genotype counts:")
    for genotype, count in genotype_counts.items():
        print(f"{genotype}: {count}")

    # 4. CERAD histogram (numeric only, if possible)
    cerad_numeric_map = {"Absent": 0, "Sparse": 1, "Moderate": 2, "Frequent": 3}

    cerad_scores = []
    for p in Patient.all_patients:
        score = p.CERAD_score.strip()
        if score in cerad_numeric_map:
            cerad_scores.append(cerad_numeric_map[score])

    if cerad_scores:
        plt.hist(cerad_scores, bins=range(min(cerad_scores), max(cerad_scores) + 2), edgecolor="black")
        plt.title("Distribution of CERAD Scores")
        plt.xlabel("CERAD Score (0=Absent, 1=Sparse, 2=Moderate, 3=Frequent)")
        plt.ylabel("Number of Patients")
        plt.show()
    else:
        print("No valid CERAD scores found for plotting.")

    # 5. Bar graph (CERAD scores vs number of patients)
    cerad_categories = ["Absent", "Sparse", "Moderate", "Frequent"]

    counts = [
        sum(1 for p in Patient.all_patients if p.CERAD_score.strip() == category)
        for category in cerad_categories
    ]

    plt.bar(cerad_categories, counts, color=["skyblue", "lightgreen", "orange", "salmon"], edgecolor="black")
    plt.ylabel("Number of Patients")
    plt.title("CERAD Score Distribution (All Categories)")
    plt.show()

    # 6. FIXED ANOVA: Test Age at Death differences across CERAD categories
    cerad_groups = {
        category: [
            int(p.age_at_death) for p in Patient.all_patients
            if p.CERAD_score.strip() == category and str(p.age_at_death).isdigit()
        ]
        for category in cerad_categories
    }

    # Keep only non-empty groups
    non_empty_groups = [values for values in cerad_groups.values() if values]

    if len(non_empty_groups) > 1:
        f_stat, p_val = stats.f_oneway(*non_empty_groups)
        print("\nANOVA Results (Age at Death across CERAD categories):")
        print(f"F-value: {f_stat:.3f}")
        print(f"P-value: {p_val:.3f}")

        if p_val < 0.05:
            print("Result: Statistically significant difference between groups (p < 0.05).")
        else:
            print("Result: No statistically significant difference (p >= 0.05).")
    else:
        print("Not enough groups with age data to run ANOVA.")

    # 7. Scatter Plot: APOE vs CERAD scores (ordered labels)
    import random
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats

    apoe_positions = {"2/3": 1, "3/3": 2, "3/4": 3, "4/4": 4}
    x, y = [], []

    # Collect scatter data
    for p in Patient.all_patients:
        genotype = p.apoe.strip().upper()
        if str(p.age_at_death).isdigit() and genotype in apoe_positions:
            # jitter x a little so points don't overlap exactly
            x.append(apoe_positions[genotype] + random.uniform(-0.1, 0.1))
            y.append(int(p.age_at_death))

    # Scatter plot
    plt.figure(figsize=(8, 6))
    plt.scatter(x, y, alpha=0.7)
    plt.xticks(list(apoe_positions.values()), list(apoe_positions.keys()))
    plt.xlabel("APOE Genotype")
    plt.ylabel("Age at Death")
    plt.title("Scatter: APOE Genotype vs Age at Death")
    plt.ylim(70, 105)  # focus on your valid range
    plt.grid(axis="y", linestyle="--", alpha=0.6)

    # Overlay mean ± std dev
    for genotype, xpos in apoe_positions.items():
        ages = [
            int(p.age_at_death) 
            for p in Patient.all_patients 
            if str(p.age_at_death).isdigit() and p.apoe.strip().upper() == genotype
        ]
        if ages:
            mean_age = np.mean(ages)
            std_age = np.std(ages)
            # mean
            plt.scatter([xpos], [mean_age], color="red", s=120, marker="_", linewidths=2)
            # std error bars
            plt.errorbar(xpos, mean_age, yerr=std_age, color="red", capsize=5, linestyle="none")

    plt.show()

    # Run ANOVA test across groups
    groups = [
        [int(p.age_at_death) for p in Patient.all_patients 
        if str(p.age_at_death).isdigit() and p.apoe.strip().upper() == genotype]
        for genotype in apoe_positions.keys()
    ]
    groups = [g for g in groups if g]  # remove empty groups

    if len(groups) > 1:
        f_stat, p_val = stats.f_oneway(*groups)
        print("\nANOVA Results (Age at Death across APOE groups):")
        print(f"F = {f_stat:.3f}, p = {p_val:.3f}")







    
    #9) GET PATIENT ATTRIBUTES THAT WE WANT TO COMPARE ON A SCATTER PLOT
    #death_age_list = []
    #ABeta42 = []

    #for patient in Patient.all_patients:
        #death_age_list.append(patient.death_age)

    #for patient in Patient.all_patients:
        #ABeta42.append(patient.ABeta42)

    #X = [death_age_list] # Independent variable
    #y = [ABeta42] # Dependent variable
    
    #10) VISUALIZE DATA ON A SCATTER PLOT
    #plt.scatter(X, y, color='blue')
    #plt.xlabel('Age of Death')
    #plt.ylabel('ABeta42')
    #plt.title('Scatter Plot of Age of Death vs ABeta42')
    #plt.show()
