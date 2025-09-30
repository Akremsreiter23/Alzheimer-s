import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from Patient import Patient 

if __name__ == "__main__":
#Load patients from UpdatedMetaData.csv
    Patient.all_patients.clear()  # avoid duplicates if re-run
    Patient.instantiate_from_csv("UpdatedMetaData.csv")

#Print all patients to verify
    print("Loaded patients:")
    for patient in Patient.all_patients:
        print(patient)

#Count patients by APOE genotype
    genotype_counts = {}
    for patient in Patient.all_patients:
        genotype_counts[patient.apoe] = genotype_counts.get(patient.apoe, 0) + 1

    print("\nAPOE Genotype counts:")
    for genotype, count in genotype_counts.items():
        print(f"{genotype}: {count}")

#Create CERAD histogram (numeric only)
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

    #Create Bar graph (CERAD scores vs number of patients)
    cerad_categories = ["Absent", "Sparse", "Moderate", "Frequent"]

    counts = [
        sum(1 for p in Patient.all_patients if p.CERAD_score.strip() == category)
        for category in cerad_categories
    ]

    plt.bar(cerad_categories, counts, color=["skyblue", "lightgreen", "orange", "salmon"], edgecolor="black")
    plt.ylabel("Number of Patients")
    plt.title("CERAD Score Distribution (All Categories)")
    plt.show()

#Making a two-way ANOVA table
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols

data = pd.DataFrame([{
    "DonorID": p.donor_id,
    "APOE": p.apoe,
    "CERAD": p.CERAD_score,
    "AgeAtDeath": float(p.age_at_death) if p.age_at_death else None
    } for p in Patient.all_patients])

# Drop rows with missing values (important for ANOVA)
data = data.dropna(subset=["APOE", "CERAD", "AgeAtDeath"])

# 3. Make sure categorical variables are categorical
data["APOE"] = data["APOE"].astype("category")
data["CERAD"] = data["CERAD"].astype("category")

# 4. Fit two-way ANOVA model
model = ols('AgeAtDeath ~ C(APOE) * C(CERAD)', data=data).fit()
anova_table = sm.stats.anova_lm(model, typ=2)

print("\nTwo-Way ANOVA Results:")
print(anova_table)

#Creating Linear Regression: Age at Death vs CERAD scores
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

cerad_numeric_map = {"Absent": 0, "Sparse": 1, "Moderate": 2, "Frequent": 3}

x_cerad, y_cerad = [], []
for p in Patient.all_patients:
    try:
        if p.CERAD_score in cerad_numeric_map and p.age_at_death:
            x_cerad.append(cerad_numeric_map[p.CERAD_score])
            y_cerad.append(float(p.age_at_death))
    except (ValueError, TypeError):
        continue

#Fit regression
slope, intercept, r_value, p_value, std_err = stats.linregress(x_cerad, y_cerad)

#Regression line
x_vals = np.array(sorted(set(x_cerad)))
y_vals = intercept + slope * x_vals

plt.figure(figsize=(8, 6))
plt.scatter(x_cerad, y_cerad, alpha=0.7, label="Patients")
plt.plot(x_vals, y_vals, color="red", label=f"Fit: y={slope:.2f}x+{intercept:.2f}")
plt.xticks(range(4), ["Absent", "Sparse", "Moderate", "Frequent"])
plt.xlabel("CERAD Score")
plt.ylabel("Age at Death")
plt.title("Linear Regression: Age at Death vs CERAD Score")
plt.ylim(70, 105)
plt.grid(axis="y", linestyle="--", alpha=0.6)
plt.legend()
plt.show()

print(f"Slope: {slope:.3f}")
print(f"Intercept: {intercept:.3f}")
print(f"R-squared: {r_value**2:.3f}")
print(f"P-value: {p_value:.3g}")

# 7. Create Scatter Plot: APOE vs CERAD scores (ordered labels)
#import random
#import numpy as np
#import matplotlib.pyplot as plt
#from scipy import stats

#apoe_positions = {"2/3": 1, "3/3": 2, "3/4": 3, "4/4": 4}
#x, y = [], []

# Collect scatter data
#for p in Patient.all_patients:
#        genotype = p.apoe.strip().upper()
#        if str(p.age_at_death) and genotype in apoe_positions:
#            # jitter x a little so points don't overlap exactly
#            x.append(apoe_positions[genotype] + random.uniform(-0.1, 0.1))
#            y.append(int(p.age_at_death))
#print(f"Collected {len(x)} data points")
#print("Unique genotypes in data:", set(p.apoe for p in Patient.all_patients))

# Collect scatter data (robust)
#apoe_positions = {"2/3": 1, "3/3": 2, "3/4": 3, "4/4": 4}
#x, y = [], []

#for p in Patient.all_patients:
    #genotype = p.apoe.strip()
    #try:
        #age = float(p.age_at_death)  # handles 85, 85.0, " 90", etc.
        #if genotype in apoe_positions:
            #x.append(apoe_positions[genotype] + random.uniform(-0.1, 0.1))
            #y.append(age)
    #except (ValueError, TypeError):
        #continue  # skip if age can't be converted


# Scatter plot
#plt.figure(figsize=(8, 6))
#plt.scatter(x, y, alpha=0.7)
#plt.xticks(list(apoe_positions.values()), list(apoe_positions.keys()))
#plt.xlabel("APOE Genotype")
#plt.ylabel("Age at Death")
#plt.title("Scatter: APOE Genotype vs Age at Death")
#plt.ylim(70, 105)  # focus on your valid range
#plt.grid(axis="y", linestyle="--", alpha=0.6)

    #Overlay mean ± std dev
#for genotype, xpos in apoe_positions.items():
        #ages = [
            #int(p.age_at_death) 
            #for p in Patient.all_patients 
            #if str(p.age_at_death).isdigit() and p.apoe.strip().upper() == genotype
        #]
        #if ages:
            #mean_age = np.mean(ages)
            #std_age = np.std(ages)
            # mean
            #plt.scatter([xpos], [mean_age], color="red", s=120, marker="_", linewidths=2)
            # std error bars
            #plt.errorbar(xpos, mean_age, yerr=std_age, color="red", capsize=5, linestyle="none")

#plt.show()









    
    #9) GET PATIENT ATTRIBUTES THAT WE WANT TO COMPARE ON A SCATTER PLOT
#death_age_list = []
#CERAD_result = []

#for patient in Patient.all_patients:
        #death_age_list.append(patient.age_at_death)

#for patient in Patient.all_patients:
        #CERAD_result.append(patient.CERAD_score)

#x = [death_age_list] # Independent variable
#y = [CERAD_result] # Dependent variable
    
    #10 VISUALIZE DATA ON A SCATTER PLOT
#plt.scatter(x, y, color='blue')
#plt.xlabel('Age of Death')
#plt.ylabel('ABeta42')
#plt.title('Scatter Plot of Age of Death vs ABeta42')
#plt.show()
