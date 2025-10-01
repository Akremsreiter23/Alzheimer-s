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