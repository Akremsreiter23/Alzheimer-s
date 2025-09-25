import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

class Patient:
    all_patients = []

    def __init__(self, DonorID, APOE_Genotype, CERAD_score, CognitiveStatus, Age_at_death):
        self.donor_id = DonorID
        self.apoe = APOE_Genotype

        # Normalize CERAD (remove spaces, fix capitalization)
        if CERAD_score:
            self.CERAD_score = CERAD_score.strip().title()
        else:
            self.CERAD_score = ""       

        self.cognitive_status = CognitiveStatus
        self.age_at_death = Age_at_death
        Patient.all_patients.append(self)

    def __repr__(self):
        return (
            f"Donor ID: {self.donor_id} | "
            f"APOE Genotype: {self.apoe} | "
            f"CERAD score: {self.CERAD_score} | "
            f"Cognitive Status: {self.cognitive_status} | "
            f"Age at death: {self.age_at_death}"
    )


    @classmethod
    def instantiate_from_csv(cls, filename: str):
        # Use utf-8-sig to handle BOM
        with open(filename, encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            
            # Normalize headers: strip spaces and lowercase them
            rows_of_patients = []
            for row in reader:
                normalized_row = {k.strip().lower(): v for k, v in row.items()}
                rows_of_patients.append(normalized_row)

            # Create Patient objects
            for row in rows_of_patients:
                Patient(
                    DonorID=row["donor id"],
                    APOE_Genotype=row["apoe genotype"],
                    CERAD_score=row["cerad score"],
                    CognitiveStatus=row["cognitive status"],
                    Age_at_death=row["age at death"]
                )
