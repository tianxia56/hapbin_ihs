import pandas as pd
import matplotlib.pyplot as plt

# Function to plot demographic history from CSV
def plot_demographic_history_from_csv(file, label):
    df = pd.read_csv(file)
    plt.plot(df['x'], df['y'], label=label)

# Load data from CSV files
files = ["analysis_IND/IND.csv", "analysis_MGN/MGN.csv", "analysis_SASI/SASI.csv", "analysis_SASP/SASP.csv", "../reference/analysis_GBR/GBR.csv","../reference/analysis_YRI/YRI.csv","../reference/analysis_KOR/KOR.csv", "../reference/analysis_PJL/PJL.csv"]
labels = ["IND", "MGN", "SASI", "SASP", "GBR", "YRI", "KOR", "PJL"]

plt.figure(figsize=(10, 6))

for file, label in zip(files, labels):
    plot_demographic_history_from_csv(file, label)

# Adding labels and title
plt.xlabel("Time (in generations)")
plt.ylabel("Effective Population Size (Ne)")
plt.title("Demographic History")
plt.legend()

# Save the plot as a PNG file
plt.savefig("all.pop.demo.png")
