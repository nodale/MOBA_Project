################################################################################################
### Computational Methods for Operator-Based Analysis
#Technical University of Munich, Professorship for Thermofluid Dynamics | Pr. Polifke, Ph.D
#
#Created: 06/2024 | P. Brokof, G. Varillon
################################################################################################
import numpy as np
from matplotlib import pyplot as plt
import csv

#######################################################
# Plot stability results
#######################################################
Re          = np.array([40,      50,       60,       70])
growth_rate = np.array([ 0,      0,         0,        0]) #<-- enter your results here.


# Load reference data.
Re_ref = []
growth_rate_ref = []
with open('reference.csv', newline ='') as file:
    reader = csv.reader(file, delimiter = '\t')
    for row in reader:
        Re_ref.append(float(row[0]))
        growth_rate_ref.append(float(row[1]))
file.close()
Re_ref          = np.array(Re_ref)
growth_rate_ref = np.array(growth_rate_ref)

# Plot against reference.
plt.plot(Re,growth_rate, "bx")
plt.plot(Re_ref, growth_rate_ref, "k-", label="reference")
plt.plot([30, 120],[0, 0], "r--")
plt.xlabel("Re")
plt.ylabel("Growth rate")
plt.legend()
plt.title("Base flow stability")
plt.savefig("stability.svg")