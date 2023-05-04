"""
This program shows an example of generating bar chart with amino acid frequencies as input
"""

import matplotlib.pyplot as plt

# Define the amino acid frequencies as a dictionary
# In later iteration the frequencies should be read into variables 
aa_frequencies = {'A': 0.08, 'R': 0.04, 'N': 0.03, 'D': 0.05,
                  'C': 0.02, 'Q': 0.03, 'E': 0.05, 'G': 0.07,
                  'H': 0.02, 'I': 0.05, 'L': 0.09, 'K': 0.05,
                  'M': 0.02, 'F': 0.04, 'P': 0.05, 'S': 0.06,
                  'T': 0.06, 'W': 0.01, 'Y': 0.03, 'V': 0.07}

# Define the x-axis and y-axis values
x_values = list(aa_frequencies.keys())
y_values = list(aa_frequencies.values())

# Create a bar chart using matplotlib
plt.bar(x_values, y_values)

# Set the title and axis labels
plt.title('Amino Acid Frequencies')
plt.xlabel('Amino Acid')
plt.ylabel('Frequency')

# Show the plot
plt.show()

