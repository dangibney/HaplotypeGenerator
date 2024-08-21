import matplotlib.pyplot as plt
from statistics import mean, mode

# Load the recombination counts from the file
with open("recombination_counts.txt", "r") as f:
    recombination_numbers = [int(line.strip()) for line in f]

# Get the unique values for the recombination numbers
unique_values = sorted(set(recombination_numbers))

# Calculate the mean and mode of recombination numbers
mean_recombinations = mean(recombination_numbers)
mode_recombinations = mode(recombination_numbers)

# Print the mean and mode
print(f"Mean number of recombinations: {mean_recombinations}")
print(f"Mode number of recombinations: {mode_recombinations}")

# Create bins where each bin is centered around the unique values
bins = [x - 0.5 for x in range(min(unique_values), max(unique_values) + 2)]

# Plot the histogram
plt.hist(recombination_numbers, bins=bins, edgecolor='black', rwidth=0.8)

# Set the x-ticks to be the unique recombination values
plt.xticks(unique_values)

# Add titles and labels
plt.title('Distribution of Number of Recombinations')
plt.xlabel('Number of Recombinations')
plt.ylabel('Frequency')

# Add a grid for better readability
plt.grid(True)

# Display the plot
plt.show()
