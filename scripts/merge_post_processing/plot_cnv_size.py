import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Read data from your file
data = pd.read_csv('cnv_size_count_final.txt', sep='\t')

# Prepare data - convert size classes to match your example format
size_classes = [
    '1-2', '2-5', '5-10', '10-20', 
    '20-50', '50-100', '100+'
]
deletion = data['Deletion']
duplication = data['Duplication']
mixed = data['Mixed']

# Create figure and axis
fig, ax = plt.subplots(figsize=(10, 6))

# Set colors as requested (red, green, blue)
colors = {
    'Deletion': '#FF0000',      # Red
    'Duplication': '#00AA00',   # Green
    'Mixed': '#0000FF'          # Blue
}

# Plot stacked bars
p1 = ax.bar(size_classes, deletion, label='Deletion', color=colors['Deletion'])
p2 = ax.bar(size_classes, duplication, bottom=deletion, label='Duplication', color=colors['Duplication'])
p3 = ax.bar(size_classes, mixed, bottom=np.array(deletion)+np.array(duplication), label='Mixed', color=colors['Mixed'])

# Add labels and title
ax.set_xlabel('Size Class (kb)', fontsize=12)
ax.set_ylabel('Number of CNVRs', fontsize=12)
ax.set_title('CNVR Distribution by Size and Type', fontsize=14, pad=20)

# Customize ticks and grid
ax.tick_params(axis='both', which='major', labelsize=10)
ax.grid(axis='y', linestyle='--', alpha=0.7)

# Add legend at upper right
ax.legend(loc='upper right', framealpha=1)

# Add value labels on top of each segment
for i, (del_val, dup_val, mix_val) in enumerate(zip(deletion, duplication, mixed)):
    total = del_val + dup_val + mix_val
    if total > 0:  # Only label if there's data
        ax.text(i, total + 0.02*total, f'{total}', 
                ha='center', va='bottom', fontsize=8)

# Save the plot with high quality
plt.tight_layout()
plt.savefig('cnvr_distribution.png', dpi=300, bbox_inches='tight')
print("Plot saved successfully as 'cnvr_distribution.png'")