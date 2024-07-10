# open the output log file
file = open("OUTLOG.runS")

n = 0                           # Count lines
mean = [0.0, 0.0, 0.0]          # Mean positions

# Open the output log file
#
while (line := file.readline()) != "":
    if n >= 6:                  # Skip the header stuff
        v = [float(x) for x in line.split('|')]
        mean[0] += v[3]         # x pos
        mean[1] += v[4]         # y pos
        mean[2] += v[5]         # z pos
    n = n + 1                   # Count lines

if n>6:                         # Sanity check
    for i in range(3):
        mean[i] = mean[i]/(n-6) - 0.5

# If the squared values are close to 0.0, assume it worked
#
if mean[0]*mean[0] > 0.03 or mean[1]*mean[1] > 0.03 or mean[2]*mean[2] > 0.03:
    exit(1)
else:
    exit(0)
