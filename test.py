import matplotlib.pyplot as plt
import matplotlib.patches as patches


fig, ax = plt.subplots()
ax.add_patch(patches.Polygon(([(0, 0),
                               (1, 0),
                               (1, 1),
                               (0, 1)])))
plt.title(f'1D stages (in meters)')
plt.show()
plt.clf()