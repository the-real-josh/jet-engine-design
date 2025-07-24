# find out the general stuff (no detail)
import matplotlib.pyplot as plt
import matplotlib.patches as patches


fig, ax = plt.subplots()
ax.add_patch(patches.Rectangle((1.0,1.0),1.0,1.0))
plt.show()