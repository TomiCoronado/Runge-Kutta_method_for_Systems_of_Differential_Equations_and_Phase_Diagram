import pandas as pd
from matplotlib import pyplot as plt 

output_P3 = pd.read_csv("output_P3.csv")
a=int(output_P3.iloc[-1,0])

for i in range(a):
    plt.plot(output_P3.iloc[2*i],output_P3.iloc[2*i+1], linewidth = 0.5)

plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.title("Phase Diagram")
plt.xlabel("y1")
plt.ylabel("y2")
plt.show()


