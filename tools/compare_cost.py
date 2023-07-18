import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

cost_true = [
    -147.0, datetime(2021, 7, 16,  9, 20),
    -100.9, datetime(2021, 7, 16,  9, 30),
    - 50.8, datetime(2021, 7, 16,  9, 40),
    - 30.8, datetime(2021, 7, 16,  9, 42),
    -  7.0, datetime(2021, 7, 16,  9, 46),
    -  5.0, datetime(2021, 7, 16,  9, 53),
    -  3.0, datetime(2021, 7, 16, 10, 13),
    -  1.0, datetime(2021, 7, 16, 10, 49),
       1.0, datetime(2021, 7, 16, 11, 34),
       3.0, datetime(2021, 7, 16, 12, 27),
       5.0, datetime(2021, 7, 16, 13, 30),
       7.0, datetime(2021, 7, 16, 14, 45),
       9.0, datetime(2021, 7, 16, 16,  7),
      11.0, datetime(2021, 7, 16, 17, 19),
      13.0, datetime(2021, 7, 16, 18, 17),
      15.0, datetime(2021, 7, 16, 19,  9),
      17.0, datetime(2021, 7, 16, 19, 56),
      19.0, datetime(2021, 7, 16, 20, 38),
      29.0, datetime(2021, 7, 16, 22, 45),
      39.0, datetime(2021, 7, 17,  0,  2),
      59.0, datetime(2021, 7, 17,  1, 46),
      89.0, datetime(2021, 7, 17,  3, 21),
     135.0, datetime(2021, 7, 17,  5,  8),
]

cost_skop = [
    -147.0, datetime(2021, 7, 19, 10, 31),
    -100.9, datetime(2021, 7, 19, 10, 40),
    - 50.8, datetime(2021, 7, 19, 10, 48),
    - 30.8, datetime(2021, 7, 19, 10, 50),
    -  7.0, datetime(2021, 7, 19, 10, 53),
    -  5.0, datetime(2021, 7, 19, 11, 00),
    -  3.0, datetime(2021, 7, 19, 11, 17),
       1.0, datetime(2021, 7, 19, 12, 28),
       3.0, datetime(2021, 7, 19, 13, 14),
       5.0, datetime(2021, 7, 19, 14,  7),
       7.0, datetime(2021, 7, 19, 15,  8),
]

def plot_cost(cost, name=""):
    c = np.array(cost).reshape((len(cost)>>1,2))
    x_arr = []
    y_arr = []
    for i in range(1,len(c)):
        x = (c[i][0] + c[i-1][0])/2
        y = (c[i][1] - c[i-1][1]).total_seconds() / (c[i][0] - c[i-1][0])
        x_arr.append(x)
        y_arr.append(y)
    plt.plot(x_arr, y_arr, ".-", label=name)
        



plot_cost(cost_true, name="Nordin")
plot_cost(cost_skop, name="Nordin_sk34op")


plt.xlabel("Deg")
plt.ylabel("Seconds per Deg")
plt.show()
