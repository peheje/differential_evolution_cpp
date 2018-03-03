import matplotlib.pyplot as plt

with open("/Users/phj/Desktop/data.txt") as f:
    data = f.read()

data = data.split('\n')

x = [float(row.split(' ')[0]) for row in data if row != ""]
y = [float(row.split(' ')[1]) for row in data if row != ""]

fig = plt.figure()

ax1 = fig.add_subplot(111)
ax1.set_yscale('log')

ax1.set_title("Differential evolution")
ax1.set_xlabel('Generation')
ax1.set_ylabel('Score')

ax1.plot(x, y, c='r', label='the data')

leg = ax1.legend()

plt.show()
