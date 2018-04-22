import matplotlib.pyplot as plt

paths = ["/Users/phj/Desktop/data0.txt", "/Users/phj/Desktop/data1.txt"]

colors = ["r", "b", "g"]
plots = []

skip_x = 10

for color, path in zip(colors, paths):
    with open(path) as f:
        data = f.read()

    data = data.split('\n')
    label = data[0]
    data.pop(0)

    for _ in range(skip_x):
        data.pop(0)

    x = [float(row.split(' ')[0]) for row in data if row != ""]
    y = [float(row.split(' ')[1]) for row in data if row != ""]

    # print(x)
    # print(y)

    item, = plt.plot(x, y, c=color, label=label)
    plots.append(item)

plt.yscale('log')
# plt.legend(handles=plots)
plt.show()
