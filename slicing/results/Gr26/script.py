count = dict()
count[-1] = 0
for i in range(7+1):
    count[2*i] = 0

print(count)

for t in range(10):
    f = open("./slicing/results/Gr26/a"+str(t)+".txt", "r")
    lines = f.read().splitlines()
    for line in lines:
        s = line.split(" ")
        if(s[1] == "0"):
            continue
        count[int(s[0])] += int(s[1])
    f.close()

b = 0
for c in count.values():
    b += c
b