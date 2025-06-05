import sys

# cutAbs = 1.0e-6
cutAbs = 0.05
cutRel = 0.05

if len(sys.argv) < 2:
    print("Использование: python script.py имя_файла")
    sys.exit(1)

filename = sys.argv[1]

with open(filename, "r") as file:
    for line in file:
        line = line.strip()
        if not line:
            continue

        parts = line.split()
        if len(parts) != 15:
            print(f"Ожидалось 15 полей, а найдено {len(parts)}: {parts}")
            continue

        index = parts[0]
        date = parts[2]
        time = parts[3]
        diffAbs = float(parts[11])
        diffRel = float(parts[12])
        name = parts[14]

        if diffAbs > cutAbs and diffRel > cutRel:
            print(index, date, time, diffAbs, diffRel, name)
