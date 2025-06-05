import sys

from dataclasses import dataclass

@dataclass
class peDissolReport:
    time: float
    timeAv: float
    numSurfaceCells: int
    numVolumeCells: int    

if len(sys.argv) < 2:
    print("Использование: python script.py имя_файла")
    sys.exit(1)

dissolReport = []

input_filename = sys.argv[1]
with open(input_filename, "r") as file:
    for line in file:
        line = line.strip()
        if not line:
            continue

        parts = line.split()
        if len(parts) != 4:
            print(f"Ожидалось 4 поля, а найдено {len(parts)}: {parts}")
            continue

        time = float(parts[1])
        numSurfaceCells = int(parts[2])
        numVolumeCells = int(parts[3])
        
        dissolReport.append(peDissolReport(time, 0.0, numSurfaceCells, numVolumeCells))

print(f"Количество записей: {len(dissolReport)}")

dissolReport.sort(key=lambda record: record.time)
print(f"Массив отсортирован по возрастанию времени")

totalTime = sum(record.time for record in dissolReport)
averageTime = totalTime / len(dissolReport)
minTime = dissolReport[0].time
maxTime = dissolReport[-1].time

print(f"   Суммарное время выполнения: {totalTime:>12.4f}")
print(f"     Среднее время выполнения: {averageTime:>12.4f}")
print(f" Минимальное время выполнения: {minTime:>12.4f}")
print(f"Максимальное время выполнения: {maxTime:>12.4f}")
print(f"        Показатель дисбаланса: {maxTime/averageTime:>12.4f}")

totalSurfaceCells = sum(record.numSurfaceCells for record in dissolReport)
averageSurfaceCells = int(totalSurfaceCells / len(dissolReport))
minSurfaceCells = min(record.numSurfaceCells for record in dissolReport)
maxSurfaceCells = max(record.numSurfaceCells for record in dissolReport)

print(f"   Суммарное число поверхностных ячеек: {totalSurfaceCells:>12}")
print(f"     Среднее число поверхностных ячеек: {averageSurfaceCells:>12}")
print(f" Минимальное число поверхностных ячеек: {minSurfaceCells:>12}")
print(f"Максимальное число поверхностных ячеек: {maxSurfaceCells:>12}")
print(f"                 Показатель дисбаланса: {float(maxSurfaceCells)/float(averageSurfaceCells):>12.4f}")

totalVolumeCells = sum(record.numVolumeCells for record in dissolReport)
averageVolumeCells = int(totalVolumeCells / len(dissolReport))
minVolumeCells = min(record.numVolumeCells for record in dissolReport)
maxVolumeCells = max(record.numVolumeCells for record in dissolReport)

print(f"   Суммарное число объемных ячеек: {totalVolumeCells:>12}")
print(f"     Среднее число объемных ячеек: {averageVolumeCells:>12}")
print(f" Минимальное число объемных ячеек: {minVolumeCells:>12}")
print(f"Максимальное число объемных ячеек: {maxVolumeCells:>12}")
print(f"            Показатель дисбаланса: {float(maxVolumeCells)/float(averageVolumeCells):>12.4f}")

for record in dissolReport:
    record.timeAv = averageTime    

output_filename = input_filename + ".report"
with open(output_filename, "w") as outfile:
    outfile.write("СТАТИСТИКА\n")
    for idx, record in enumerate(dissolReport, start=1):
        outfile.write(f"{idx:>4} {record.time:>8.4f} {record.timeAv:>8.4f} {record.numSurfaceCells:>8} {record.numVolumeCells:>8}\n")
    outfile.write("КОНЕЦ СТАТИСТИКИ\n")
print(f"Результаты записаны в файл: {output_filename}")
