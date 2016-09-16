
numTickString = "Total number of CTS: "
numTickPacketString = "Total number of RTS: "
numTicksCoding = 0
numTicksNoCoding = 0
numTickPacketCoding = 0
numTickPacketNoCoding = 0
writeData = open('tickdata.csv', 'w')
for x in range(1,1000):
	filename = "smallMesh/smallMeshLog" + str(x)
	dataFile = open(filename, 'r')
	first = True
	for line in dataFile:
		if numTickString in line and first:
			numTicksCoding += int(line[len(numTickString):])
		elif numTickString in line:
			numTicksNoCoding += int(line[len(numTickString):])
		elif numTickPacketString in line and first:
			numTickPacketCoding += int(line[len(numTickPacketString):])
			first = False
		elif numTickPacketString in line:
			numTickPacketNoCoding += int(line[len(numTickPacketString):])
	dataFile.close()
avgTicksCoding = numTicksCoding / 1000
avgTicksNoCoding = numTicksNoCoding / 1000
avgTickPacketCoding = float(numTickPacketCoding / 1000)
avgTickPacketNoCoding = float(numTickPacketNoCoding / 1000)

avgEfficiencyCoding = avgTickPacketCoding / avgTicksCoding
avgEfficiencyNoCoding = avgTickPacketNoCoding / avgTicksNoCoding
print avgTickPacketNoCoding
print avgTicksNoCoding
writeData.write("10, " + str(avgEfficiencyCoding) + ", " + '10' + ", " + str(avgEfficiencyNoCoding) + "\n")

numTicksCoding = 0
numTicksNoCoding = 0
numTickPacketCoding = 0
numTickPacketNoCoding = 0
for x in range(1,1000):
	filename = "medMesh/medMeshLog" + str(x)
	dataFile = open(filename, 'r')
	first = True
	for line in dataFile:
		if numTickString in line and first:
			numTicksCoding += int(line[len(numTickString):])
		elif numTickString in line:
			numTicksNoCoding += int(line[len(numTickString):])
		elif numTickPacketString in line and first:
			numTickPacketCoding += int(line[len(numTickPacketString):])
			first = False
		elif numTickPacketString in line:
			numTickPacketNoCoding += int(line[len(numTickPacketString):])
	dataFile.close()
avgTicksCoding = numTicksCoding / 1000
avgTicksNoCoding = numTicksNoCoding / 1000
avgTickPacketCoding = float(numTickPacketCoding / 1000)
avgTickPacketNoCoding = float(numTickPacketNoCoding / 1000)

print avgTickPacketNoCoding
print avgTicksNoCoding
avgEfficiencyCoding = avgTickPacketCoding / avgTicksCoding
avgEfficiencyNoCoding = avgTickPacketNoCoding / avgTicksNoCoding
writeData.write("25, " + str(avgEfficiencyCoding) + ", " + '25' + ", " + str(avgEfficiencyNoCoding) + "\n")


numTicksCoding = 0
numTicksNoCoding = 0
numTickPacketCoding = 0
numTickPacketNoCoding = 0
for x in range(1,1000):
	filename = "largeMesh/largeMeshLog" + str(x)
	dataFile = open(filename, 'r')
	first = True
	for line in dataFile:
		if numTickString in line and first:
			numTicksCoding += int(line[len(numTickString):])
		elif numTickString in line:
			numTicksNoCoding += int(line[len(numTickString):])
		elif numTickPacketString in line and first:
			numTickPacketCoding += int(line[len(numTickPacketString):])
			first = False
		elif numTickPacketString in line:
			numTickPacketNoCoding += int(line[len(numTickPacketString):])
	dataFile.close()
avgTicksCoding = numTicksCoding / 1000
avgTicksNoCoding = numTicksNoCoding / 1000
avgTickPacketCoding = float(numTickPacketCoding / 1000)
avgTickPacketNoCoding = float(numTickPacketNoCoding / 1000)

print avgTickPacketNoCoding
print avgTicksNoCoding
avgEfficiencyCoding = avgTickPacketCoding / avgTicksCoding
avgEfficiencyNoCoding = avgTickPacketNoCoding / avgTicksNoCoding
writeData.write("75, " + str(avgEfficiencyCoding) + ", " + '75' + ", " + str(avgEfficiencyNoCoding) + "\n")

