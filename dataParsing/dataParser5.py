tickArriveString = "has reached destination"
numTicksCoding = 0
numTicksNoCoding = 0

counterCoding = 0
counterNoCoding = 0
writeData = open('tickdata.csv', 'w')
for x in range(1,1000):
	filename = "smallMesh/smallMeshLog" + str(x)
	dataFile = open(filename, 'r')
	first = True
	for line in dataFile:
		if tickArriveString in line and first:
			start = line.find(' in ')
			end = line.find(' ticks')
			numTicksCoding += int(line[start+4:end])
			counterCoding+=1
		elif tickArriveString in line:
			start = line.find(' in ')
			end = line.find(' ticks')
			numTicksNoCoding += int(line[start+4:end])
			counterNoCoding+=1
		elif "End of simulator" in line:
			first = False
	dataFile.close()
avgTicksCoding = numTicksCoding / counterCoding
avgTicksNoCoding = numTicksNoCoding / counterNoCoding

writeData.write("10, " + str(counterCoding) + ", " + '10' + ", " + str(counterNoCoding) + "\n")

numTicksCoding = 0
numTicksNoCoding = 0

counterCoding = 0
counterNoCoding = 0
for x in range(1,1000):
	filename = "medMesh/medMeshLog" + str(x)
	dataFile = open(filename, 'r')
	first = True
	for line in dataFile:
		if tickArriveString in line and first:
			start = line.find(' in ')
			end = line.find(' ticks')
			numTicksCoding += int(line[start+4:end])
			counterCoding+=1
		elif tickArriveString in line:
			start = line.find(' in ')
			end = line.find(' ticks')
			numTicksNoCoding += int(line[start+4:end])
			counterNoCoding+=1
		elif "End of simulator" in line:
			first = False
	dataFile.close()
avgTicksCoding = numTicksCoding / counterCoding
avgTicksNoCoding = numTicksNoCoding / counterNoCoding

writeData.write("25, " + str(counterCoding) + ", " + '25' + ", " + str(counterNoCoding) + "\n")
numTicksCoding = 0
numTicksNoCoding = 0

counterCoding = 0
counterNoCoding = 0
for x in range(1,1000):
	filename = "largeMesh/largeMeshLog" + str(x)
	dataFile = open(filename, 'r')
	first = True
	for line in dataFile:
		if tickArriveString in line and first:
			start = line.find(' in ')
			end = line.find(' ticks')
			numTicksCoding += int(line[start+4:end])
			counterCoding+=1
		elif tickArriveString in line:
			start = line.find(' in ')
			end = line.find(' ticks')
			numTicksNoCoding += int(line[start+4:end])
			counterNoCoding+=1
		elif "End of simulator" in line:
			first = False
	dataFile.close()
avgTicksCoding = numTicksCoding / counterCoding
avgTicksNoCoding = numTicksNoCoding / counterNoCoding

writeData.write("75, " + str(counterCoding) + ", " + '75' + ", " + str(counterNoCoding) + "\n")
