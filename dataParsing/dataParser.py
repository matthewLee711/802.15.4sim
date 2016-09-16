import sys

lncString = "with network coding"
noLncString = "Starting simulator without network coding"
tickString = "Total number of ticks: "
tickPacketString = "Number of ticks during which a packet reached a destination: "
ctsEmitString = "Total number of CTS Emit: "
rtsEmitString = "Total number of RTS Emit: "
ctsReceiveString = "Total number of CTS Receive: "
rtsReceiveString = "Total number of RTS Receive: "
endSimulatorString = "End of simulator"
numDestString = "Number of Destinations: "
latencyString = "has reached destination node"
activeTicksString = "Total number of active ticks: "

iterations = 1000		#Num of iterations of log loops

usingLnc = {}
noLnc = {}

usingLnc["totalTicks"] = []	#[value, # of times this value was taken]
usingLnc["tickPacket"] = []
usingLnc["ctsEmit"] = []
usingLnc["rtsEmit"] = []
usingLnc["ctsReceive"] = []
usingLnc["rtsReceive"] = []
usingLnc["numDest"] = []
usingLnc["latency"] = []
usingLnc["activeTicks"] = []

noLnc["totalTicks"] = []
noLnc["tickPacket"] = []
noLnc["ctsEmit"] = []
noLnc["rtsEmit"] = []
noLnc["ctsReceive"] = []
noLnc["rtsReceive"] = []
noLnc["numDest"] = []
noLnc["latency"] = []
noLnc["activeTicks"] = []

fileType = sys.argv[1]
for x in range(1, iterations+1):
	filename = fileType + str(x)
	#filename = "logfile" 	#for testing
	dataFile = open(filename, 'r')
	lnc = True

	for line in dataFile:
		if lncString in line:
			lnc = True
		elif noLncString in line:
			lnc = False
		if lnc:
			if tickString in line:
				usingLnc["totalTicks"].append(int(line[len(tickString):]))
			elif latencyString in line:
				usingLnc["latency"].append(int(line[line.find("in ")+3:line.find("ticks")]))
			elif numDestString in line:
				usingLnc["numDest"].append(int(line[len(numDestString):]))
			elif tickPacketString in line:
				usingLnc["tickPacket"].append(int(line[len(tickPacketString):]))
			elif ctsEmitString in line:
				usingLnc["ctsEmit"].append(int(line[len(ctsEmitString):]))
			elif rtsEmitString in line:
				usingLnc["rtsEmit"].append(int(line[len(rtsEmitString):]))
			elif ctsReceiveString in line:
				usingLnc["ctsReceive"].append(int(line[len(ctsReceiveString):]))
			elif rtsReceiveString in line:
				usingLnc["rtsReceive"].append(int(line[len(ctsReceiveString):]))
			elif activeTicksString in line:
				usingLnc["activeTicks"].append(int(line[len(activeTicksString):]))
			elif endSimulatorString in line:
					lnc = False
		else:
			if tickString in line:
				noLnc["totalTicks"].append(int(line[len(tickString):]))
			elif latencyString in line:
				noLnc["latency"].append(int(line[line.find("in ")+3:line.find("ticks")]))
			elif numDestString in line:
				noLnc["numDest"].append(int(line[len(numDestString):]))
			elif tickPacketString in line:
				noLnc["tickPacket"].append(int(line[len(tickPacketString):]))
			elif ctsEmitString in line:
				noLnc["ctsEmit"].append(int(line[len(ctsEmitString):]))
			elif rtsEmitString in line:
				noLnc["rtsEmit"].append(int(line[len(rtsEmitString):]))
			elif ctsReceiveString in line:
				noLnc["ctsReceive"].append(int(line[len(ctsReceiveString):]))
			elif rtsReceiveString in line:
				noLnc["rtsReceive"].append(int(line[len(ctsReceiveString):]))
			elif activeTicksString in line:
				noLnc["activeTicks"].append(int(line[len(activeTicksString):]))
			elif endSimulatorString in line:
				lnc = True
	dataFile.close()

mean = lambda MyList : reduce(lambda x, y: x + y, MyList) / float(len(MyList))
stdv = lambda MyList : (reduce(lambda x,y : x + y , map(lambda x: (x-mean(MyList))**2 , MyList)) / float(len(MyList)))**.5
writeFile = open(fileType + "Results.csv", 'w')
writeFile.write("Using LNC,Min,Max,Mean,STDV\n")
for x in usingLnc.keys():
	write =  x + "," + str(min(usingLnc[x])) + "," +  str(max(usingLnc[x])) + "," + str(mean(usingLnc[x])) + "," + str(stdv(usingLnc[x]))
	writeFile.write(write+"\n")

writeFile.write("\nNo LNC,Min,Max,Mean,STDV\n")
for x in noLnc.keys():
	write =  x + ", " + str(min(noLnc[x])) + "," + str(max(noLnc[x])) + "," +str(mean(noLnc[x])) + "," + str(stdv(noLnc[x]))
	writeFile.write(write+"\n")
writeFile.close()