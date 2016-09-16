#ifndef SIMULATOR_SIMULATOR_HPP
#define SIMULATOR_SIMULATOR_HPP

#include "Node.hpp"
#include "Queue.hpp"
#include <string.h>
#include <iostream>
#include <fstream>
#include <pthread.h>
#include <unistd.h>
#include <iterator>
#include <sstream>
#include <array>

class Simulator {
private:
    bool simulating = true;
    std::thread simulatorThread;

    unsigned int sleepTime = 1;
    unsigned int numDestinations = 0;
    std::ofstream out;
    Queue<std::string> queue;
    std::string logString;
    Node* nodes;
    int nodeCount;
    unsigned int currentTick;
    std::vector<Packet> unaddedPackets; // Must be presorted
    int packetIndex = 0;
    bool tickWasActive = false;
    int countActiveTicks = 0;

    void runTick(); // Packets that finish transmitting during this tick will be added to Simulator::transmittedPackets

public:
    Simulator(Node* nodes, int nodeCount, std::vector<Packet> & packets,std::string logFilePath);
    std::queue<std::pair<unsigned short,Packet>> transmittedPackets;

    ~Simulator();

    void handler();
    void log(std::string logString);
    void log(std::vector<std::string> logVector);

    void start(bool withNetworkCoding = true); 

};


#endif //SIMULATOR_SIMULATOR_HPP
