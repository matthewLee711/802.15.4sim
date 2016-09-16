#ifndef SIMULATOR_NODE_HPP
#define SIMULATOR_NODE_HPP

#include <vector>
#include <unordered_map>
#include <queue>
#include <set>
#include <map>
#include <iostream>
#include <algorithm>
#include <cstdint>

#include "Packet.hpp"

class Simualtor;
struct PacketComparison;

class Node {
private:

    unsigned int collisionTick = 0;
    unsigned int delayEmitCTS = 0;
    unsigned short sourceIDRTS = 0;
    unsigned short sourceIDCTS = 0;
    unsigned int backoffCounter = 0;
    uint8_t expCounter = 0;
    bool canSend = false;
    unsigned int outQueueCount = 0;
    unsigned int lastSuccessfulRTSTick = 0;
    int queueDelayTick = 0;
    int alternateDelayTick = 0;
    bool sentRTS = false;
    int lastTickCTSDelay = 0;
    int lastTickReceivedRTS = 0;

    static unsigned short sequenceID;
    unsigned short uniqueID;
    std::set<Node*> neighbors;
    std::queue<Packet> inputBuffer;
    std::vector<Packet> outputBuffer;
    std::unordered_map<unsigned short, Node*> routingTable;

    unsigned int queueCount; // The number of packets that have been added to the queue during this tick
                             // This is needed to prevent receiving and processing a packet in the same tick
                             // It will always be set to 0 at the end of Node::slotAction() and will only be incremented
                             // by Node::receivePacket() and only if Node::lastTickActed < tick

    unsigned int lastTickActed; // Last tick that the node acted on, always updated by Node::slotAction()

    void sendPacket(const Packet & packet, const unsigned int &tick);
    void emitCTS(unsigned short sourceID, const unsigned int & tick, bool & tickWasActive);
    void emitRTS(unsigned short sourceID, std::set<unsigned short> destinationID, const unsigned int &tick, bool & tickWasActive);
    std::set<Node*> buildTopology(); // fills allNodes to create topology of network

public:
    static int MAX_DELAY_FOR_LOW_PRIORITY;
    static bool NETWORK_CODING;
    //Counters for CTS and RTS
    static unsigned int countRTS;
    static unsigned int countCTS;
    static unsigned int countReceiveRTS;
    static unsigned int countReceiveCTS;
    unsigned int numSends = 0;//number of sends for specific node

    Node();
    Node(unsigned short uniqueID);

    int getNumPacketsSent();

    void setNeighbors(std::set<Node*> & neighbors);
    void receivePacket(const Packet & packet, const unsigned int & tick); // Called by neighbor nodes when they send a packet
    void queuePacket(const Packet & p, const unsigned int & tick); // Called by simulator when a packet is "created" for the node to send
    void slotAction(const unsigned int & tick, std::queue<std::pair<unsigned short,Packet>> & transmittedPackets, bool & tickWasActive);
    unsigned short getUniqueID();
    // Called by simulator to run the node's actions during the current time slot (tick)
    void receiveRTS(unsigned short sourceID, std::set<unsigned short> destinationID, const unsigned int &tick);
    void receiveCTS(unsigned short rstSourceID, const unsigned int & tick);
    void buildRoutes(); // Use Dijkstra's algorithm to build the routing table
    void printRoutingTable();
};

#endif //SIMULATOR_NODE_HPP
