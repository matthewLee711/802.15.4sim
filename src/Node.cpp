#include <algorithm>
#include <climits>
#include <stack>

#include "Node.hpp"

unsigned short Node::sequenceID = 0;

unsigned int Node::countCTS = 0;

unsigned int Node::countRTS = 0;

unsigned int Node::countReceiveRTS = 0;

unsigned int Node::countReceiveCTS = 0;

void Node::sendPacket(const Packet & packet, const unsigned int &tick) {
    if( packet.getDestination().size() == 1){
        this->routingTable[*packet.getDestination().begin()]->receivePacket(packet, tick);
    }
    else{
        std::multimap<Node*, unsigned short> data;

        // Place destinations in the map
        for( auto &d : packet.getDestination())
            data.emplace(this->routingTable[d], d);

        Packet temp = packet;

        // Send packets
        for( auto &n : data){
            // Send modified version of the packet with only 1 destination
            temp.setDestination(n.second);
            // Route the packet
            n.first->receivePacket(temp, tick);
            this->numSends++;
        }
    }
}

// Needed to create array of <Node> , 0 is reserved to single blank node
Node::Node() { this->uniqueID = ++(Node::sequenceID); }

Node::Node(unsigned short uniqueID)
        : uniqueID{uniqueID}
        { }

int Node::getNumPacketsSent() {
    return numSends;
}

void Node::setNeighbors(std::set<Node *> &neighbors) {
    this->neighbors = neighbors;
}

void Node::receivePacket(const Packet &packet, const unsigned int &tick) {
    this->inputBuffer.push(packet);

    if(this->lastTickActed < tick) {
        this->queueCount++;

        if (packet.isHighPriority()) {
            if (this->alternateDelayTick < 0 || this->alternateDelayTick > tick)
                this->alternateDelayTick = tick;
        }
        else {
            if (this->alternateDelayTick < 0 || this->alternateDelayTick > tick + Node::MAX_DELAY_FOR_LOW_PRIORITY)
                this->alternateDelayTick = tick + Node::MAX_DELAY_FOR_LOW_PRIORITY;
        }
    }
    else{
        if(packet.isHighPriority()) {
            if( this->queueDelayTick < 0 || this->queueDelayTick > tick)
                this->queueDelayTick = tick;
        }
        else{
            if( this->queueDelayTick < 0 || this->queueDelayTick > tick + Node::MAX_DELAY_FOR_LOW_PRIORITY)
                this->queueDelayTick = tick + Node::MAX_DELAY_FOR_LOW_PRIORITY;
        }
    }
}

void Node::queuePacket(const Packet &p, const unsigned int & tick) {
    if(p.isHighPriority()) {
        if( this->queueDelayTick < 0 || this->queueDelayTick > tick)
            this->queueDelayTick = tick;
    }
    else{
        if( this->queueDelayTick < 0 || this->queueDelayTick > tick + Node::MAX_DELAY_FOR_LOW_PRIORITY)
            this->queueDelayTick = tick + Node::MAX_DELAY_FOR_LOW_PRIORITY;
    }

    this->outputBuffer.push_back(p);

    if(this->lastSuccessfulRTSTick != 0)
        this->outQueueCount++;
}

void Node::slotAction(const unsigned int &tick, std::queue<std::pair<unsigned short,Packet>> & transmittedPackets, bool & tickWasActive) {
    ////////////////// READ //////////////////

    Packet temp;
    while( this->inputBuffer.size() - this->queueCount > 0 ) {
        temp = this->inputBuffer.front();
        this->inputBuffer.pop();

        if(temp.findAndRemove(this->uniqueID))
            transmittedPackets.push(std::pair<unsigned short,Packet>(this->uniqueID,temp));

        if(temp.getDestination().size() > 0)
            this->queuePacket(temp, tick);
    }

    // Apply alternate tick if needed
    if(this->alternateDelayTick >= 0) {
        this->queueDelayTick = this->alternateDelayTick;
        this->alternateDelayTick = -1;
    }

    ////////////////// SEND + SCHEDULE //////////////////
    if( this->sourceIDRTS != 0 && collisionTick -1 < tick) {
        // std::cout << "EMITING CTS FROM NODE: " << uniqueID << " TO NODE: " << this->sourceIDRTS << " ON TICK: " << tick << std::endl;
        this->emitCTS(this->sourceIDRTS, tick, tickWasActive);
        this->sourceIDRTS = 0;
    }
    else if( this->canSend ){
        // std::cout << "SENDING PACKET FROM NODE: " << uniqueID << " ON TICK: " << tick << std::endl;
        if(this->outputBuffer.size() - this->outQueueCount > 0) {
            this->sendPacket(this->outputBuffer.front(), tick);
            this->outputBuffer.erase(this->outputBuffer.begin());
        }

        if(Node::NETWORK_CODING) {
            while(this->outputBuffer.size() - this->outQueueCount > 0) {
                this->sendPacket(this->outputBuffer.front(), tick);
                this->outputBuffer.erase(this->outputBuffer.begin());
            }
        }

        this->outQueueCount = 0;

        this->canSend = false;
        this->queueDelayTick = -1;
        this->lastSuccessfulRTSTick = 0;
        this->expCounter = 1;
    }
    else if( this->sentRTS && backoffCounter == 0) {
        backoffCounter = rand() % 2u << expCounter;
        if(expCounter < 10)
            expCounter++;
        sentRTS = false;
    }
    else if( this->backoffCounter == 0 ){
        if( this->outputBuffer.size() - this->outQueueCount > 0 &&
                (this->queueDelayTick < 0 || this->queueDelayTick < tick) ) {
            std::set<unsigned short> tempset;

            // Union all recipients of the buffer
            for( auto itr = this->outputBuffer.begin(); itr + this->outQueueCount < this->outputBuffer.end(); itr++)
                std::set_union(
                        tempset.begin(),
                        tempset.end(),
                        (*itr).getDestination().begin(),
                        (*itr).getDestination().end(),
                        std::inserter(tempset, tempset.begin())
                );

            std::set<unsigned short> tempset2;

            for(auto itr = tempset.begin(); itr != tempset.end(); itr++)
                tempset2.insert(routingTable[*itr]->uniqueID);

            //std::cout << "EMITING RTS FROM NODE: " << uniqueID << " TO NODE(s): ";
            //for(auto itr = tempset2.begin(); itr != tempset2.end(); itr++)
              //  std::cout << *itr << " ";
            //std::cout << "ON TICK: " << tick << std::endl;
            this->emitRTS(this->uniqueID, tempset2, tick, tickWasActive);
        }
    }
    else{
        this->backoffCounter--;
    }

    this->lastTickActed = tick;
    this->queueCount = 0;
}

void Node::emitCTS(unsigned short sourceID, const unsigned int &tick, bool & tickWasActive) {
    Node::countCTS++;
    tickWasActive = true;
    //call receive cts on all neighbors
    for(auto &n : neighbors) {
        n->receiveCTS(sourceID, tick);
    }
    backoffCounter++;
}

void Node::emitRTS(unsigned short sourceID, std::set<unsigned short> destinationID, const unsigned int &tick, bool & tickWasActive) {
    Node::countRTS++;
    tickWasActive = true;
    //call receive rts on all neighbors
    for(auto &n : neighbors) {
        n->receiveRTS(sourceID, destinationID, tick);
    }
    sentRTS = true;
    backoffCounter++;
}

void Node::receiveCTS(unsigned short rstSourceID, const unsigned int &tick) {
    if(rstSourceID == uniqueID) {
        canSend = true;
        lastSuccessfulRTSTick = tick - 1;
        sentRTS = false;
        Node::countReceiveCTS++;
    }
    else if(this->lastTickActed < tick && lastTickCTSDelay < tick) {
        this->backoffCounter++;
        this->lastTickCTSDelay = tick;
    }
}

void Node::receiveRTS(unsigned short sourceID, std::set<unsigned short> destinationID, const unsigned int &tick) {
   //if(collisionTick == tick) {
        if(collisionTick < tick && sourceIDRTS == 0 && !sentRTS) {
            if(destinationID.find(uniqueID) != destinationID.end()) {
                sourceIDRTS = sourceID;
                Node::countReceiveRTS++;
            }
            else {
                if(this->lastTickActed < tick) {
                    delayEmitCTS = tick + 1;
                }
                if( lastTickReceivedRTS < tick ) {
                    backoffCounter++;
                    lastTickReceivedRTS = tick;
                }
                else{
                    //std::cout << "RTS COLLISION A" << std::endl;
                    collisionTick = tick;
                    sourceIDRTS = 0;
                }

            }
        }
        else {
            //std::cout << "RTS COLLISION B ON NODE: " << uniqueID <<
              //  " WITH RTS FROM NODE " << sourceID << std::endl;
            collisionTick = tick;
            sourceIDRTS = 0;
        }
    //}
}

//written by Eric Smith
//Dijkstra's Algorithm
void Node::buildRoutes() {
    std::set<Node*> allNodes = buildTopology();

    std::unordered_map<Node*, std::tuple<unsigned short, unsigned int, Node*> >initialRouting;

    //initialize routingTable
    for(auto itr = allNodes.begin(); itr != allNodes.end(); itr++) {
        //destination, # hops, 1st hop
        std::tuple<unsigned short, unsigned int, Node*> tempTuple((*itr)->uniqueID, UINT_MAX, *itr);
        initialRouting.insert({*itr,tempTuple});
    }

    initialRouting.at(this) = std::make_tuple(this->uniqueID, 0, this); //set this as root node

    std::set<Node*> pathKnown; //list of nodes whose least cost path is known
    pathKnown.insert(this);

    //for all nodes adjacent to root (this), set their distances to 1
    for(auto itr = allNodes.begin(); itr != allNodes.end(); itr++){
        Node* v = *itr;

        //current node is adjacent to root (this)
        if( this->neighbors.find(v) != this->neighbors.end() )
            initialRouting.at(v) = std::make_tuple(v->uniqueID, 1, v);
    }

    //loop until all least cost paths are known
    while( pathKnown.size() != allNodes.size() ){
        Node* w = nullptr;
        unsigned int leastCost = UINT_MAX;

        //find a node not in pathKnown with the minimum current cost
        for(auto itr = allNodes.begin(); itr != allNodes.end(); itr++){
            std::tuple<unsigned short, unsigned int, Node*> tempTuple = initialRouting.at( *itr );
            unsigned int curValue = std::get<1>(tempTuple);

            if( pathKnown.find(*itr) == pathKnown.end() && curValue < leastCost ){
                w = *itr;
                leastCost = curValue;
            }
        }

        //add w to list of known shortest paths
        pathKnown.insert(w);

        //update cost of path for all v adjacent to w and not in pathKnown
        for(auto itr = w->neighbors.begin(); itr != w->neighbors.end(); itr++){
            Node* v = *itr;

            if( pathKnown.find(v) == pathKnown.end() ){
                std::tuple<unsigned short, unsigned int, Node*> vTuple = initialRouting.at(v);
                unsigned int vDist = std::get<1>(vTuple);

                std::tuple<unsigned short, unsigned int, Node*> wTuple = initialRouting.at(w);
                unsigned int wDist = std::get<1>(wTuple);

                unsigned int minDist = std::min(vDist, (wDist + 1));

                Node* firstHop = std::get<2>(vTuple);

                //if this node is beyond depth 1, go back through nodes until firstHop is found
                //firstHop is a neighbor of the root (this)
                if( (wDist + 1) < vDist ){

                    Node* targetFirstHop = w;
                    while( this->neighbors.find(targetFirstHop) == this->neighbors.end() ){
                        std::tuple<unsigned short, unsigned int, Node*> target = initialRouting.at(targetFirstHop);
                        targetFirstHop = std::get<2>(target);
                    }

                    firstHop = targetFirstHop;
                }

                initialRouting.at(v) = std::make_tuple(std::get<0>(vTuple), minDist, firstHop);
            }
        }
    }

    //populate routingTable
    for(auto& x : initialRouting ){
        std::tuple<unsigned short, unsigned int, Node*> tempT = x.second;
        //int dist = std::get<1>(tempT); //distance
        routingTable.insert({x.first->uniqueID, std::get<2>(tempT)});
    }
}

//BFS to find all the nodes in the network
std::set<Node*> Node::buildTopology(){
    std::queue<Node*> frontier;

    std::set<Node*> allNodes;

    //mark root visited
    allNodes.insert(this);
    frontier.push(this);

    while( !frontier.empty() ){
        for(auto &n : frontier.front()->neighbors ){
            if(allNodes.insert(n).second){
                frontier.push(n);
            }
        }
        frontier.pop();
    }

    return allNodes;
}

void Node::printRoutingTable(){
    std::cout << " ---- Routing Table ----" << std::endl;
    std::cout << "| Dest ID | First Hop |" << std::endl;
    std::cout << "------------------------" << std::endl;
    for(auto& x: routingTable){
        if(x.first == this->uniqueID)
            std::cout << "  " << "root" << "\t|\t" << x.second->uniqueID << std::endl;
        else
            std::cout << "  " << x.first << "\t\t|\t" << x.second->uniqueID << std::endl;
    }
    std::cout << "------------------------" << std::endl;

}

unsigned short Node::getUniqueID() {
    return uniqueID;
}
