#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "matrix.hpp"
#include "Simulator.hpp"

Node* starGen(int numNodes);
Node* gridGen(int numNodes, int numColumns);
Node* meshGen(int numNodes, std::vector<std::vector<std::string>>& neighbors);

int main( int argc, char * argv[] ) {

    // Check that required command line args were supplied
    if( argc == 4 ){
        // Load config
        std::string configPath = argv[1];
        std::ifstream configFile;
        configFile.open(configPath, std::ios::in);

        std::string topologyType;        //type of topology
        int numNodes;                   //number of nodes
        getline(configFile, topologyType);
        configFile >> numNodes;
        configFile.ignore();
        Node* nodes;
        if(topologyType.compare("Grid") == 0) {
            int numCol;
            configFile >> numCol;       //get number of columns

            //CREATE TOPOLOGY HERE
            nodes = gridGen(numNodes,numCol);

        }
        else if(topologyType.compare("Star") == 0){
            nodes = starGen(numNodes);
        }
        else if(topologyType.compare("Mesh") == 0) {
            std::vector<std::vector<std::string>> neighbors(numNodes);    //vector for neighbors; ID used as index
            for(int i = 0; i < numNodes; i++) {
                std::string buffer;
                int nodeid;             //gets id for node
                configFile >> nodeid;
                configFile.ignore();
                getline(configFile, buffer);        //get whole line of neighbors for parsing
                std::stringstream stream(buffer);
                std::string temp;
                while(stream >> temp) {         //parse whole line and push to vector
                    neighbors[nodeid-1].push_back(temp);
                }
            }
            //CREATE TOPOLOGY HERE
             nodes = meshGen(numNodes, neighbors);
        }
        configFile.close();

        //only for mesh networks

        // Load messages
        //std::ifstream messageFile(argv[1], std::ios::in);
        std::string messagePath = argv[2];
        std::ifstream messageFile;
        messageFile.open(messagePath, std::ios::in);

        std::vector<Packet> packets;

        std::string line;
        unsigned short id = 1;
        while(getline(messageFile, line)) {
            unsigned short source;
            std::set<unsigned short> destinations;
            unsigned int tick;
            bool priority;

            int priorityInput;
            std::stringstream stream(line);
            stream >> source;
            stream.ignore();        //ignore space
            stream.ignore();        //ignore '['
            unsigned short temp;
            while(stream.peek() != ']') {
                stream >> temp;
                destinations.insert(temp);
            }
            stream.ignore();        //ignore space
            stream >> tick;
            stream.ignore();
            stream >> priorityInput;
            priority = (bool) priorityInput;
            packets.push_back(Packet(id, source, destinations, tick, priority));
            id++;


        }

        messageFile.close();
        // Create random network
        std::string logFilePath = argv[3];
        std::cout << "Simulation with network coding" << std::endl;
        Simulator simulation(nodes,numNodes,packets,logFilePath);
        simulation.start();
        std::cout << "Simulation without network coding (control)" << std::endl;
        Simulator simulation2(nodes,numNodes,packets,logFilePath);
        simulation2.start(false);
        // Create & start simulator
        //Simulator exampleSimulator;
        //exampleSimulator.log("testdata1,testdata2");
        return 0;
    }
    // Invalid command line arguments, alert user and halt program
    else{
        std::cout << "ERROR: INVALID COMMAND LINE ARGS..." << std::endl
                  << "       Proper Usage: <exec_cmd> <topology_conf_path> <message_file_path> <log_file_path>" << std::endl;
    }
}
