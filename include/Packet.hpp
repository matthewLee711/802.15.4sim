#ifndef SIMULATOR_PACKET_HPP
#define SIMULATOR_PACKET_HPP

#include <cstdint>
#include <vector>
#include <set>

class Packet {
private:

    unsigned int uniqueID;
    unsigned int creationTick;
    unsigned short source; // Unique ID of the source
    std::set<unsigned short> destination; // Unique ID of the destination
    bool highPriority; // Signals message must be sent ASAP (true) or it can be delayed for linear combination (false)

public:

    Packet();
    Packet(unsigned int uniqueID, unsigned short source, std::set<unsigned short> destination, unsigned int creationTick,
           bool highPriority=false);

    unsigned int getUniqueID();
    unsigned int getCreationTick() const;
    unsigned short getSource();
    const std::set<unsigned short>& getDestination() const;
    bool getPriority() const;

    bool isHighPriority() const;
    bool isLowPriority();

    void setUniqueID(unsigned int uniqueID);
    void setSource(unsigned short source);
    void setDestination(std::set<unsigned short> destination);
    void setDestination(unsigned short destination);
    void setPriority(bool priority);

    void setHighPriority();
    void setLowPriority();

    bool operator<(const Packet &rhs);
    bool operator==(const Packet &rhs);

    bool findAndRemove(unsigned short destination); // true or false if found. if found then removes
};


#endif //SIMULATOR_PACKET_HPP
