#include "Packet.hpp"

Packet::Packet() {
    this->source = 0;
    this->highPriority = false;
}

Packet::Packet(unsigned int uniqueID, unsigned short source, std::set<unsigned short> destination, unsigned int creationTick,
               bool highPriority)
        : uniqueID{uniqueID}
        , source{source}
        , destination{destination}
        , creationTick{creationTick}
        , highPriority{highPriority}
    {}

unsigned short Packet::getSource() { return this->source; }

const std::set<unsigned short>& Packet::getDestination() const { return this->destination; }

bool Packet::getPriority() const { return this->highPriority; }

unsigned int Packet::getUniqueID() { return this->uniqueID; }

unsigned int Packet::getCreationTick() const { return this->creationTick; }

bool Packet::isHighPriority() const { return this->highPriority; }

bool Packet::isLowPriority() { return !this->highPriority; }

void Packet::setUniqueID(unsigned int uniqueID) { this->uniqueID = uniqueID; }

void Packet::setSource(unsigned short source) { this->source = source; }

void Packet::setDestination(std::set<unsigned short> destination) { this->destination = destination; }

void Packet::setDestination(unsigned short destination) { this->destination = {destination}; }

void Packet::setPriority(bool highPriority) { this->highPriority = highPriority; }

void Packet::setHighPriority() { this->highPriority = true; }

void Packet::setLowPriority() { this->highPriority = false; }

bool Packet::findAndRemove(unsigned short destination) {
    auto temp = this->destination.find(destination);

    if( temp == this->destination.end()) {
        return false;
    }
    else {
        this->destination.erase(temp);
        return true;
    }
}

bool Packet::operator<(const Packet &rhs) {
    return this->source < rhs.source;
}

bool Packet::operator==(const Packet &rhs) {
    return this->source == rhs.source;
}
