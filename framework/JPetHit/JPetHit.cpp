#include "./JPetHit.h"

Vector3::Vector3() : 
    fCords{0.0f,0.0f,0.0f} {
}

Vector3::Vector3(float x, float y, float z) :
    fCords{x,y,z} {
}

Vector3::~Vector3() {
}

JPetSignals::JPetSignals() :
    fLeft(NULL), fRight(NULL) {
}

JPetSignals::JPetSignals (JPetSignal& leftSignal, JPetSignal& rightSignal) :
    fLeft(&leftSignal), fRight(&rightSignal) {
}

JPetSignals::~JPetSignals() {
}

JPetHit::JPetHit() :
    fEnergy(0.0f), fQualityOfEnergy(0.0f), fTime(0.0f), fQualityOfTime(0.0f),
    fPos(NULL), fSignals(NULL), fBarrelSlot(NULL), fScintillator(NULL), fTSlot(NULL) {
}

JPetHit::JPetHit(float e, float qe, float t, float qt, Vector3& pos, JPetSignals& sigs, 
                  JPetBarrelSlot& bslot, JPetScin& scin, JPetTSlot& tslot) : 
    fEnergy(e), fQualityOfEnergy(qe), fTime(t), fQualityOfTime(qt),
    fPos(&pos), fSignals(&sigs), fBarrelSlot(&bslot), fScintillator(&scin), fTSlot(&tslot) {
}

JPetHit::~JPetHit() {
}
