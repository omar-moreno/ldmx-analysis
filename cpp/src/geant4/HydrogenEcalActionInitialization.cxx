
#include <HydrogenEcalActionInitialization.h>

HydrogenEcalActionInitialization::HydrogenEcalActionInitialization()
    : G4VUserActionInitialization() {     
}


HydrogenEcalActionInitialization::~HydrogenEcalActionInitialization() {
}

void HydrogenEcalActionInitialization::BuildForMaster() const {
;
}

void HydrogenEcalActionInitialization::Build() const {

    SetUserAction(new SimpleLdmxEcalPrimaryGeneratorAction());
    SetUserAction(new HydrogenEcalRunAction()); 
    SetUserAction(new HydrogenEcalSteppingAction());
}

