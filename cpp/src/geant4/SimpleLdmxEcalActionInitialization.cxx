
#include <SimpleLdmxEcalActionInitialization.h>

SimpleLdmxEcalActionInitialization::SimpleLdmxEcalActionInitialization()
    : G4VUserActionInitialization() {     
}


SimpleLdmxEcalActionInitialization::~SimpleLdmxEcalActionInitialization() {
}

void SimpleLdmxEcalActionInitialization::BuildForMaster() const {
;
}

void SimpleLdmxEcalActionInitialization::Build() const {

    SetUserAction(new SimpleLdmxEcalPrimaryGeneratorAction());
}

