#ifndef GPDICT_LINKDEF_H
#define GPDICT_LINKDEF_H

#ifdef __ROOTCLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link off all typedefs;

#pragma link C++ class std::vector<UParticle>+;


#pragma link C++ class std::vector<TVector3>+;
#pragma link C++ class std::vector<ROOT::Math::XYZVector>+;


#endif // __ROOTCLING__

#endif // GPDICT_LINKDEF_H