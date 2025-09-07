// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__GlobalPolarization
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "/lhep/users/dflusova/lambda/afterburner/release/build/dict/GPDictTypes.hh"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *vectorlEUParticlegR_Dictionary();
   static void vectorlEUParticlegR_TClassManip(TClass*);
   static void *new_vectorlEUParticlegR(void *p = nullptr);
   static void *newArray_vectorlEUParticlegR(Long_t size, void *p);
   static void delete_vectorlEUParticlegR(void *p);
   static void deleteArray_vectorlEUParticlegR(void *p);
   static void destruct_vectorlEUParticlegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<UParticle>*)
   {
      vector<UParticle> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<UParticle>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<UParticle>", -2, "vector", 425,
                  typeid(vector<UParticle>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEUParticlegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<UParticle>) );
      instance.SetNew(&new_vectorlEUParticlegR);
      instance.SetNewArray(&newArray_vectorlEUParticlegR);
      instance.SetDelete(&delete_vectorlEUParticlegR);
      instance.SetDeleteArray(&deleteArray_vectorlEUParticlegR);
      instance.SetDestructor(&destruct_vectorlEUParticlegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<UParticle> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<UParticle>","std::vector<UParticle, std::allocator<UParticle> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<UParticle>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEUParticlegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<UParticle>*>(nullptr))->GetClass();
      vectorlEUParticlegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEUParticlegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEUParticlegR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<UParticle> : new vector<UParticle>;
   }
   static void *newArray_vectorlEUParticlegR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<UParticle>[nElements] : new vector<UParticle>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEUParticlegR(void *p) {
      delete (static_cast<vector<UParticle>*>(p));
   }
   static void deleteArray_vectorlEUParticlegR(void *p) {
      delete [] (static_cast<vector<UParticle>*>(p));
   }
   static void destruct_vectorlEUParticlegR(void *p) {
      typedef vector<UParticle> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<UParticle>

namespace ROOT {
   static TClass *vectorlETVector3gR_Dictionary();
   static void vectorlETVector3gR_TClassManip(TClass*);
   static void *new_vectorlETVector3gR(void *p = nullptr);
   static void *newArray_vectorlETVector3gR(Long_t size, void *p);
   static void delete_vectorlETVector3gR(void *p);
   static void deleteArray_vectorlETVector3gR(void *p);
   static void destruct_vectorlETVector3gR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TVector3>*)
   {
      vector<TVector3> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TVector3>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TVector3>", -2, "vector", 425,
                  typeid(vector<TVector3>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETVector3gR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TVector3>) );
      instance.SetNew(&new_vectorlETVector3gR);
      instance.SetNewArray(&newArray_vectorlETVector3gR);
      instance.SetDelete(&delete_vectorlETVector3gR);
      instance.SetDeleteArray(&deleteArray_vectorlETVector3gR);
      instance.SetDestructor(&destruct_vectorlETVector3gR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TVector3> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<TVector3>","std::vector<TVector3, std::allocator<TVector3> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<TVector3>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETVector3gR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<TVector3>*>(nullptr))->GetClass();
      vectorlETVector3gR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETVector3gR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETVector3gR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TVector3> : new vector<TVector3>;
   }
   static void *newArray_vectorlETVector3gR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TVector3>[nElements] : new vector<TVector3>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETVector3gR(void *p) {
      delete (static_cast<vector<TVector3>*>(p));
   }
   static void deleteArray_vectorlETVector3gR(void *p) {
      delete [] (static_cast<vector<TVector3>*>(p));
   }
   static void destruct_vectorlETVector3gR(void *p) {
      typedef vector<TVector3> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<TVector3>

namespace ROOT {
   static TClass *vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR_Dictionary();
   static void vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR_TClassManip(TClass*);
   static void *new_vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR(void *p = nullptr);
   static void *newArray_vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR(Long_t size, void *p);
   static void delete_vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR(void *p);
   static void deleteArray_vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR(void *p);
   static void destruct_vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >*)
   {
      vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >", -2, "vector", 425,
                  typeid(vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >) );
      instance.SetNew(&new_vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR);
      instance.SetNewArray(&newArray_vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR);
      instance.SetDelete(&delete_vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR);
      instance.SetDestructor(&destruct_vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >","std::vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>, std::allocator<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >*>(nullptr))->GetClass();
      vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> > : new vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >;
   }
   static void *newArray_vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >[nElements] : new vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR(void *p) {
      delete (static_cast<vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >*>(p));
   }
   static void deleteArray_vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR(void *p) {
      delete [] (static_cast<vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >*>(p));
   }
   static void destruct_vectorlEROOTcLcLMathcLcLDisplacementVector3DlEROOTcLcLMathcLcLCartesian3DlEdoublegRcOROOTcLcLMathcLcLDefaultCoordinateSystemTaggRsPgR(void *p) {
      typedef vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >

namespace {
  void TriggerDictionaryInitialization_libgp_dict_Impl() {
    static const char* headers[] = {
"dict/GPDictTypes.hh",
nullptr
    };
    static const char* includePaths[] = {
"{PDROOT}/include",
"/lhep/users/dflusova/lambda/afterburner/release",
"/lhep/users/dflusova/lambda/afterburner/release/{PDROOT}/include",
"/lhep/users/dflusova/mpdroot/install/include",
"/lhep/users/dflusova/lambda/afterburner/release",
"/lhep/users/dflusova/lambda/afterburner/release/{PDROOT}/include",
"/lhep/users/dflusova/mpdroot/install/include",
"/cvmfs/nica.jinr.ru/sw/202309/slc7_x86-64/ROOT/v6.32.06-1/include",
"/cvmfs/nica.jinr.ru/sw/202309/slc7_x86-64/ROOT/v6.32.06-1/include/",
"/lhep/users/dflusova/lambda/afterburner/release/build/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libgp_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
namespace ROOT{namespace Math{template <class T = double> class __attribute__((annotate("$clingAutoload$Math/GenVector/Cartesian3D.h")))  __attribute__((annotate("$clingAutoload$dict/GPDictTypes.hh")))  Cartesian3D;
}}
namespace ROOT{namespace Math{class __attribute__((annotate("$clingAutoload$Math/GenVector/CoordinateSystemTags.h")))  __attribute__((annotate("$clingAutoload$dict/GPDictTypes.hh")))  DefaultCoordinateSystemTag;}}
namespace ROOT{namespace Math{template <class CoordSystem, class Tag = ROOT::Math::DefaultCoordinateSystemTag> class __attribute__((annotate("$clingAutoload$Math/GenVector/DisplacementVector3D.h")))  __attribute__((annotate("$clingAutoload$dict/GPDictTypes.hh")))  DisplacementVector3D;
}}
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$TVector3.h")))  __attribute__((annotate("$clingAutoload$dict/GPDictTypes.hh")))  TVector3;
class __attribute__((annotate("$clingAutoload$UParticle.h")))  __attribute__((annotate("$clingAutoload$dict/GPDictTypes.hh")))  UParticle;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libgp_dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "dict/GPDictTypes.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libgp_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libgp_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libgp_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libgp_dict() {
  TriggerDictionaryInitialization_libgp_dict_Impl();
}
