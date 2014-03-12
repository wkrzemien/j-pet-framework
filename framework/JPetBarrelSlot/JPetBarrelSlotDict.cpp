//
// File generated by rootcint at Tue Mar 11 21:21:17 2014

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME JPetBarrelSlotDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "JPetBarrelSlotDict.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"

// START OF SHADOWS

namespace ROOT {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void JPetBarrelSlot_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_JPetBarrelSlot(void *p = 0);
   static void *newArray_JPetBarrelSlot(Long_t size, void *p);
   static void delete_JPetBarrelSlot(void *p);
   static void deleteArray_JPetBarrelSlot(void *p);
   static void destruct_JPetBarrelSlot(void *p);
   static void streamer_JPetBarrelSlot(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::JPetBarrelSlot*)
   {
      ::JPetBarrelSlot *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::JPetBarrelSlot >(0);
      static ::ROOT::TGenericClassInfo 
         instance("JPetBarrelSlot", ::JPetBarrelSlot::Class_Version(), "./JPetBarrelSlot.h", 7,
                  typeid(::JPetBarrelSlot), DefineBehavior(ptr, ptr),
                  &::JPetBarrelSlot::Dictionary, isa_proxy, 0,
                  sizeof(::JPetBarrelSlot) );
      instance.SetNew(&new_JPetBarrelSlot);
      instance.SetNewArray(&newArray_JPetBarrelSlot);
      instance.SetDelete(&delete_JPetBarrelSlot);
      instance.SetDeleteArray(&deleteArray_JPetBarrelSlot);
      instance.SetDestructor(&destruct_JPetBarrelSlot);
      instance.SetStreamerFunc(&streamer_JPetBarrelSlot);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::JPetBarrelSlot*)
   {
      return GenerateInitInstanceLocal((::JPetBarrelSlot*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::JPetBarrelSlot*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *JPetBarrelSlot::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *JPetBarrelSlot::Class_Name()
{
   return "JPetBarrelSlot";
}

//______________________________________________________________________________
const char *JPetBarrelSlot::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::JPetBarrelSlot*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int JPetBarrelSlot::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::JPetBarrelSlot*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void JPetBarrelSlot::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::JPetBarrelSlot*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *JPetBarrelSlot::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::JPetBarrelSlot*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void JPetBarrelSlot::Streamer(TBuffer &R__b)
{
   // Stream an object of class JPetBarrelSlot.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      R__b >> fSlotID;
      R__b >> fLayerID;
      R__b >> fLayerRad;
      R__b >> fSlotTheta;
      R__b.CheckByteCount(R__s, R__c, JPetBarrelSlot::IsA());
   } else {
      R__c = R__b.WriteVersion(JPetBarrelSlot::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      R__b << fSlotID;
      R__b << fLayerID;
      R__b << fLayerRad;
      R__b << fSlotTheta;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

//______________________________________________________________________________
void JPetBarrelSlot::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class JPetBarrelSlot.
      TClass *R__cl = ::JPetBarrelSlot::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fSlotID", &fSlotID);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fLayerID", &fLayerID);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fLayerRad", &fLayerRad);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fSlotTheta", &fSlotTheta);
      TNamed::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_JPetBarrelSlot(void *p) {
      return  p ? new(p) ::JPetBarrelSlot : new ::JPetBarrelSlot;
   }
   static void *newArray_JPetBarrelSlot(Long_t nElements, void *p) {
      return p ? new(p) ::JPetBarrelSlot[nElements] : new ::JPetBarrelSlot[nElements];
   }
   // Wrapper around operator delete
   static void delete_JPetBarrelSlot(void *p) {
      delete ((::JPetBarrelSlot*)p);
   }
   static void deleteArray_JPetBarrelSlot(void *p) {
      delete [] ((::JPetBarrelSlot*)p);
   }
   static void destruct_JPetBarrelSlot(void *p) {
      typedef ::JPetBarrelSlot current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_JPetBarrelSlot(TBuffer &buf, void *obj) {
      ((::JPetBarrelSlot*)obj)->::JPetBarrelSlot::Streamer(buf);
   }
} // end of namespace ROOT for class ::JPetBarrelSlot

/********************************************************
* JPetBarrelSlotDict.cpp
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtableJPetBarrelSlotDict();

extern "C" void G__set_cpp_environmentJPetBarrelSlotDict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("JPetBarrelSlot.h");
  G__cpp_reset_tagtableJPetBarrelSlotDict();
}
#include <new>
extern "C" int G__cpp_dllrevJPetBarrelSlotDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* JPetBarrelSlot */
static int G__JPetBarrelSlotDict_177_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   JPetBarrelSlot* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new JPetBarrelSlot[n];
     } else {
       p = new((void*) gvp) JPetBarrelSlot[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new JPetBarrelSlot;
     } else {
       p = new((void*) gvp) JPetBarrelSlot;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_JPetBarrelSlot));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JPetBarrelSlotDict_177_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   JPetBarrelSlot* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 4
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new JPetBarrelSlot(
(int) G__int(libp->para[0]), (int) G__int(libp->para[1])
, (int) G__int(libp->para[2]), (float) G__double(libp->para[3]));
   } else {
     p = new((void*) gvp) JPetBarrelSlot(
(int) G__int(libp->para[0]), (int) G__int(libp->para[1])
, (int) G__int(libp->para[2]), (float) G__double(libp->para[3]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_JPetBarrelSlot));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JPetBarrelSlotDict_177_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((const JPetBarrelSlot*) G__getstructoffset())->getSlotID());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JPetBarrelSlotDict_177_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((const JPetBarrelSlot*) G__getstructoffset())->getLayerID());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JPetBarrelSlotDict_177_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((const JPetBarrelSlot*) G__getstructoffset())->getLayerRad());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JPetBarrelSlotDict_177_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 102, (double) ((const JPetBarrelSlot*) G__getstructoffset())->getSlotTheta());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JPetBarrelSlotDict_177_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((JPetBarrelSlot*) G__getstructoffset())->setSlotID((int) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JPetBarrelSlotDict_177_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((JPetBarrelSlot*) G__getstructoffset())->setLayerID((int) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JPetBarrelSlotDict_177_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((JPetBarrelSlot*) G__getstructoffset())->setLayerRad((int) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JPetBarrelSlotDict_177_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((JPetBarrelSlot*) G__getstructoffset())->setSlotTheta((float) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JPetBarrelSlotDict_177_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) JPetBarrelSlot::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JPetBarrelSlotDict_177_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) JPetBarrelSlot::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JPetBarrelSlotDict_177_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) JPetBarrelSlot::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JPetBarrelSlotDict_177_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      JPetBarrelSlot::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JPetBarrelSlotDict_177_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((JPetBarrelSlot*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JPetBarrelSlotDict_177_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) JPetBarrelSlot::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JPetBarrelSlotDict_177_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) JPetBarrelSlot::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JPetBarrelSlotDict_177_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) JPetBarrelSlot::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__JPetBarrelSlotDict_177_0_22(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) JPetBarrelSlot::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__JPetBarrelSlotDict_177_0_23(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   JPetBarrelSlot* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new JPetBarrelSlot(*(JPetBarrelSlot*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_JPetBarrelSlot));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef JPetBarrelSlot G__TJPetBarrelSlot;
static int G__JPetBarrelSlotDict_177_0_24(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 1
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (JPetBarrelSlot*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((JPetBarrelSlot*) (soff+(sizeof(JPetBarrelSlot)*i)))->~G__TJPetBarrelSlot();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (JPetBarrelSlot*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((JPetBarrelSlot*) (soff))->~G__TJPetBarrelSlot();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__JPetBarrelSlotDict_177_0_25(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   JPetBarrelSlot* dest = (JPetBarrelSlot*) G__getstructoffset();
   *dest = *(JPetBarrelSlot*) libp->para[0].ref;
   const JPetBarrelSlot& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* JPetBarrelSlot */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncJPetBarrelSlotDict {
 public:
  G__Sizep2memfuncJPetBarrelSlotDict(): p(&G__Sizep2memfuncJPetBarrelSlotDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncJPetBarrelSlotDict::*p)();
};

size_t G__get_sizep2memfuncJPetBarrelSlotDict()
{
  G__Sizep2memfuncJPetBarrelSlotDict a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceJPetBarrelSlotDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_JPetBarrelSlot))) {
     JPetBarrelSlot *G__Lderived;
     G__Lderived=(JPetBarrelSlot*)0x1000;
     {
       TNamed *G__Lpbase=(TNamed*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_JPetBarrelSlot),G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_TNamed),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_JPetBarrelSlot),G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableJPetBarrelSlotDict() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* JPetBarrelSlot */
static void G__setup_memvarJPetBarrelSlot(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_JPetBarrelSlot));
   { JPetBarrelSlot *p; p=(JPetBarrelSlot*)0x1000; if (p) { }
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,4,"fSlotID=",0,(char*)NULL);
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,4,"fLayerID=",0,(char*)NULL);
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,4,"fLayerRad=",0,(char*)NULL);
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"fSlotTheta=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarJPetBarrelSlotDict() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncJPetBarrelSlot(void) {
   /* JPetBarrelSlot */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_JPetBarrelSlot));
   G__memfunc_setup("JPetBarrelSlot",1389,G__JPetBarrelSlotDict_177_0_1, 105, G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_JPetBarrelSlot), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("JPetBarrelSlot",1389,G__JPetBarrelSlotDict_177_0_2, 105, G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_JPetBarrelSlot), -1, 0, 4, 1, 1, 0, 
"i - - 0 - slotId i - - 0 - layerID "
"i - - 0 - layerRad f - - 0 - slotTheta", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getSlotID",879,G__JPetBarrelSlotDict_177_0_3, 105, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getLayerID",970,G__JPetBarrelSlotDict_177_0_4, 105, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getLayerRad",1108,G__JPetBarrelSlotDict_177_0_5, 105, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getSlotTheta",1240,G__JPetBarrelSlotDict_177_0_6, 102, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setSlotID",891,G__JPetBarrelSlotDict_177_0_7, 121, -1, -1, 0, 1, 1, 1, 0, "i - - 0 - id", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setLayerID",982,G__JPetBarrelSlotDict_177_0_8, 121, -1, -1, 0, 1, 1, 1, 0, "i - - 0 - id", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setLayerRad",1120,G__JPetBarrelSlotDict_177_0_9, 121, -1, -1, 0, 1, 1, 1, 0, "i - - 0 - layerRad", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setSlotTheta",1252,G__JPetBarrelSlotDict_177_0_10, 121, -1, -1, 0, 1, 1, 1, 0, "f - - 0 - slotTheta", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__JPetBarrelSlotDict_177_0_11, 85, G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&JPetBarrelSlot::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__JPetBarrelSlotDict_177_0_12, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&JPetBarrelSlot::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__JPetBarrelSlotDict_177_0_13, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&JPetBarrelSlot::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__JPetBarrelSlotDict_177_0_14, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&JPetBarrelSlot::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - insp", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__JPetBarrelSlotDict_177_0_18, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__JPetBarrelSlotDict_177_0_19, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&JPetBarrelSlot::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__JPetBarrelSlotDict_177_0_20, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&JPetBarrelSlot::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__JPetBarrelSlotDict_177_0_21, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&JPetBarrelSlot::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__JPetBarrelSlotDict_177_0_22, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&JPetBarrelSlot::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("JPetBarrelSlot", 1389, G__JPetBarrelSlotDict_177_0_23, (int) ('i'), G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_JPetBarrelSlot), -1, 0, 1, 1, 1, 0, "u 'JPetBarrelSlot' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~JPetBarrelSlot", 1515, G__JPetBarrelSlotDict_177_0_24, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 0);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__JPetBarrelSlotDict_177_0_25, (int) ('u'), G__get_linked_tagnum(&G__JPetBarrelSlotDictLN_JPetBarrelSlot), -1, 1, 1, 1, 1, 0, "u 'JPetBarrelSlot' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncJPetBarrelSlotDict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalJPetBarrelSlotDict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {
}

static void G__cpp_setup_func13() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcJPetBarrelSlotDict() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
  G__cpp_setup_func13();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__JPetBarrelSlotDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__JPetBarrelSlotDictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__JPetBarrelSlotDictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__JPetBarrelSlotDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__JPetBarrelSlotDictLN_TNamed = { "TNamed" , 99 , -1 };
G__linked_taginfo G__JPetBarrelSlotDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__JPetBarrelSlotDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__JPetBarrelSlotDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__JPetBarrelSlotDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__JPetBarrelSlotDictLN_JPetBarrelSlot = { "JPetBarrelSlot" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableJPetBarrelSlotDict() {
  G__JPetBarrelSlotDictLN_TClass.tagnum = -1 ;
  G__JPetBarrelSlotDictLN_TBuffer.tagnum = -1 ;
  G__JPetBarrelSlotDictLN_TMemberInspector.tagnum = -1 ;
  G__JPetBarrelSlotDictLN_TObject.tagnum = -1 ;
  G__JPetBarrelSlotDictLN_TNamed.tagnum = -1 ;
  G__JPetBarrelSlotDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__JPetBarrelSlotDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__JPetBarrelSlotDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__JPetBarrelSlotDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__JPetBarrelSlotDictLN_JPetBarrelSlot.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableJPetBarrelSlotDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__JPetBarrelSlotDictLN_TClass);
   G__get_linked_tagnum_fwd(&G__JPetBarrelSlotDictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__JPetBarrelSlotDictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__JPetBarrelSlotDictLN_TObject);
   G__get_linked_tagnum_fwd(&G__JPetBarrelSlotDictLN_TNamed);
   G__get_linked_tagnum_fwd(&G__JPetBarrelSlotDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__JPetBarrelSlotDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__JPetBarrelSlotDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__JPetBarrelSlotDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__JPetBarrelSlotDictLN_JPetBarrelSlot),sizeof(JPetBarrelSlot),-1,61696,(char*)NULL,G__setup_memvarJPetBarrelSlot,G__setup_memfuncJPetBarrelSlot);
}
extern "C" void G__cpp_setupJPetBarrelSlotDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupJPetBarrelSlotDict()");
  G__set_cpp_environmentJPetBarrelSlotDict();
  G__cpp_setup_tagtableJPetBarrelSlotDict();

  G__cpp_setup_inheritanceJPetBarrelSlotDict();

  G__cpp_setup_typetableJPetBarrelSlotDict();

  G__cpp_setup_memvarJPetBarrelSlotDict();

  G__cpp_setup_memfuncJPetBarrelSlotDict();
  G__cpp_setup_globalJPetBarrelSlotDict();
  G__cpp_setup_funcJPetBarrelSlotDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncJPetBarrelSlotDict();
  return;
}
class G__cpp_setup_initJPetBarrelSlotDict {
  public:
    G__cpp_setup_initJPetBarrelSlotDict() { G__add_setup_func("JPetBarrelSlotDict",(G__incsetup)(&G__cpp_setupJPetBarrelSlotDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initJPetBarrelSlotDict() { G__remove_setup_func("JPetBarrelSlotDict"); }
};
G__cpp_setup_initJPetBarrelSlotDict G__cpp_setup_initializerJPetBarrelSlotDict;

