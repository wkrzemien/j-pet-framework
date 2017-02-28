//
// File generated by /home/klara/root-system/bin/rootcint at Fri Feb 17 19:33:32 2017

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME HitdOh_Dict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "Hit.h_Dict.h"

#include "TClass.h"
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

// Direct notice to TROOT of the dictionary's loading.
namespace {
   static struct DictInit {
      DictInit() {
         ROOT::RegisterModule();
      }
   } __TheDictionaryInitializer;
}

// START OF SHADOWS

namespace ROOTShadow {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOTShadow
// END OF SHADOWS

namespace ROOTDict {
   void Hit_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_Hit(void *p = 0);
   static void *newArray_Hit(Long_t size, void *p);
   static void delete_Hit(void *p);
   static void deleteArray_Hit(void *p);
   static void destruct_Hit(void *p);
   static void streamer_Hit(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static ROOT::TGenericClassInfo *GenerateInitInstanceLocal(const ::Hit*)
   {
      ::Hit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Hit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Hit", ::Hit::Class_Version(), "/home/klara/Studia/III_rok/licencjat/j-pet-framework/Unpacker2/Unpacker2/Hit.h", 7,
                  typeid(::Hit), ::ROOT::DefineBehavior(ptr, ptr),
                  &::Hit::Dictionary, isa_proxy, 0,
                  sizeof(::Hit) );
      instance.SetNew(&new_Hit);
      instance.SetNewArray(&newArray_Hit);
      instance.SetDelete(&delete_Hit);
      instance.SetDeleteArray(&deleteArray_Hit);
      instance.SetDestructor(&destruct_Hit);
      instance.SetStreamerFunc(&streamer_Hit);
      return &instance;
   }
   ROOT::TGenericClassInfo *GenerateInitInstance(const ::Hit*)
   {
      return GenerateInitInstanceLocal((::Hit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::Hit*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOTDict

//______________________________________________________________________________
atomic_TClass_ptr Hit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Hit::Class_Name()
{
   return "Hit";
}

//______________________________________________________________________________
const char *Hit::ImplFileName()
{
   return ::ROOTDict::GenerateInitInstanceLocal((const ::Hit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Hit::ImplFileLine()
{
   return ::ROOTDict::GenerateInitInstanceLocal((const ::Hit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void Hit::Dictionary()
{
   fgIsA = ::ROOTDict::GenerateInitInstanceLocal((const ::Hit*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *Hit::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gCINTMutex); if(!fgIsA) {fgIsA = ::ROOTDict::GenerateInitInstanceLocal((const ::Hit*)0x0)->GetClass();} }
   return fgIsA;
}

//______________________________________________________________________________
void Hit::Streamer(TBuffer &R__b)
{
   // Stream an object of class Hit.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, Hit::IsA());
   } else {
      R__c = R__b.WriteVersion(Hit::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

//______________________________________________________________________________
void Hit::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class Hit.
      TClass *R__cl = ::Hit::IsA();
      if (R__cl || R__insp.IsA()) { }
      TObject::ShowMembers(R__insp);
}

namespace ROOTDict {
   // Wrappers around operator new
   static void *new_Hit(void *p) {
      return  p ? new(p) ::Hit : new ::Hit;
   }
   static void *newArray_Hit(Long_t nElements, void *p) {
      return p ? new(p) ::Hit[nElements] : new ::Hit[nElements];
   }
   // Wrapper around operator delete
   static void delete_Hit(void *p) {
      delete ((::Hit*)p);
   }
   static void deleteArray_Hit(void *p) {
      delete [] ((::Hit*)p);
   }
   static void destruct_Hit(void *p) {
      typedef ::Hit current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_Hit(TBuffer &buf, void *obj) {
      ((::Hit*)obj)->::Hit::Streamer(buf);
   }
} // end of namespace ROOTDict for class ::Hit

/********************************************************
* Hit.h_Dict.cxx
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

extern "C" void G__cpp_reset_tagtableHitdOh_Dict();

extern "C" void G__set_cpp_environmentHitdOh_Dict() {
  G__cpp_reset_tagtableHitdOh_Dict();
}
#include <new>
extern "C" int G__cpp_dllrevHitdOh_Dict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* Hit */
static int G__HitdOh_Dict_168_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   Hit* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new Hit[n];
     } else {
       p = new((void*) gvp) Hit[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new Hit;
     } else {
       p = new((void*) gvp) Hit;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__HitdOh_DictLN_Hit));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HitdOh_Dict_168_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) Hit::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HitdOh_Dict_168_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) Hit::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HitdOh_Dict_168_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) Hit::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HitdOh_Dict_168_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      Hit::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HitdOh_Dict_168_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((Hit*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HitdOh_Dict_168_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) Hit::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HitdOh_Dict_168_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) Hit::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HitdOh_Dict_168_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) Hit::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HitdOh_Dict_168_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) Hit::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__HitdOh_Dict_168_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   Hit* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new Hit(*(Hit*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__HitdOh_DictLN_Hit));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef Hit G__THit;
static int G__HitdOh_Dict_168_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (Hit*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((Hit*) (soff+(sizeof(Hit)*i)))->~G__THit();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (Hit*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((Hit*) (soff))->~G__THit();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__HitdOh_Dict_168_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   Hit* dest = (Hit*) G__getstructoffset();
   *dest = *(Hit*) libp->para[0].ref;
   const Hit& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* Hit */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncHitdOh_Dict {
 public:
  G__Sizep2memfuncHitdOh_Dict(): p(&G__Sizep2memfuncHitdOh_Dict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncHitdOh_Dict::*p)();
};

size_t G__get_sizep2memfuncHitdOh_Dict()
{
  G__Sizep2memfuncHitdOh_Dict a;
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
extern "C" void G__cpp_setup_inheritanceHitdOh_Dict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__HitdOh_DictLN_Hit))) {
     Hit *G__Lderived;
     G__Lderived=(Hit*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__HitdOh_DictLN_Hit),G__get_linked_tagnum(&G__HitdOh_DictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableHitdOh_Dict() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__HitdOh_DictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__HitdOh_DictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__HitdOh_DictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__HitdOh_DictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__HitdOh_DictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__HitdOh_DictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__HitdOh_DictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__HitdOh_DictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__HitdOh_DictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__HitdOh_DictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* Hit */
static void G__setup_memvarHit(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__HitdOh_DictLN_Hit));
   { Hit *p; p=(Hit*)0x1000; if (p) { }
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__HitdOh_DictLN_TClass),G__defined_typename("atomic_TClass_ptr"),-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarHitdOh_Dict() {
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
static void G__setup_memfuncHit(void) {
   /* Hit */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__HitdOh_DictLN_Hit));
   G__memfunc_setup("Hit",293,G__HitdOh_Dict_168_0_1, 105, G__get_linked_tagnum(&G__HitdOh_DictLN_Hit), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__HitdOh_Dict_168_0_2, 85, G__get_linked_tagnum(&G__HitdOh_DictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&Hit::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__HitdOh_Dict_168_0_3, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&Hit::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__HitdOh_Dict_168_0_4, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&Hit::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__HitdOh_Dict_168_0_5, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&Hit::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__HitdOh_DictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__HitdOh_Dict_168_0_9, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - ClassDef_StreamerNVirtual_b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__HitdOh_Dict_168_0_10, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&Hit::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__HitdOh_Dict_168_0_11, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&Hit::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__HitdOh_Dict_168_0_12, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&Hit::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__HitdOh_Dict_168_0_13, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&Hit::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("Hit", 293, G__HitdOh_Dict_168_0_14, (int) ('i'), G__get_linked_tagnum(&G__HitdOh_DictLN_Hit), -1, 0, 1, 1, 1, 0, "u 'Hit' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~Hit", 419, G__HitdOh_Dict_168_0_15, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__HitdOh_Dict_168_0_16, (int) ('u'), G__get_linked_tagnum(&G__HitdOh_DictLN_Hit), -1, 1, 1, 1, 1, 0, "u 'Hit' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncHitdOh_Dict() {
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
extern "C" void G__cpp_setup_globalHitdOh_Dict() {
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

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcHitdOh_Dict() {
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
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__HitdOh_DictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__HitdOh_DictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__HitdOh_DictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__HitdOh_DictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__HitdOh_DictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__HitdOh_DictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__HitdOh_DictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__HitdOh_DictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__HitdOh_DictLN_Hit = { "Hit" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableHitdOh_Dict() {
  G__HitdOh_DictLN_TClass.tagnum = -1 ;
  G__HitdOh_DictLN_TBuffer.tagnum = -1 ;
  G__HitdOh_DictLN_TMemberInspector.tagnum = -1 ;
  G__HitdOh_DictLN_TObject.tagnum = -1 ;
  G__HitdOh_DictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__HitdOh_DictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__HitdOh_DictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__HitdOh_DictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__HitdOh_DictLN_Hit.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableHitdOh_Dict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__HitdOh_DictLN_TClass);
   G__get_linked_tagnum_fwd(&G__HitdOh_DictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__HitdOh_DictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__HitdOh_DictLN_TObject);
   G__get_linked_tagnum_fwd(&G__HitdOh_DictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__HitdOh_DictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__HitdOh_DictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__HitdOh_DictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__HitdOh_DictLN_Hit),sizeof(Hit),-1,29952,(char*)NULL,G__setup_memvarHit,G__setup_memfuncHit);
}
extern "C" void G__cpp_setupHitdOh_Dict(void) {
  G__check_setup_version(30051515,"G__cpp_setupHitdOh_Dict()");
  G__set_cpp_environmentHitdOh_Dict();
  G__cpp_setup_tagtableHitdOh_Dict();

  G__cpp_setup_inheritanceHitdOh_Dict();

  G__cpp_setup_typetableHitdOh_Dict();

  G__cpp_setup_memvarHitdOh_Dict();

  G__cpp_setup_memfuncHitdOh_Dict();
  G__cpp_setup_globalHitdOh_Dict();
  G__cpp_setup_funcHitdOh_Dict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncHitdOh_Dict();
  return;
}
class G__cpp_setup_initHitdOh_Dict {
  public:
    G__cpp_setup_initHitdOh_Dict() { G__add_setup_func("HitdOh_Dict",(G__incsetup)(&G__cpp_setupHitdOh_Dict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initHitdOh_Dict() { G__remove_setup_func("HitdOh_Dict"); }
};
G__cpp_setup_initHitdOh_Dict G__cpp_setup_initializerHitdOh_Dict;

