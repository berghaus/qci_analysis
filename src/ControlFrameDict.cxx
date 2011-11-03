//
// File generated by rootcint at Thu Nov  3 14:11:45 2011

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME srcdIControlFrameDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "ControlFrameDict.hpp"

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
   void ControlFrame_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void delete_ControlFrame(void *p);
   static void deleteArray_ControlFrame(void *p);
   static void destruct_ControlFrame(void *p);
   static void streamer_ControlFrame(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ControlFrame*)
   {
      ::ControlFrame *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ControlFrame >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ControlFrame", ::ControlFrame::Class_Version(), "./inc/ControlFrame.hpp", 8,
                  typeid(::ControlFrame), DefineBehavior(ptr, ptr),
                  &::ControlFrame::Dictionary, isa_proxy, 0,
                  sizeof(::ControlFrame) );
      instance.SetDelete(&delete_ControlFrame);
      instance.SetDeleteArray(&deleteArray_ControlFrame);
      instance.SetDestructor(&destruct_ControlFrame);
      instance.SetStreamerFunc(&streamer_ControlFrame);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ControlFrame*)
   {
      return GenerateInitInstanceLocal((::ControlFrame*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::ControlFrame*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *ControlFrame::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *ControlFrame::Class_Name()
{
   return "ControlFrame";
}

//______________________________________________________________________________
const char *ControlFrame::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ControlFrame*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ControlFrame::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ControlFrame*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void ControlFrame::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ControlFrame*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *ControlFrame::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ControlFrame*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void ControlFrame::Streamer(TBuffer &R__b)
{
   // Stream an object of class ControlFrame.

   TGMainFrame::Streamer(R__b);
}

//______________________________________________________________________________
void ControlFrame::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class ControlFrame.
      TClass *R__cl = ::ControlFrame::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*fCframe", &fCframe);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*fStart", &fStart);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*fPause", &fPause);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*fExit", &fExit);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "start", &start);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "pause", &pause);
      TGMainFrame::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_ControlFrame(void *p) {
      delete ((::ControlFrame*)p);
   }
   static void deleteArray_ControlFrame(void *p) {
      delete [] ((::ControlFrame*)p);
   }
   static void destruct_ControlFrame(void *p) {
      typedef ::ControlFrame current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_ControlFrame(TBuffer &buf, void *obj) {
      ((::ControlFrame*)obj)->::ControlFrame::Streamer(buf);
   }
} // end of namespace ROOT for class ::ControlFrame

/********************************************************
* src/ControlFrameDict.cxx
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

extern "C" void G__cpp_reset_tagtableControlFrameDict();

extern "C" void G__set_cpp_environmentControlFrameDict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("inc/ControlFrame.hpp");
  G__cpp_reset_tagtableControlFrameDict();
}
#include <new>
extern "C" int G__cpp_dllrevControlFrameDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* ControlFrame */
static int G__ControlFrameDict_301_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   ControlFrame* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 3
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new ControlFrame(
(TGWindow*) G__int(libp->para[0]), (UInt_t) G__int(libp->para[1])
, (UInt_t) G__int(libp->para[2]));
   } else {
     p = new((void*) gvp) ControlFrame(
(TGWindow*) G__int(libp->para[0]), (UInt_t) G__int(libp->para[1])
, (UInt_t) G__int(libp->para[2]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__ControlFrameDictLN_ControlFrame));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ControlFrameDict_301_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((ControlFrame*) G__getstructoffset())->ChangeStartLabel();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ControlFrameDict_301_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((ControlFrame*) G__getstructoffset())->ChangePauseLabel();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ControlFrameDict_301_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ControlFrame::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ControlFrameDict_301_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) ControlFrame::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ControlFrameDict_301_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) ControlFrame::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ControlFrameDict_301_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ControlFrame::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ControlFrameDict_301_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((ControlFrame*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ControlFrameDict_301_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) ControlFrame::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ControlFrameDict_301_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ControlFrame::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ControlFrameDict_301_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) ControlFrame::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__ControlFrameDict_301_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ControlFrame::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef ControlFrame G__TControlFrame;
static int G__ControlFrameDict_301_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (ControlFrame*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((ControlFrame*) (soff+(sizeof(ControlFrame)*i)))->~G__TControlFrame();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (ControlFrame*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((ControlFrame*) (soff))->~G__TControlFrame();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* ControlFrame */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncControlFrameDict {
 public:
  G__Sizep2memfuncControlFrameDict(): p(&G__Sizep2memfuncControlFrameDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncControlFrameDict::*p)();
};

size_t G__get_sizep2memfuncControlFrameDict()
{
  G__Sizep2memfuncControlFrameDict a;
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
extern "C" void G__cpp_setup_inheritanceControlFrameDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__ControlFrameDictLN_ControlFrame))) {
     ControlFrame *G__Lderived;
     G__Lderived=(ControlFrame*)0x1000;
     {
       TGMainFrame *G__Lpbase=(TGMainFrame*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__ControlFrameDictLN_ControlFrame),G__get_linked_tagnum(&G__ControlFrameDictLN_TGMainFrame),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
     {
       TGCompositeFrame *G__Lpbase=(TGCompositeFrame*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__ControlFrameDictLN_ControlFrame),G__get_linked_tagnum(&G__ControlFrameDictLN_TGCompositeFrame),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
     {
       TGFrame *G__Lpbase=(TGFrame*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__ControlFrameDictLN_ControlFrame),G__get_linked_tagnum(&G__ControlFrameDictLN_TGFrame),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
     {
       TGWindow *G__Lpbase=(TGWindow*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__ControlFrameDictLN_ControlFrame),G__get_linked_tagnum(&G__ControlFrameDictLN_TGWindow),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
     {
       TGObject *G__Lpbase=(TGObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__ControlFrameDictLN_ControlFrame),G__get_linked_tagnum(&G__ControlFrameDictLN_TGObject),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__ControlFrameDictLN_ControlFrame),G__get_linked_tagnum(&G__ControlFrameDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
     {
       TQObject *G__Lpbase=(TQObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__ControlFrameDictLN_ControlFrame),G__get_linked_tagnum(&G__ControlFrameDictLN_TQObject),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableControlFrameDict() {

   /* Setting up typedef entry */
   G__search_typename2("UInt_t",104,-1,0,-1);
   G__setnewtype(-1,"Unsigned integer 4 bytes (unsigned int)",0);
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__ControlFrameDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__ControlFrameDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__ControlFrameDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__ControlFrameDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__ControlFrameDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__ControlFrameDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__ControlFrameDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__ControlFrameDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__ControlFrameDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__ControlFrameDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<std::bidirectional_iterator_tag,TObject*,std::ptrdiff_t,const TObject**,const TObject*&>",117,G__get_linked_tagnum(&G__ControlFrameDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*>",117,G__get_linked_tagnum(&G__ControlFrameDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long>",117,G__get_linked_tagnum(&G__ControlFrameDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long,const TObject**>",117,G__get_linked_tagnum(&G__ControlFrameDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* ControlFrame */
static void G__setup_memvarControlFrame(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__ControlFrameDictLN_ControlFrame));
   { ControlFrame *p; p=(ControlFrame*)0x1000; if (p) { }
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__ControlFrameDictLN_TGCompositeFrame),-1,-1,4,"fCframe=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__ControlFrameDictLN_TGTextButton),-1,-1,4,"fStart=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__ControlFrameDictLN_TGTextButton),-1,-1,4,"fPause=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__ControlFrameDictLN_TGTextButton),-1,-1,4,"fExit=",0,(char*)NULL);
   G__memvar_setup((void*)0,103,0,0,-1,G__defined_typename("Bool_t"),-1,4,"start=",0,(char*)NULL);
   G__memvar_setup((void*)0,103,0,0,-1,G__defined_typename("Bool_t"),-1,4,"pause=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__ControlFrameDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarControlFrameDict() {
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
static void G__setup_memfuncControlFrame(void) {
   /* ControlFrame */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__ControlFrameDictLN_ControlFrame));
   G__memfunc_setup("ControlFrame",1228,G__ControlFrameDict_301_0_1, 105, G__get_linked_tagnum(&G__ControlFrameDictLN_ControlFrame), -1, 0, 3, 1, 1, 0, 
"U 'TGWindow' - 10 - p h - 'UInt_t' 0 - w "
"h - 'UInt_t' 0 - h", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("ChangeStartLabel",1588,G__ControlFrameDict_301_0_2, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("ChangePauseLabel",1572,G__ControlFrameDict_301_0_3, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__ControlFrameDict_301_0_4, 85, G__get_linked_tagnum(&G__ControlFrameDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&ControlFrame::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__ControlFrameDict_301_0_5, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&ControlFrame::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__ControlFrameDict_301_0_6, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&ControlFrame::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__ControlFrameDict_301_0_7, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&ControlFrame::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__ControlFrameDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - insp", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__ControlFrameDict_301_0_11, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__ControlFrameDict_301_0_12, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&ControlFrame::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__ControlFrameDict_301_0_13, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&ControlFrame::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__ControlFrameDict_301_0_14, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&ControlFrame::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__ControlFrameDict_301_0_15, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&ControlFrame::DeclFileLine) ), 0);
   // automatic destructor
   G__memfunc_setup("~ControlFrame", 1354, G__ControlFrameDict_301_0_16, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncControlFrameDict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {
}

static void G__cpp_setup_global2() {
}

static void G__cpp_setup_global3() {
}

static void G__cpp_setup_global4() {
}

static void G__cpp_setup_global5() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalControlFrameDict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
  G__cpp_setup_global2();
  G__cpp_setup_global3();
  G__cpp_setup_global4();
  G__cpp_setup_global5();
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

extern "C" void G__cpp_setup_funcControlFrameDict() {
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
G__linked_taginfo G__ControlFrameDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__ControlFrameDictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__ControlFrameDictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__ControlFrameDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__ControlFrameDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__ControlFrameDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__ControlFrameDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__ControlFrameDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__ControlFrameDictLN_TQObject = { "TQObject" , 99 , -1 };
G__linked_taginfo G__ControlFrameDictLN_TGWindow = { "TGWindow" , 99 , -1 };
G__linked_taginfo G__ControlFrameDictLN_TGObject = { "TGObject" , 99 , -1 };
G__linked_taginfo G__ControlFrameDictLN_TGFrame = { "TGFrame" , 99 , -1 };
G__linked_taginfo G__ControlFrameDictLN_TGCompositeFrame = { "TGCompositeFrame" , 99 , -1 };
G__linked_taginfo G__ControlFrameDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR = { "iterator<bidirectional_iterator_tag,TObject*,long,const TObject**,const TObject*&>" , 115 , -1 };
G__linked_taginfo G__ControlFrameDictLN_TGTextButton = { "TGTextButton" , 99 , -1 };
G__linked_taginfo G__ControlFrameDictLN_TGMainFrame = { "TGMainFrame" , 99 , -1 };
G__linked_taginfo G__ControlFrameDictLN_ControlFrame = { "ControlFrame" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableControlFrameDict() {
  G__ControlFrameDictLN_TClass.tagnum = -1 ;
  G__ControlFrameDictLN_TBuffer.tagnum = -1 ;
  G__ControlFrameDictLN_TMemberInspector.tagnum = -1 ;
  G__ControlFrameDictLN_TObject.tagnum = -1 ;
  G__ControlFrameDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__ControlFrameDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__ControlFrameDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__ControlFrameDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__ControlFrameDictLN_TQObject.tagnum = -1 ;
  G__ControlFrameDictLN_TGWindow.tagnum = -1 ;
  G__ControlFrameDictLN_TGObject.tagnum = -1 ;
  G__ControlFrameDictLN_TGFrame.tagnum = -1 ;
  G__ControlFrameDictLN_TGCompositeFrame.tagnum = -1 ;
  G__ControlFrameDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR.tagnum = -1 ;
  G__ControlFrameDictLN_TGTextButton.tagnum = -1 ;
  G__ControlFrameDictLN_TGMainFrame.tagnum = -1 ;
  G__ControlFrameDictLN_ControlFrame.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableControlFrameDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__ControlFrameDictLN_TClass);
   G__get_linked_tagnum_fwd(&G__ControlFrameDictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__ControlFrameDictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__ControlFrameDictLN_TObject);
   G__get_linked_tagnum_fwd(&G__ControlFrameDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__ControlFrameDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__ControlFrameDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__ControlFrameDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__ControlFrameDictLN_TQObject);
   G__get_linked_tagnum_fwd(&G__ControlFrameDictLN_TGWindow);
   G__get_linked_tagnum_fwd(&G__ControlFrameDictLN_TGObject);
   G__get_linked_tagnum_fwd(&G__ControlFrameDictLN_TGFrame);
   G__get_linked_tagnum_fwd(&G__ControlFrameDictLN_TGCompositeFrame);
   G__get_linked_tagnum_fwd(&G__ControlFrameDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR);
   G__get_linked_tagnum_fwd(&G__ControlFrameDictLN_TGTextButton);
   G__get_linked_tagnum_fwd(&G__ControlFrameDictLN_TGMainFrame);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__ControlFrameDictLN_ControlFrame),sizeof(ControlFrame),-1,62464,(char*)NULL,G__setup_memvarControlFrame,G__setup_memfuncControlFrame);
}
extern "C" void G__cpp_setupControlFrameDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupControlFrameDict()");
  G__set_cpp_environmentControlFrameDict();
  G__cpp_setup_tagtableControlFrameDict();

  G__cpp_setup_inheritanceControlFrameDict();

  G__cpp_setup_typetableControlFrameDict();

  G__cpp_setup_memvarControlFrameDict();

  G__cpp_setup_memfuncControlFrameDict();
  G__cpp_setup_globalControlFrameDict();
  G__cpp_setup_funcControlFrameDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncControlFrameDict();
  return;
}
class G__cpp_setup_initControlFrameDict {
  public:
    G__cpp_setup_initControlFrameDict() { G__add_setup_func("ControlFrameDict",(G__incsetup)(&G__cpp_setupControlFrameDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initControlFrameDict() { G__remove_setup_func("ControlFrameDict"); }
};
G__cpp_setup_initControlFrameDict G__cpp_setup_initializerControlFrameDict;

