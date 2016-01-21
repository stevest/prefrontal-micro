/* Created by Language version: 6.2.0 */
/* VECTORIZED */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
 
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gnafbar _p[0]
#define ar2 _p[1]
#define ina _p[2]
#define gna _p[3]
#define m _p[4]
#define h _p[5]
#define s _p[6]
#define Dm _p[7]
#define Dh _p[8]
#define Ds _p[9]
#define NAFX_minf _p[10]
#define NAFX_hinf _p[11]
#define NAFX_sinf _p[12]
#define NAFX_mtau _p[13]
#define NAFX_htau _p[14]
#define NAFX_stau _p[15]
#define ena _p[16]
#define v _p[17]
#define _g _p[18]
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_NAFX_betr(void);
 static void _hoc_NAFX_alpr(void);
 static void _hoc_NAFX_alpv(void);
 static void _hoc_NAFX_hbet(void);
 static void _hoc_NAFX_half(void);
 static void _hoc_NAFX_mbet(void);
 static void _hoc_NAFX_malf(void);
 static void _hoc_rate(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_Nafx", _hoc_setdata,
 "NAFX_betr_Nafx", _hoc_NAFX_betr,
 "NAFX_alpr_Nafx", _hoc_NAFX_alpr,
 "NAFX_alpv_Nafx", _hoc_NAFX_alpv,
 "NAFX_hbet_Nafx", _hoc_NAFX_hbet,
 "NAFX_half_Nafx", _hoc_NAFX_half,
 "NAFX_mbet_Nafx", _hoc_NAFX_mbet,
 "NAFX_malf_Nafx", _hoc_NAFX_malf,
 "rate_Nafx", _hoc_rate,
 0, 0
};
#define NAFX_betr NAFX_betr_Nafx
#define NAFX_alpr NAFX_alpr_Nafx
#define NAFX_alpv NAFX_alpv_Nafx
#define NAFX_hbet NAFX_hbet_Nafx
#define NAFX_half NAFX_half_Nafx
#define NAFX_mbet NAFX_mbet_Nafx
#define NAFX_malf NAFX_malf_Nafx
 extern double NAFX_betr( _threadargsprotocomma_ double );
 extern double NAFX_alpr( _threadargsprotocomma_ double );
 extern double NAFX_alpv( _threadargsprotocomma_ double );
 extern double NAFX_hbet( _threadargsprotocomma_ double );
 extern double NAFX_half( _threadargsprotocomma_ double );
 extern double NAFX_mbet( _threadargsprotocomma_ double );
 extern double NAFX_malf( _threadargsprotocomma_ double );
 /* declare global and static user variables */
#define NAFX_gmr NAFX_gmr_Nafx
 double NAFX_gmr = 0.2;
#define NAFX_zetas NAFX_zetas_Nafx
 double NAFX_zetas = 12;
#define NAFX_zetar NAFX_zetar_Nafx
 double NAFX_zetar = 12;
#define NAFX_b0r NAFX_b0r_Nafx
 double NAFX_b0r = 0.0003;
#define NAFX_a0r NAFX_a0r_Nafx
 double NAFX_a0r = 0.0003;
#define NAFX_vvs NAFX_vvs_Nafx
 double NAFX_vvs = 2;
#define NAFX_vvh NAFX_vvh_Nafx
 double NAFX_vvh = -58;
#define NAFX_vhalfr NAFX_vhalfr_Nafx
 double NAFX_vhalfr = -60;
#define NAFX_taumin NAFX_taumin_Nafx
 double NAFX_taumin = 30;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "NAFX_taumin_Nafx", "ms",
 "NAFX_vhalfr_Nafx", "mV",
 "NAFX_vvh_Nafx", "mV",
 "NAFX_vvs_Nafx", "mV",
 "NAFX_a0r_Nafx", "/ms",
 "NAFX_b0r_Nafx", "/ms",
 "gnafbar_Nafx", "mho/cm2",
 "ina_Nafx", "mA/cm2",
 "gna_Nafx", "mho/cm2",
 0,0
};
 static double delta_t = 1;
 static double h0 = 0;
 static double m0 = 0;
 static double s0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "NAFX_taumin_Nafx", &NAFX_taumin_Nafx,
 "NAFX_vhalfr_Nafx", &NAFX_vhalfr_Nafx,
 "NAFX_vvh_Nafx", &NAFX_vvh_Nafx,
 "NAFX_vvs_Nafx", &NAFX_vvs_Nafx,
 "NAFX_a0r_Nafx", &NAFX_a0r_Nafx,
 "NAFX_b0r_Nafx", &NAFX_b0r_Nafx,
 "NAFX_zetar_Nafx", &NAFX_zetar_Nafx,
 "NAFX_zetas_Nafx", &NAFX_zetas_Nafx,
 "NAFX_gmr_Nafx", &NAFX_gmr_Nafx,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[3]._i
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "6.2.0",
"Nafx",
 "gnafbar_Nafx",
 "ar2_Nafx",
 0,
 "ina_Nafx",
 "gna_Nafx",
 0,
 "m_Nafx",
 "h_Nafx",
 "s_Nafx",
 0,
 0};
 static Symbol* _na_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 19, _prop);
 	/*initialize range parameters*/
 	gnafbar = 0;
 	ar2 = 1;
 	_prop->param = _p;
 	_prop->param_size = 19;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*f)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _nafx_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("na", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_dparam_size(_mechtype, 4);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 Nafx /cluster/stefanos/Documents/GitHub/prefrontal-micro/mechanism_simple/x86_64/nafx.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rate(_threadargsprotocomma_ double, double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[3], _dlist1[3];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   rate ( _threadargscomma_ v , ar2 ) ;
   Dm = ( NAFX_minf - m ) / NAFX_mtau ;
   Dh = ( NAFX_hinf - h ) / NAFX_htau ;
   Ds = ( NAFX_sinf - s ) / NAFX_stau ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 rate ( _threadargscomma_ v , ar2 ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / NAFX_mtau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / NAFX_htau )) ;
 Ds = Ds  / (1. - dt*( ( ( ( - 1.0 ) ) ) / NAFX_stau )) ;
 return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   rate ( _threadargscomma_ v , ar2 ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / NAFX_mtau)))*(- ( ( ( NAFX_minf ) ) / NAFX_mtau ) / ( ( ( ( - 1.0) ) ) / NAFX_mtau ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / NAFX_htau)))*(- ( ( ( NAFX_hinf ) ) / NAFX_htau ) / ( ( ( ( - 1.0) ) ) / NAFX_htau ) - h) ;
    s = s + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / NAFX_stau)))*(- ( ( ( NAFX_sinf ) ) / NAFX_stau ) / ( ( ( ( - 1.0) ) ) / NAFX_stau ) - s) ;
   }
  return 0;
}
 
double NAFX_malf ( _threadargsprotocomma_ double _lv ) {
   double _lNAFX_malf;
 double _lva ;
 _lva = _lv + 28.0 ;
   if ( fabs ( _lva ) < 1e-04 ) {
     _lNAFX_malf = - 0.2816 * ( - 9.3 + _lva * 0.5 ) ;
     }
   else {
     _lNAFX_malf = - 0.2816 * ( _lv + 28.0 ) / ( - 1.0 + exp ( - _lva / 9.3 ) ) ;
     }
   
return _lNAFX_malf;
 }
 
static void _hoc_NAFX_malf(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  NAFX_malf ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double NAFX_mbet ( _threadargsprotocomma_ double _lv ) {
   double _lNAFX_mbet;
 double _lvb ;
 _lvb = _lv + 1.0 ;
   if ( fabs ( _lvb ) < 1e-04 ) {
     _lNAFX_mbet = 0.2464 * ( 6.0 + _lvb * 0.5 ) ;
     }
   else {
     _lNAFX_mbet = 0.2464 * ( _lv + 1.0 ) / ( - 1.0 + exp ( _lvb / 6.0 ) ) ;
     }
   
return _lNAFX_mbet;
 }
 
static void _hoc_NAFX_mbet(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  NAFX_mbet ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double NAFX_half ( _threadargsprotocomma_ double _lv ) {
   double _lNAFX_half;
 double _lvc ;
 _lvc = _lv + 40.1 ;
   if ( fabs ( _lvc ) < 1e-04 ) {
     _lNAFX_half = 0.098 * ( 20.0 + _lvc * 0.5 ) ;
     }
   else {
     _lNAFX_half = 0.098 / exp ( _lvc + 43.1 / 20.0 ) ;
     }
   
return _lNAFX_half;
 }
 
static void _hoc_NAFX_half(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  NAFX_half ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double NAFX_hbet ( _threadargsprotocomma_ double _lv ) {
   double _lNAFX_hbet;
 double _lvd ;
 _lvd = _lv + 13.1 ;
   if ( fabs ( _lvd ) < 1e-04 ) {
     _lNAFX_hbet = 1.4 * ( 10.0 + _lvd * 0.5 ) ;
     }
   else {
     _lNAFX_hbet = 1.4 / ( 1.0 + exp ( - ( _lvd - 13.1 ) / 10.0 ) ) ;
     }
   
return _lNAFX_hbet;
 }
 
static void _hoc_NAFX_hbet(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  NAFX_hbet ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double NAFX_alpv ( _threadargsprotocomma_ double _lv ) {
   double _lNAFX_alpv;
 _lNAFX_alpv = 1.0 / ( 1.0 + exp ( ( _lv - NAFX_vvh ) / NAFX_vvs ) ) ;
   
return _lNAFX_alpv;
 }
 
static void _hoc_NAFX_alpv(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  NAFX_alpv ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double NAFX_alpr ( _threadargsprotocomma_ double _lv ) {
   double _lNAFX_alpr;
 _lNAFX_alpr = exp ( 1.e-3 * NAFX_zetar * ( _lv - NAFX_vhalfr ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lNAFX_alpr;
 }
 
static void _hoc_NAFX_alpr(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  NAFX_alpr ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double NAFX_betr ( _threadargsprotocomma_ double _lv ) {
   double _lNAFX_betr;
 _lNAFX_betr = exp ( 1.e-3 * NAFX_zetar * NAFX_gmr * ( _lv - NAFX_vhalfr ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lNAFX_betr;
 }
 
static void _hoc_NAFX_betr(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  NAFX_betr ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int  rate ( _threadargsprotocomma_ double _lv , double _lNAFX_a2 ) {
   double _lNAFX_q10 , _lNAFX_msum , _lNAFX_hsum , _lNAFX_ma , _lNAFX_mb , _lNAFX_ha , _lNAFX_hb , _lNAFX_c ;
 _lNAFX_ma = NAFX_malf ( _threadargscomma_ _lv ) ;
   _lNAFX_mb = NAFX_mbet ( _threadargscomma_ _lv ) ;
   _lNAFX_ha = NAFX_half ( _threadargscomma_ _lv ) ;
   _lNAFX_hb = NAFX_hbet ( _threadargscomma_ _lv ) ;
   _lNAFX_msum = _lNAFX_ma + _lNAFX_mb ;
   NAFX_minf = _lNAFX_ma / _lNAFX_msum ;
   NAFX_mtau = 1.0 / ( _lNAFX_msum ) ;
   _lNAFX_hsum = _lNAFX_ha + _lNAFX_hb ;
   NAFX_hinf = _lNAFX_ha / _lNAFX_hsum ;
   NAFX_htau = 1.0 / ( _lNAFX_hsum ) ;
   NAFX_stau = NAFX_betr ( _threadargscomma_ _lv ) / ( NAFX_a0r * ( 1.0 + NAFX_alpr ( _threadargscomma_ _lv ) ) ) ;
   if ( NAFX_stau < NAFX_taumin ) {
     NAFX_stau = NAFX_taumin ;
     }
   _lNAFX_c = NAFX_alpv ( _threadargscomma_ _lv ) ;
   NAFX_sinf = _lNAFX_c + _lNAFX_a2 * ( 1.0 - _lNAFX_c ) ;
    return 0; }
 
static void _hoc_rate(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rate ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 3;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  h = h0;
  m = m0;
  s = s0;
 {
   rate ( _threadargscomma_ v , ar2 ) ;
   m = NAFX_minf ;
   h = NAFX_hinf ;
   s = NAFX_sinf ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ena = _ion_ena;
 initmodel(_p, _ppvar, _thread, _nt);
 }}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   gna = gnafbar * m * m * m * h * s ;
   ina = gna * ( v - 55.0 ) ;
   }
 _current += ina;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ena = _ion_ena;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
 double _break, _save;
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _break = t + .5*dt; _save = t;
 v=_v;
{
  ena = _ion_ena;
 { {
 for (; t < _break; t += dt) {
   states(_p, _ppvar, _thread, _nt);
  
}}
 t = _save;
 } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
 _slist1[2] = &(s) - _p;  _dlist1[2] = &(Ds) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif
